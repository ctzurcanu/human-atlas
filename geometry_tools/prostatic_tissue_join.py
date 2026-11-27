"""Conform the urinary coats to one prostate bore in their source planes.

The gland and urinary wall share actual Boolean boundary fragments, rather
than superimposed independent contact faces. Original exterior UVs/normals
are interpolated only within their own source triangles.
"""
from collections import defaultdict
from pathlib import Path
import json

import numpy as np

from .algorithms.finite_tissue_boolean import barycentric
from .audit import topology
from .boolean_cells import BooleanCells
from .common import atomic_json, atomic_npz, digest_file
from .model import Mesh
from .spatial import BVH


def key(points):
    return tuple(sorted(tuple(p)for p in points))


def build(arrangement_dir, coordinates_dir, outlet_dir, outlet_sources, output, *, resume=False):
    arrangement_dir,coordinates_dir,outlet_dir,outlet_sources,output=map(lambda p:Path(p).resolve(),
        (arrangement_dir,coordinates_dir,outlet_dir,outlet_sources,output))
    if output.exists() and not resume:
        raise FileExistsError(output)
    spec=json.loads((arrangement_dir/'arrangement.json').read_text())
    coord=json.loads((coordinates_dir/'coordinates.json').read_text())
    wall=Mesh.load(coord['wall']);prostate=Mesh.load(outlet_sources/'prostate.npz')
    if digest_file(coord['wall'])!=spec['source_wall_sha256']:
        raise ValueError('The whole urinary source wall changed')
    with np.load(Path(coord['wall']).parent/'carriers.npz',allow_pickle=False)as z:
        ov,of=z['vertices'].copy(),z['outer_faces'].copy()
    original_outer={key(ov[t]):i for i,t in enumerate(of)}
    wall_to_outer={i:original_outer[key(wall.vertices[t])]for i,t in enumerate(wall.faces)if key(wall.vertices[t])in original_outer}
    with np.load(arrangement_dir/'arrangement.npz',allow_pickle=False)as z:
        v=z['positions'].copy();sf={sid:z[a].copy()for sid,a in spec['surfaces'].items()}
        tags={sid:z[a].copy()for sid,a in spec['source_face_arrays'].items()}
    bindings=dict(source_arrangement_sha256=digest_file(arrangement_dir/'arrangement.npz'),
        source_wall_sha256=coord['wall_sha256'],source_prostate_sha256=digest_file(outlet_sources/'prostate.npz'))
    if resume:
        checkpoint=json.loads((output/'canonical-geometry.json').read_text())
        if any(checkpoint.get(k)!=value for k,value in bindings.items()) or \
                digest_file(output/'canonical-geometry.npz')!=checkpoint['geometry_sha256']:
            raise ValueError('A checkpoint source or geometry changed; refusing stale attribute transfer')
        with np.load(output/'canonical-geometry.npz',allow_pickle=False)as z:
            positions=z['positions'].copy();surfaces={sid:z[key].copy()for sid,key in checkpoint['surfaces'].items()}
        cells=checkpoint['cells'];r=dict(evidence=checkpoint['evidence'])
    else:
        graph=BooleanCells();outer=graph.source(ov,of,'actual urinary outer envelope');outer_oid=max(graph.primitives)
        gland=graph.source(prostate.vertices,prostate.faces,'source prostate');solids={'prostate':gland-outer}
        for name,uses in spec['cells'].items():
            faces=[];ancestry=[]
            for use in uses:
                faces.append(sf[use['surface']]if use['side']==1 else sf[use['surface']][:,[0,2,1]])
                ancestry.extend(map(int,tags[use['surface']]))
            fs=np.concatenate(faces);solid=graph.literal_cell(v,fs,name);oid=max(graph.primitives)
            for i,parent in enumerate(ancestry):
                if parent in wall_to_outer:
                    graph.aliases[(oid,i)]=(outer_oid,wall_to_outer[parent])
            solids[name]=solid
        output.mkdir(parents=True)
        graph.save(solids,output/'boolean-source.npz')
        r=graph.assemble(solids,progress=lambda i,n,o:print('source-plane overlay',i,n,'overlays',o,flush=True))
        positions,surfaces,cells=r['positions'],r['surfaces'],r['cells']
        checkpoint={'positions':positions};checkpoint_keys={}
        for i,(sid,fs)in enumerate(surfaces.items()):
            checkpoint['surface_'+str(i)]=fs;checkpoint_keys[sid]='surface_'+str(i)
        atomic_npz(output/'canonical-geometry.npz',**checkpoint)
        atomic_json(output/'canonical-geometry.json',dict(cells=cells,surfaces=checkpoint_keys,
            frame=wall.frame,evidence=r['evidence'],**bindings,
            geometry_sha256=digest_file(output/'canonical-geometry.npz'),installed=False,anatomical_acceptance=False))
    ownership=defaultdict(list)
    for name,uses in cells.items():
        for use in uses:
            ownership[use['surface']].append(name)
    arrays={'positions':positions};surface_keys={};tag_keys={};uv_keys={};word_keys={}
    wall_index=BVH(wall);prostate_index=BVH(prostate);query_evidence=[]
    source_normal=prostate.attributes['normals'].astype(np.int32)
    if 'normalSigns'in prostate.attributes:
        source_normal*=prostate.attributes['normalSigns'][:,None]
    check_rows=[];root_original={key(wall.vertices[wall.faces[int(i)]]):tuple(role)for i,role in coord['source_root_faces'].items()}
    retained_roots=set()
    for number,(sid,faces)in enumerate(surfaces.items()):
        owners=ownership[sid];q=positions[faces];centers=q.mean(1)
        field=np.full(len(faces),-1,int)
        if any(n!='prostate'for n in owners)and (len(owners)==1 or 'prostate'in owners):
            nearest,evidence=wall_index.closest_points(centers,backend='auto');query_evidence.append(evidence);parents=nearest['face_ids']
            for i,parent in enumerate(parents):
                triangle=wall.vertices[wall.faces[parent]];b=barycentric(q[i],triangle)
                if np.max(np.linalg.norm(b@triangle-q[i],axis=1))>4e-12 or b.min()<-1e-7:
                    raise ValueError('A urinary exterior fragment left its bound original face')
                field[i]=int(parent)
                k=key(q[i])
                if k in root_original:
                    retained_roots.add(k)
        if 'prostate'in owners:
            uv=np.empty((len(faces),3,2),dtype='<f4')
            if len(owners)==1:
                words=np.empty((len(faces),3,3),dtype='<i2')
                nearest,evidence=prostate_index.closest_points(centers,backend='auto');query_evidence.append(evidence);parents=nearest['face_ids']
                for i,parent in enumerate(parents):
                    tri=prostate.faces[parent];b=barycentric(q[i],prostate.vertices[tri])
                    if np.max(np.linalg.norm(b@prostate.vertices[tri]-q[i],axis=1))>4e-12 or b.min()<-1e-7:
                        failure=dict(surface=sid,owners=owners,triangle=q[i].tolist(),source_face=int(parent),
                            source_triangle=prostate.vertices[tri].tolist(),barycentric=b.tolist(),
                            maximum_plane_residual_m=float(np.max(np.linalg.norm(b@prostate.vertices[tri]-q[i],axis=1))),
                            area_m2=float(np.linalg.norm(np.cross(q[i,1]-q[i,0],q[i,2]-q[i,0]))/2),
                            installed=False,anatomical_acceptance=False)
                        atomic_json(output/'uv-lineage-failure.json',failure)
                        raise ValueError('A gland exterior fragment lost its original UV chart; see uv-lineage-failure.json')
                    uv[i]=b@prostate.attributes['uvs'][tri]
                    words[i]=np.rint(b@source_normal[tri]).astype('<i2')
                name='normal_'+str(number);arrays[name]=words;word_keys[sid]=name
            else:
                nearest,evidence=prostate_index.closest_points(q.reshape(-1,3),backend='auto');query_evidence.append(evidence)
                parents,b=nearest['face_ids'],nearest['barycentric']
                uv[:]=np.einsum('ij,ijk->ik',b,prostate.attributes['uvs'][prostate.faces[parents]]).reshape(-1,3,2)
            name='uv_'+str(number);arrays[name]=uv;uv_keys[sid]=name
        name='surface_'+str(number);arrays[name]=faces;surface_keys[sid]=name
        name='source_face_'+str(number);arrays[name]=field;tag_keys[sid]=name
    if retained_roots!=set(root_original):
        raise ValueError('The coupled gland/wall arrangement changed a renal-root facet')
    for name,uses in cells.items():
        fs=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
        check=topology(Mesh(positions,fs,wall.frame))
        if not check['closed_oriented_manifold']:
            Mesh(positions,fs,wall.frame).save(output/'topology-failure.npz')
            atomic_json(output/'topology-failure.json',dict(cell=name,topology=check))
            raise ValueError('The conformed urinary/prostate tissue has invalid topology: '+name)
        check_rows.append(dict(cell=name,topology=check))
    atomic_npz(output/'arrangement.npz',**arrays)
    result=dict(schema='continuous-urinary-wall-arrangement/v1',frame=wall.frame,
        surfaces=surface_keys,source_face_arrays=tag_keys,cells=cells,
        part_ids=dict(spec['part_ids'],prostate='ZA:Prostate'),
        corner_uvs=uv_keys,corner_word_overrides=word_keys,
        root_faces_preserved=len(retained_roots),all_inner_faces_mucosal=spec['all_inner_faces_mucosal'],
        source_wall_sha256=spec['source_wall_sha256'],coordinates_sha256=spec['coordinates_sha256'],
        source_catalogue_sha256=coord['catalogue_sha256'],
        source_arrangement_sha256=digest_file(arrangement_dir/'arrangement.npz'),
        feature_plan=spec.get('feature_plan'),feature_plan_sha256=spec.get('feature_plan_sha256'),
        opening_curves=spec.get('opening_curves'),opening_surround_width_m=spec.get('opening_surround_width_m'),
        opening_support_method=spec.get('opening_support_method'),opening_windows=spec.get('opening_windows'),
        prostate_source_sha256=digest_file(outlet_sources/'prostate.npz'),
        boolean_evidence=r['evidence'],topology=check_rows,
        original_surface_query_evidence=query_evidence,
        source_exterior_uvs_retained=True,gland_wall_interfaces_stored_once=True,
        gland_internal_uv_projection_illustrative=True,anatomical_boundaries_provisional=True,
        installed=False,anatomical_acceptance=False)
    atomic_json(output/'arrangement.json',result)
    return result


if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['arrangement-dir','coordinates-dir','outlet-dir','outlet-sources','output']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--resume',action='store_true')
    r=build(**vars(p.parse_args()));print('conformed',len(r['cells']),'finite urinary/prostate tissues')
