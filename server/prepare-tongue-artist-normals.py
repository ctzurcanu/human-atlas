"""One presentation-only normal field; exact existing vertices and indices."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file


def minimum_incident_cosine(normals,faces,geometric):
    values=np.einsum('fki,fi->fk',normals[faces],geometric)
    minimum=np.full(len(normals),np.inf)
    np.minimum.at(minimum,faces.ravel(),values.ravel())
    return minimum,values


def main(base,out):
    began=time.monotonic();out.mkdir(parents=True,exist_ok=True)
    construction_sha=digest_file(base/'construction.json')
    record=json.loads((base/'construction.json').read_text());assert record['numerical_and_resource_passed']
    native=Mesh.load(base.parent/'source-v1/meshes/native-00.npz')
    source_words=native.attributes['normals'];source_field=source_words.astype(float)
    source_field/=np.linalg.norm(source_field,axis=1)[:,None]
    rows=[];proofs=[]
    for row in record['meshes']:
        assert digest_file(row['mesh'])==row['sha256'];mesh=Mesh.load(row['mesh'])
        own=mesh.attributes['normals'].astype(float);own/=np.linalg.norm(own,axis=1)[:,None]
        _,ids,points=igl.point_mesh_squared_distance(mesh.vertices,native.vertices,native.faces)
        triangles=native.triangles[ids];a=triangles[:,1]-triangles[:,0];b=triangles[:,2]-triangles[:,0];c=points-triangles[:,0]
        aa=np.einsum('ij,ij->i',a,a);ab=np.einsum('ij,ij->i',a,b);bb=np.einsum('ij,ij->i',b,b)
        ac=np.einsum('ij,ij->i',a,c);bc=np.einsum('ij,ij->i',b,c);den=aa*bb-ab*ab
        assert np.all(den>0)
        u=(bb*ac-ab*bc)/den;v=(aa*bc-ab*ac)/den;bary=np.c_[1-u-v,u,v]
        sampled=np.einsum('ni,nij->nj',bary,source_field[native.faces[ids]])
        assert np.all(np.linalg.norm(sampled,axis=1)>0)
        sampled/=np.linalg.norm(sampled,axis=1)[:,None]
        cosine=np.einsum('ij,ij->i',own,sampled);eligible=np.abs(cosine)>=.7
        sampled*=np.where(cosine<0,-1.,1.)[:,None]
        tri=mesh.triangles;geometric=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
        geometric/=np.linalg.norm(geometric,axis=1)[:,None]
        own_min,own_corners=minimum_incident_cosine(own,mesh.faces,geometric)
        result=own.copy();weights=np.zeros(len(own))
        remaining=eligible.copy()
        # A finite presentation blend sequence enforces side orientation;
        # no control, triangle or geometric defect is moved or concealed.
        for alpha in [1.,.75,.5,.25]:
            candidate=(1-alpha)*own+alpha*sampled
            candidate/=np.linalg.norm(candidate,axis=1)[:,None]
            lowest,_=minimum_incident_cosine(candidate,mesh.faces,geometric)
            accept=remaining&(lowest>0)
            result[accept]=candidate[accept];weights[accept]=alpha;remaining[accept]=False
        # Area-averaged igl normals can oppose a small face at a sharp crease.
        # A strict geometric normal cone corrects presentation only; geometry
        # and its independent self-intersection checks stay untouched.
        preliminary,corners=minimum_incident_cosine(result,mesh.faces,geometric)
        bad_vertices=np.flatnonzero(preliminary<=0);cone=[]
        for vertex in bad_vertices:
            incident=np.flatnonzero(np.any(mesh.faces==vertex,axis=1))
            lp=linprog([0.,0.,0.,-1.],A_ub=np.c_[-geometric[incident],np.ones(len(incident))],
                b_ub=np.zeros(len(incident)),bounds=[(-1.,1.)]*3+[(0.,1.)],method='highs')
            if not lp.success or lp.x[3]<=1e-10:
                raise ValueError(f'{row["key"]}: vertex {vertex} needs a corner-normal chart')
            result[vertex]=lp.x[:3]/np.linalg.norm(lp.x[:3]);weights[vertex]=0.
            cone.append(dict(vertex=int(vertex),minimum_corner_cosine=float((geometric[incident]@result[vertex]).min())))
        words=np.rint(result*32767).astype('<i2')
        decoded=words.astype(float);decoded/=np.linalg.norm(decoded,axis=1)[:,None]
        packed_min,_=minimum_incident_cosine(decoded,mesh.faces,geometric)
        # Test finite nearest integer words only, without moving geometry or
        # accepting a nonpositive corner. Normally the geometric cone is ample.
        word_adjustments=[]
        for vertex in np.flatnonzero(packed_min<=0):
            incident=np.flatnonzero(np.any(mesh.faces==vertex,axis=1))
            options=[]
            for x in [-1,0,1]:
                for y in [-1,0,1]:
                    for z in [-1,0,1]:
                        w=words[vertex].astype(np.int32)+[x,y,z]
                        if np.any(w<-32767)or np.any(w>32767)or not np.any(w):continue
                        n=w/np.linalg.norm(w);margin=float((geometric[incident]@n).min())
                        if margin>0:options.append((float(np.linalg.norm(n-result[vertex])),w,n,margin))
            if not options:raise ValueError(f'{row["key"]}: no permitted adjacent Int16 normal word at vertex {vertex}')
            _,w,n,margin=min(options,key=lambda t:t[0]);words[vertex]=w;decoded[vertex]=n
            word_adjustments.append(dict(vertex=int(vertex),minimum_corner_cosine=margin))
        result=decoded
        minimum,corners=minimum_incident_cosine(result,mesh.faces,geometric)
        unit_error=float(np.max(np.abs(np.linalg.norm(result,axis=1)-1)))
        changed=np.linalg.norm(result-own,axis=1)>1e-12
        attributes={**mesh.attributes,'normals':words}
        changed_mesh=Mesh(mesh.vertices,mesh.faces,mesh.frame,attributes,{**mesh.source,
            'presentation_normal_variant':'native packed artist-field barycentric closest-face sampling with sign alignment and incident-face side guards',
            'geometry_changed':False})
        assert changed_mesh.geometry_sha==mesh.geometry_sha
        assert np.array_equal(changed_mesh.vertices,mesh.vertices) and np.array_equal(changed_mesh.faces,mesh.faces)
        path=out/(row['key']+'.npz');changed_mesh.save(path)
        p=dict(key=row['key'],base_mesh_sha256=row['sha256'],base_geometry_sha256=mesh.geometry_sha,
            variant_mesh_sha256=digest_file(path),variant_geometry_sha256=changed_mesh.geometry_sha,
            source_tongue_geometry_sha256=native.geometry_sha,source_normal_word_dtype=str(source_words.dtype),
            barycentric_query_count=len(mesh.vertices),minimum_barycentric=float(bary.min()),
            maximum_barycentric=float(bary.max()),eligibility_cosine_threshold=.7,
            eligible_vertices=int(eligible.sum()),changed_vertices=int(changed.sum()),
            full_artist_vertices=int(np.count_nonzero(weights==1.)),
            blended_vertices=int(np.count_nonzero((weights>0)&(weights<1))),
            retained_geometric_vertices=int(np.count_nonzero(weights==0)),
            geometric_cone_corrections=cone,geometric_cone_corrected_vertices=len(cone),
            actual_normal_storage='Int16 words; rendered and side-checked after decoding and unit normalization',
            normal_word_adjustments=word_adjustments,
            blend_counts={str(alpha):int(np.count_nonzero(weights==alpha))for alpha in [0.,.25,.5,.75,1.]},
            original_nonpositive_face_corners=int(np.count_nonzero(own_corners<=0)),
            variant_nonpositive_face_corners=int(np.count_nonzero(corners<=0)),
            minimum_original_corner_cosine=float(own_corners.min()),minimum_variant_corner_cosine=float(corners.min()),
            maximum_unit_length_error=unit_error,finite_nonzero_normals=bool(np.isfinite(result).all()and np.all(np.linalg.norm(result,axis=1)>0)),
            side_orientation_proven=bool(np.all(corners>0)),
            vertices_unchanged=True,indices_unchanged=True,geometry_sha_exact=True,
            maximum_normal_delta=float(np.linalg.norm(result-own,axis=1).max()))
        proofs.append(p)
        rows.append({**row,'mesh':str(path),'sha256':digest_file(path),'geometry_sha256':changed_mesh.geometry_sha})
    passed=all(p['finite_nonzero_normals']and p['side_orientation_proven']and p['maximum_unit_length_error']<1e-12 for p in proofs)
    evidence=dict(schema='tongue-normal-only-presentation-variant/v1',base_construction_sha256=construction_sha,
        source_catalogue_sha256=record['source_catalogue_sha256'],source_tongue_mesh_sha256=digest_file(base.parent/'source-v1/meshes/native-00.npz'),
        geometry_unchanged=True,source_geometry_unchanged=True,styloglossus_bone_caps_untouched=True,
        all_original_source_normal_words_unchanged=True,geometry_defects_repaired_or_hidden=False,
        meshes=proofs,passed=passed,elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'normal-only-proof.json',evidence)
    variant={**record,'meshes':rows,'base_construction_sha256':construction_sha,'presentation_normal_proof_sha256':digest_file(out/'normal-only-proof.json'),
        'presentation_normal_variant_passed':passed,'numerical_and_resource_passed':record['numerical_and_resource_passed']and passed,
        'geometry_checks_reused_only_because_all_geometry_hashes_are_exact':True,'visual_review':'pending',
        'construction_script_sha256':digest_file(__file__)}
    atomic_json(out/'construction.json',variant)
    views=[dict(name='right',direction=[1,0,0]),dict(name='top',direction=[0,1,0])]
    before=dict(frame=record['frame'],resolution=640,wall_seconds=60,views=views,
        meshes=[dict(mesh=r['mesh'],name=r['key'],color=r['color'])for r in record['meshes']])
    after={**before,'meshes':[dict(mesh=r['mesh'],name=r['key'],color=r['color'])for r in rows]}
    atomic_json(out/'review-before.json',before);atomic_json(out/'review-after.json',after)
    assert digest_file(base/'construction.json')==construction_sha
    print(json.dumps(dict(passed=passed,geometry_unchanged=True,
        changed_vertices=sum(p['changed_vertices']for p in proofs),
        nonpositive_corners=sum(p['variant_nonpositive_face_corners']for p in proofs),
        elapsed_seconds=evidence['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--base',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args();main(args.base.resolve(),args.output.resolve())
