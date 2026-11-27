"""Literal-geometry normal charts with complete-face facing constraints."""
from pathlib import Path
from collections import defaultdict
import argparse
import json
import sys
import time
import numpy as np
import igl

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from tongue_whole_input import whole_input,new_output


def guarded_words(desired,geometric):
    """Finite blends; each decoded corner obeys linear face-facing tests."""
    for weight in [1.,.75,.5,.25,0.]:
        n=weight*desired+(1-weight)*geometric;n/=np.linalg.norm(n)
        words=np.rint(n*32767).astype('<i2')
        for axis in [0,1]:
            if abs(geometric[axis])>1e-5 and words[axis]==0:words[axis]=1 if geometric[axis]>0 else-1
        decoded=words.astype(float);decoded/=np.linalg.norm(decoded)
        if decoded@geometric<.7:continue
        if any(abs(geometric[a])>1e-5 and decoded[a]*geometric[a]<=0 for a in [0,1]):continue
        return words,weight
    raise ValueError('No permitted finite face-facing normal word')


def main(base,out,whole=None,fit_result=None,validate_inputs_only=False):
    began=time.monotonic();out=new_output(out)
    base_hash=digest_file(base/'construction.json');record=json.loads((base/'construction.json').read_text())
    prior=record.get('whole_input',{})
    if whole is None and fit_result is None and prior.get('mode')=='explicit-reviewed-corrected-whole':
        whole=Path(prior['whole_mesh']);fit_result=Path(prior['fit_result_path'])
    native,whole_path,binding=whole_input(whole,fit_result)
    if prior:
        if prior['whole_geometry_sha256']!=native.geometry_sha or prior['whole_mesh_sha256']!=digest_file(whole_path):
            raise ValueError('Normal charts must sample the same literal whole used to regenerate these cells')
        if prior.get('fit_result_sha256')!=binding.get('fit_result_sha256'):
            raise ValueError('Normal-chart corrected-whole fit binding differs from construction')
    elif whole is not None:
        raise ValueError('Legacy intrinsic cells cannot be rebound to a corrected whole; regenerate them first')
    if whole is not None and not(record.get('numerical_and_resource_passed')is True and record.get('fixed_bone_crossings_zero')is True):
        raise ValueError('Regenerated cells must pass their actual numerical and fixed-bone checks before normal charts')
    for row in record['meshes']:
        if digest_file(row['mesh'])!=row['sha256']:
            raise ValueError('Base intrinsic NPZ hash changed')
    if validate_inputs_only:
        print(json.dumps(dict(validated=True,generation_started=False,base_construction_sha256=base_hash,
            whole_input=binding,artist_normals_from_selected_whole=True)),flush=True)
        return
    out.mkdir(parents=True)
    artist=native.attributes['normals'].astype(float);artist/=np.linalg.norm(artist,axis=1)[:,None]
    rows=[];proofs=[]
    for row in record['meshes']:
        assert digest_file(row['mesh'])==row['sha256'];mesh=Mesh.load(row['mesh'])
        tri=mesh.triangles;cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1)
        g=cross/area[:,None];incident=defaultdict(set);edges=defaultdict(list);adjacent=defaultdict(set);sharp=[]
        for face,t in enumerate(mesh.faces):
            for vertex in t:incident[int(vertex)].add(face)
            for k in range(3):edges[tuple(sorted((int(t[k]),int(t[(k+1)%3]))))].append(face)
        sharp_corners=set()
        for edge,faces in edges.items():
            assert len(faces)==2
            if g[faces[0]]@g[faces[1]]<np.cos(np.pi/4):
                sharp.append((edge,faces))
                for face in faces:
                    for vertex in edge:sharp_corners.add((face,vertex))
            else:
                adjacent[faces[0]].add(faces[1]);adjacent[faces[1]].add(faces[0])
        _,ids,p=igl.point_mesh_squared_distance(mesh.vertices,native.vertices,native.faces)
        st=native.triangles[ids];e1=st[:,1]-st[:,0];e2=st[:,2]-st[:,0];q=p-st[:,0]
        aa=np.einsum('ij,ij->i',e1,e1);ab=np.einsum('ij,ij->i',e1,e2);bb=np.einsum('ij,ij->i',e2,e2)
        aq=np.einsum('ij,ij->i',e1,q);bq=np.einsum('ij,ij->i',e2,q);den=aa*bb-ab*ab
        u=(bb*aq-ab*bq)/den;v=(aa*bq-ab*aq)/den;bary=np.c_[1-u-v,u,v]
        sampled=np.einsum('ni,nij->nj',bary,artist[native.faces[ids]])
        sampled/=np.linalg.norm(sampled,axis=1)[:,None]
        words=np.zeros((len(mesh.faces),3,3),dtype='<i2');chart_ids=np.zeros((len(mesh.faces),3),dtype=np.int32)
        chart=0;artist_charts=0;corner_blends=0
        for vertex,face_set in incident.items():
            remaining=set(face_set)
            while remaining:
                seed=min(remaining);component={seed};stack=[seed];remaining.remove(seed)
                while stack:
                    current=stack.pop()
                    for other in adjacent[current]&remaining:
                        if other in face_set:component.add(other);stack.append(other);remaining.remove(other)
                fs=sorted(component);geometric=np.sum(cross[fs],axis=0)
                if np.linalg.norm(geometric)==0:geometric=g[fs[0]]
                geometric/=np.linalg.norm(geometric)
                source=sampled[vertex];cosine=source@geometric
                source=source if cosine>=0 else-source
                use_artist=abs(cosine)>=.7;artist_charts+=int(use_artist)
                desired=source if use_artist else geometric
                for face in fs:
                    corner=int(np.flatnonzero(mesh.faces[face]==vertex)[0])
                    # Crease-side geometric charts remain separate; source
                    # artist smoothing is reserved for compatible cap fields.
                    wanted=geometric if (face,vertex)in sharp_corners else desired
                    word,weight=guarded_words(wanted,g[face]);words[face,corner]=word
                    chart_ids[face,corner]=chart;corner_blends+=weight<1.
                chart+=1
        # Preserve a visible normal discontinuity at every >45-degree edge.
        for edge,faces in sharp:
            for vertex in edge:
                slots=[(f,int(np.flatnonzero(mesh.faces[f]==vertex)[0]))for f in faces]
                values=[words[f,k].astype(float)/np.linalg.norm(words[f,k])for f,k in slots]
                if values[0]@values[1]>np.cos(np.deg2rad(5)):
                    for face,k in slots:words[face,k]=guarded_words(g[face],g[face])[0]
        decoded=words.astype(float);decoded/=np.linalg.norm(decoded,axis=2)[:,:,None]
        plane=np.einsum('fki,fi->fk',decoded,g)
        facing=[]
        for axis,name in [(0,'right-left'),(1,'top-bottom')]:
            relevant=np.abs(g[:,axis])>1e-5
            signed=decoded[relevant,:,axis]*np.sign(g[relevant,axis,None])
            facing.append(dict(axis=name,relevant_faces=int(relevant.sum()),grazing_faces=int((~relevant).sum()),
                minimum_corner_facing_component=float(signed.min()),failed_corner_count=int(np.count_nonzero(signed<=0)),
                complete_face_proof='Every barycentric interpolant is a positive weighted sum of corners; unit normalization preserves the signs'))
        discontinuities=[]
        for edge,faces in sharp:
            for vertex in edge:
                corners=[int(np.flatnonzero(mesh.faces[f]==vertex)[0])for f in faces]
                value=float(np.degrees(np.arccos(np.clip(decoded[faces[0],corners[0]]@decoded[faces[1],corners[1]],-1,1))))
                discontinuities.append(value)
        changed=Mesh(mesh.vertices,mesh.faces,mesh.frame,mesh.attributes,
            {**mesh.source,'normal_chart_presentation':'Crease-separated per-corner artist/geometric fields with complete right/top face-facing constraints'})
        assert changed.geometry_sha==mesh.geometry_sha and np.array_equal(changed.vertices,mesh.vertices)and np.array_equal(changed.faces,mesh.faces)
        path=out/(row['key']+'.npz')
        # Mesh attributes remain compact point fields. Explicit corner charts
        # are a separate normal bank and do not alter the position/index bank.
        atomic_npz(path,vertices=changed.vertices,faces=changed.faces,frame=np.array(changed.frame),
            source_json=np.array(json.dumps(changed.source,sort_keys=True)),
            **{'attribute_'+k:v for k,v in changed.attributes.items()},corner_normals=words,normal_chart_ids=chart_ids)
        passed=bool(np.all(plane>=.7)and all(r['failed_corner_count']==0 for r in facing)and np.isfinite(decoded).all()
            and min(discontinuities,default=180)>=4.9)
        proofs.append(dict(key=row['key'],base_mesh_sha256=row['sha256'],variant_mesh_sha256=digest_file(path),
            geometry_sha256=mesh.geometry_sha,vertices_unchanged=True,indices_unchanged=True,
            charts=chart,artist_eligible_charts=artist_charts,sharp_edges=len(sharp),sharp_threshold_degrees=45,
            minimum_sharp_edge_normal_discontinuity_degrees=min(discontinuities,default=180),
            normal_array='corner_normals',word_dtype='int16',word_shape=list(words.shape),
            minimum_complete_face_geometric_dot=float(plane.min()),normalization_nonzero_proof='All corner-plane dots >=0.7, so every barycentric field has positive dot and cannot be zero',
            facing_constraints=facing,guarded_corner_blends=corner_blends,passed=passed))
        rows.append({**row,'mesh':str(path),'sha256':digest_file(path),
            'corner_normal_source':dict(path=str(path),sha256=digest_file(path),array='corner_normals',chart_ids_array='normal_chart_ids')})
    passed=all(r['passed']for r in proofs)
    proof=dict(schema='tongue-literal-geometry-corner-normal-charts/v1',base_construction_path=str(base/'construction.json'),
        base_construction_sha256=base_hash,source_tongue_mesh_sha256=digest_file(whole_path),whole_input=binding,
        meshes=proofs,geometry_unchanged=True,styloglossus_source_words_untouched=True,passed=passed,
        source_bone_crossing_candidate_not_accepted_or_installed=whole is None,
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(out/'normal-chart-proof.json',proof)
    atomic_json(out/'construction.json',{**record,'meshes':rows,'base_construction_sha256':base_hash,
        'normal_chart_proof_sha256':digest_file(out/'normal-chart-proof.json'),'normal_chart_variant_passed':passed,
        'visual_review':'pending','installed':False,'construction_script_sha256':digest_file(__file__)})
    views=[dict(name='right',direction=[1,0,0]),dict(name='top',direction=[0,1,0])]
    atomic_json(out/'review-spec.json',dict(frame=record['frame'],resolution=640,wall_seconds=60,views=views,
        meshes=[dict(mesh=r['mesh'],name=r['key'],color=r['color'],corner_normal_source=r['corner_normal_source'])for r in rows]))
    print(json.dumps(dict(passed=passed,geometry_unchanged=True,meshes=[dict(key=p['key'],charts=p['charts'],
        minimum_plane_dot=p['minimum_complete_face_geometric_dot'])for p in proofs],elapsed_seconds=proof['elapsed_seconds'])),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--base',required=True,type=Path);ap.add_argument('--output',required=True,type=Path)
    ap.add_argument('--whole-tongue',type=Path);ap.add_argument('--fit-result',type=Path)
    ap.add_argument('--validate-inputs-only',action='store_true');args=ap.parse_args()
    main(args.base.resolve(),args.output,args.whole_tongue,args.fit_result,args.validate_inputs_only)
