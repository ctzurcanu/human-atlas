#!/usr/bin/env python3
"""One source-bound rigid registration; guides are never installation anatomy."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
from geometry_tools.common import ROOT, atomic_json, atomic_npz, digest_file
from geometry_tools.model import AtlasReader

DONORS=['FJ3430','FJ3444','FJ3433','FJ3419']
GUIDES=['FJ1895','FJ2562','FJ3420','FJ3544','FJ3640','FJ3421']
TARGETS=['TA98PART:ZA:Pancreas:'+x for x in ['head','uncinate-process','neck','body','tail']]
S=1.0704662233322269
T=np.array([0.,.00040938289203434165,0.])

def sample(v,f,n,seed):
    tri=v[f];areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
    if not np.isfinite(areas).all() or areas.sum()<=0:raise ValueError('No finite source surface area')
    rng=np.random.default_rng(seed);u=(np.arange(n)+rng.random(n))/n
    ids=np.searchsorted(np.cumsum(areas)/areas.sum(),u,side='right');b=rng.random((n,2));r=np.sqrt(b[:,0]);w=np.c_[1-r,r*(1-b[:,1]),r*b[:,1]]
    return np.einsum('ni,nij->nj',w,tri[ids])

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);p.add_argument('--contract');p.add_argument('--fit',action='store_true');a=p.parse_args()
    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    source=ROOT/'public/models/atlas.json';current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    public=json.loads(source.read_text());native=json.loads(current.read_text());pr={x['id']:x for x in public['parts']};nr={x['id']:x for x in native['parts']};reader=AtlasReader()
    meshes={k:reader.read_part(public,pr[k],frame='public:BP3D:metres',catalogue_path=source,catalogue_sha=digest_file(source)) for k in DONORS+GUIDES}
    replay={}
    for k in DONORS:
        old=meshes[k];now=reader.read_part(native,nr[k],frame='model:local-ta98:assembled-metres',catalogue_path=current,catalogue_sha=digest_file(current))
        if old.faces.shape!=now.faces.shape or not np.array_equal(old.faces,now.faces):raise ValueError(k+' current/source face identity differs')
        error=np.linalg.norm(old.vertices*S+T-now.vertices,axis=1)
        replay[k]={'max_error_metres':float(error.max()),'public_geometry_sha256':old.geometry_sha,'current_geometry_sha256':now.geometry_sha,'public_bank_sha256':old.source['bank_sha256'],'current_bank_sha256':now.source['bank_sha256']}
        if error.max()>2e-7:raise ValueError(k+' does not replay declared existing height frame')
    # Remove only exact duplicated physical facets across the five finite cells.
    target=[reader.read_part(native,nr[k],frame='model:local-ta98:assembled-metres',catalogue_path=current,catalogue_sha=digest_file(current)) for k in TARGETS]
    records={}
    for mi,m in enumerate(target):
        for fi,tri in enumerate(m.triangles):
            key=tuple(sorted(tuple(x) for x in tri));records.setdefault(key,[]).append((mi,fi))
    arrays={};exterior=[];duplicate=0
    for refs in records.values():
        if len(refs)==1:
            mi,fi=refs[0];exterior.append(target[mi].triangles[fi])
        elif len(refs)==2:duplicate+=2
        else:raise ValueError('Ambiguous target shared face multiplicity')
    tv=np.asarray(exterior).reshape(-1,3);tf=np.arange(len(tv)).reshape(-1,3)
    prep={'schema':'pancreatic-bp3d-cohort-registration/v1','source_catalogue_sha256':digest_file(source),'current_catalogue_sha256':digest_file(current),'source_height_scale':S,'source_height_translation':T,'donor_replay':replay,'target_ids':TARGETS,'target_external_faces':len(tf),'target_internal_duplicate_faces_removed':duplicate,'target_geometry_sha256':{k:m.geometry_sha for k,m in zip(TARGETS,target)},'guides_only':GUIDES,'donor_ids':DONORS,'one_transform':True,'radius_changes':False,'anatomical_acceptance':False}
    atomic_json(out/'source-replay.json',prep)
    atomic_npz(out/'registration-surfaces.npz',source_vertices=meshes['FJ1895'].vertices*S+T,source_faces=meshes['FJ1895'].faces,target_vertices=tv,target_faces=tf)
    if not a.fit:return
    if not a.contract:raise ValueError('Parent frozen construction contract required before fit')
    contract=Path(a.contract);cx=json.loads(contract.read_text())
    # Contract is retained literally and source binding is checked wherever declared.
    for key in ['current_catalogue_sha256','source_catalogue_sha256']:
        if key in cx and cx[key]!=prep['current_catalogue_sha256']:raise ValueError('Contract '+key+' mismatch')
    import open3d as o3d
    sv=meshes['FJ1895'].vertices*S+T;sf=meshes['FJ1895'].faces
    fit=cx['fit'];count=int(fit['samples']);seed=int(fit['seed'])
    sp=sample(sv,sf,count,seed);tp=sample(tv,tf,count,seed+1)
    pc=lambda q:o3d.geometry.PointCloud(o3d.utility.Vector3dVector(q))
    transform=np.eye(4);levels=[]
    for voxel in fit['voxel_sizes_m']:
        sc=pc(sp).voxel_down_sample(voxel);tc=pc(tp).voxel_down_sample(voxel)
        result=o3d.pipelines.registration.registration_icp(sc,tc,float(fit['max_correspondence_m']),transform,o3d.pipelines.registration.TransformationEstimationPointToPoint(False),o3d.pipelines.registration.ICPConvergenceCriteria(relative_fitness=1e-7,relative_rmse=1e-7,max_iteration=int(fit['iterations_per_level'])))
        transform=np.asarray(result.transformation)
        levels.append({'voxel_metres':voxel,'source_points':len(sc.points),'target_points':len(tc.points),'fitness':result.fitness,'inlier_rmse_metres':result.inlier_rmse,'max_iterations':int(fit['iterations_per_level'])})
    R=transform[:3,:3];translation=transform[:3,3];angle=float(np.arccos(np.clip((np.trace(R)-1)/2,-1,1)))
    tri=sv[sf];areas=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    pivot=(tri.mean(1)*areas[:,None]).sum(0)/areas.sum();pivot_translation=pivot@R.T+translation-pivot
    max_motion={k:float(np.linalg.norm((m.vertices*S+T)@R.T+translation-(m.vertices*S+T),axis=1).max()) for k,m in meshes.items() if k!='FJ1895'}
    passed=bool(np.linalg.det(R)>0 and np.max(abs(R.T@R-np.eye(3)))<1e-8 and np.linalg.norm(pivot_translation)<=.03 and max(max_motion.values())<=.03 and angle<=np.deg2rad(20) and levels[-1]['fitness']>=.5)
    proof=dict(prep,contract_path=str(contract.resolve()),contract_sha256=digest_file(contract),script_sha256=digest_file(Path(__file__)),transform=transform,rotation_degrees=np.rad2deg(angle),world_origin_affine_translation_metres=translation,area_centroid_pivot_metres=pivot,pivot_translation_metres=pivot_translation,pivot_translation_norm_metres=np.linalg.norm(pivot_translation),max_arterial_control_displacement_metres=max_motion,sampling={'kind':'deterministic area-stratified barycentric','count_per_surface':count,'source_seed':seed,'target_seed':seed+1},levels=levels,passed=passed,limits={'max_correspondence_metres':.03,'max_pivot_translation_metres':.03,'max_arterial_control_displacement_metres':.03,'max_rotation_degrees':20,'minimum_final_sample_fitness':.5},anatomical_pending=['actual parent-mouth correspondence','fixed-context crossings and occupied volume','organ containment','render and section review'])
    atomic_json(out/'rigid-registration-proof.json',proof)
    if not passed:return
    for k,m in meshes.items():
        attrs={n:v for n,v in m.attributes.items()};arrays={'vertices':(m.vertices*S+T)@R.T+translation,'faces':m.faces,'source_vertices_height_frame':m.vertices*S+T,'source_faces':m.faces,'rigid_rotation':R,'rigid_translation':translation,'source_json':np.array(json.dumps(m.source,sort_keys=True)),'guide_only':np.array(k in GUIDES)}
        for name,v in attrs.items():arrays['source_attribute_'+name]=v
        if 'normals'in attrs:
            n=attrs['normals'].astype(float);length=np.linalg.norm(n,axis=1);arrays['normals_world_float']=(n/np.maximum(length[:,None],1e-30))@R.T
        atomic_npz(out/(k+'-rigid-candidate.npz'),**arrays)

if __name__=='__main__':main()
