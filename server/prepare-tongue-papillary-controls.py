"""Bind actual papillary subfacets and preserve the six earlier typed controls."""
from pathlib import Path
import argparse
import copy
import json
import sys
import time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
WORK=ROOT/'.local-models/geometry-workbench/tongue-completion'


def main(candidate,output):
    began=time.monotonic()
    if output.exists()or not output.is_relative_to(WORK):raise ValueError('Use a new local controls stage')
    report=json.loads((candidate/'construction.json').read_text())
    if not report['passed']:raise ValueError('No controls for failed actual geometry')
    muco_path=candidate/'packed/TA98TONGUE_male_mucosa.npz';muco=Mesh.load(muco_path)
    old_path=WORK/'surface-controls-shared-coverings-v1/surface-controls.json'
    previous=json.loads(old_path.read_text())
    old_mesh=Mesh.load(previous['carrierMesh'])
    if digest_file(previous['carrierMesh'])!=previous['carrierMeshSha256']:raise ValueError('Changed source controls')
    changed=report['old_source_face'];n=report['replacement_faces']
    ancestry=np.r_[np.arange(changed),np.full(n,changed),np.arange(changed+1,len(old_mesh.faces))]
    if len(ancestry)!=len(muco.faces):raise ValueError('Actual mucosal face order/count changed')
    with np.load(candidate/'relief-correspondence.npz',allow_pickle=False)as saved:
        local=saved['local_parameters'].copy();tri=saved['local_triangles'].copy()
        height=saved['height'].copy();owners=saved['papilla_type'].copy()
    plan=json.loads(Path(report['plan']['path']).read_text());basis=np.asarray(plan['basis']);center=np.asarray(plan['center'])
    groups=copy.deepcopy(previous['termGroups'])
    for row in groups:
        if 'surface'in row:
            ids=row['surface']['faces'];row['surface']['faces']=np.flatnonzero(np.isin(ancestry,ids)).tolist()
            row['surface']['sourceMeshSha256']=digest_file(muco_path)
        if 'anchor'in row:
            a=row['anchor'];old=a['face']
            if old!=changed:
                new=int(np.flatnonzero(ancestry==old)[0]);weights=np.asarray(a['barycentric'])
            else:
                uv=(np.asarray(a['point'])-center)@basis.T;new=None
                for i,t in enumerate(local[tri]):
                    b=np.linalg.solve((t[1:]-t[0]).T,uv-t[0]);w=np.r_[1-b.sum(),b]
                    if w.min()>=-1e-10:new=changed+i;weights=w;break
                if new is None:raise ValueError('Old anatomical anchor left the actual source facet')
            a.update(face=new,barycentric=weights.tolist(),point=(weights@muco.triangles[new]).tolist())
        for curve in row.get('curves',[]):
            if 'sourceMeshSha256'in curve:curve['sourceMeshSha256']=digest_file(muco_path)
        row['papillaryReliefFaceAncestryRemapped']=True
    for row in copy.deepcopy(report['termGroups']):
        chosen=np.asarray(row['surface']['faces'],int)
        best=int(max(chosen,key=lambda i:height[tri[int(i)-changed]].mean()))
        row['surface']['sourceMeshSha256']=digest_file(muco_path)
        row['anchor']=dict(partId='TA98TONGUE:male:mucosa',face=best,barycentric=[1/3]*3,
                           point=muco.triangles[best].mean(0).tolist())
        row['boundaryProvisional']=False
        row['coverageLimitation']='A correctly scaled original anterior patch; full dorsal distribution remains incomplete.'
        groups.append(row)
    output.mkdir(parents=True)
    atomic_npz(output/'source-face-ancestry.npz',source_mucosa_face_of_new_face=ancestry)
    result={**previous,'carrierMesh':str(muco_path),'carrierMeshSha256':digest_file(muco_path),
        'carrierGeometrySha256':muco.geometry_sha,'carrierBankBindingPending':True,'carrierBankHashes':{},
        'termGroups':groups,'exposedCarrierFacets':report['new_outer_faces'],
        'sourceFaceAncestry':dict(path=str(output/'source-face-ancestry.npz'),
             sha256=digest_file(output/'source-face-ancestry.npz'),array='source_mucosa_face_of_new_face'),
        'papillaryConstruction':dict(path=str(candidate/'construction.json'),sha256=digest_file(candidate/'construction.json')),
        'fullDorsalPapillaryDistribution':False,'additionalSolidMeshes':0,'additionalImageAssets':0,
        'anatomicalAcceptance':False,'references':list(dict.fromkeys(previous['references']+plan['primary_references']))}
    atomic_json(output/'surface-controls.json',result)
    errors=[]
    for row in groups:
        if 'surface'in row:
            ids=row['surface']['faces']
            if not ids or min(ids)<0 or max(ids)>=len(muco.faces):errors.append('surface:'+row['key'])
        a=row.get('anchor')
        if a:
            p=np.asarray(a['barycentric'])@muco.triangles[a['face']]
            if np.linalg.norm(p-a['point'])>1e-12:errors.append('anchor:'+row['key'])
    atomic_json(output/'checks.json',dict(passed=not errors,errors=errors,typed_controls=len(groups),
        original_controls_remapped=6,real_new_relief_controls=3,physical_meshes_added=0,
        construction_sha256=digest_file(candidate/'construction.json'),
        controls_sha256=digest_file(output/'surface-controls.json'),anatomical_acceptance=False,
        installed=False,elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(passed=not errors,typed_controls=len(groups),physical_meshes_added=0,
                         elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--candidate',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();main(a.candidate.resolve(),a.output.resolve())
