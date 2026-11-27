"""Constrain storage-only slivers without moving native source controls."""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.shared_edge_collapse import collapse_short_edges

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--precheck',required=True)
ap.add_argument('--output',required=True);ap.add_argument('--origin',type=float,nargs=3,required=True)
args=ap.parse_args();started=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
d=json.loads((source/'cells.json').read_text());z=np.load(source/'cells.npz');origin=np.array(args.origin)
record=json.loads(Path(args.precheck).read_text());assert record['source_sha256']==digest_file(source/'cells.npz')
encoded,first,mapping=np.unique((z['positions']-origin).astype('<f4').astype(float),axis=0,return_index=True,return_inverse=True)
encoded+=origin;desired=z['positions'][first];surfaces={}
for sid in d['surfaces']:
    f=mapping[z[sid]];q=encoded[f];surfaces[sid]=f[np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)>0]
offset={};count=0
for sid,f in surfaces.items():offset[sid]=count;count+=len(f)
pairs=[];diagnoses=[]
for name,a in record['checks']['cells'].items():
    if a['passed']:continue
    assert a['failures']==['self_intersection_check']
    uses=d['cells'][name]
    f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
    ids=np.concatenate([offset[u['surface']]+np.arange(len(surfaces[u['surface']]))for u in uses])
    mesh=Mesh(encoded,f,'model:local-female-ta98:assembled-metres');assert mesh.geometry_sha==a['geometry_sha256']
    summary=a['self_intersections']
    if summary['defect_pairs']!=len(summary['examples']):
        summary=self_contact_summary(mesh,max_pairs=500000,max_tests=2000000,example_limit=2000)
    assert not summary['truncated']and summary['defect_pairs']==len(summary['examples'])
    pairs.extend(ids[np.asarray(summary['examples'],int)]);diagnoses.append(dict(name=name,defects=summary['defect_pairs']))
pairs=np.unique(np.sort(np.asarray(pairs,int),axis=1),axis=0);assert len(pairs)
normals={}
for sid,f in surfaces.items():
    q=desired[f];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    normals[sid]=np.repeat(n[:,None,:],3,axis=1)
# The existing native coordinates are immutable, including grey/white and
# pial interfaces. Only newly derived, sub-resolution edge controls can go.
native=np.load(ROOT/'.local-models/geometry-workbench/spinal-cord-completion/female-31-segments-v2/cells.npz')['positions']
distance,_=cKDTree(native).query(desired);locked=np.flatnonzero(distance<2.1e-12)
result=collapse_short_edges(encoded,surfaces,d['cells'],normals,desired_positions=desired,
    maximum_displacement_m=5e-8,maximum_steps=200,wall_seconds=20.,candidate_face_rings=1,
    defect_face_pairs=pairs,preserve_material_membership=False,locked=locked)
atomic_npz(out/'cells.npz',positions=encoded,**result['surfaces'])
metadata=dict(d,surfaces={sid:len(f)for sid,f in result['surfaces'].items()},construction={**d['construction'],
    'storage_sliver_repair':str(out/'conditioning.json'),'native_controls_immutable':True,'storage_origin':origin.tolist()})
atomic_json(out/'cells.json',metadata)
atomic_json(out/'conditioning.json',dict(result['evidence'],diagnoses=diagnoses,source_sha256=digest_file(source/'cells.npz'),
    immutable_native_controls=len(locked),origin=origin.tolist(),complete_surface_audit_required=True,
    anatomical_acceptance=False,installed=False,total_elapsed_seconds=time.monotonic()-started))
print(json.dumps(dict(passed=result['evidence']['passed'],contractions=len(result['evidence']['contractions']),
    remaining_pairs=len(result['evidence']['remaining_declared_defect_pairs']),
    stop=result['evidence']['budget_or_stop_reason'],elapsed_seconds=time.monotonic()-started)),flush=True)
