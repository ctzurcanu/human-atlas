"""Bounded body expansion; actual source mouths are held fixed."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.gastric_body_field import expand_body
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/'male-gastric-volume-v268';out.mkdir(exist_ok=False);r=AtlasReader();start=time.monotonic();a=json.loads((b/'final-ready-v265/male.json').read_text());parts={p['id']:p for p in a['parts']};tri=[]
for region in ['007','009','012','014']:
 p=parts['GASTRIC:male:019:'+region];m=r.read_part(a,p,frame='model:local-ta98:assembled-metres');tri.append(m.triangles[:p['tissueDrawGroups'][0]['count']//3])
v,inv=np.unique(np.concatenate(tri).reshape(-1,3),axis=0,return_inverse=True);f=inv.reshape(-1,3);loops=boundary_loops(f)
assert len(loops)==2,('Expected two computational volume mouths',len(loops))
def volume(vertices):
 triangles=[vertices[f]]
 for loop in loops:
  rim=vertices[loop];centre=rim.mean(0);triangles.append(np.array([[centre,rim[(i+1)%len(rim)],rim[i]]for i in range(len(rim))]))
 t=np.concatenate(triangles)-vertices.mean(0);return float(abs(np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()/6))
baseline=volume(v);targets=[];hidden=set(a.get('defaultHidden',[]));liver=parts['TA98PART:ZA:Liver:capsule'];lm=r.read_part(a,liver,frame='model:local-ta98:assembled-metres');targets.append((liver['id'],Mesh(lm.vertices,lm.faces[:liver['liverReferenceTriangleCount']],lm.frame)))
for p in a['parts']:
 if p.get('suppressed')or p.get('defaultHidden')or p['id']in hidden:continue
 if p['id'].startswith('GASTRIC:')or p['id'].startswith('TA98PART:ZA:Liver:')or p['id'].startswith('TA98HEPATIC:')or p['id'].startswith('TA98LIVER:'):continue
 if p.get('system')not in ['digestive','urinary','lymphatic','skeletal']and p['id']not in a['finiteDiaphragmIntegration']['partIds']:continue
 bounds=np.asarray(p['bounds']);expanded=np.array([v.min(0),v.max(0)]);expanded[1,0]+=.030
 if np.any(bounds[1]<expanded[0])or np.any(bounds[0]>expanded[1]):continue
 targets.append((p['id'],r.read_part(a,p,frame='model:local-ta98:assembled-metres')))
targets=[(pid,m,BVH(m))for pid,m in targets];rows=[];chosen=None
for strength in [0.,.20,.35,.50]:
 q,J=expand_body(v,strength);m=Mesh(q,f,'model:local-ta98:assembled-metres');bv=BVH(m);checks=[]
 for pid,peer,pb in targets:
  _,proof=pair_contacts(m,peer,0,first_bvh=bv,second_bvh=pb,max_pairs=300000,max_tests=6000000);checks.append(dict(part=pid,crossings=proof['crossing_face_pairs'],coplanar=proof['coplanar_overlap_face_pairs'],truncated=proof['truncated']))
 report=dict(strength=strength,volume_ml=volume(q)*1e6,volume_ratio=volume(q)/baseline,checks=checks,maximum_movement_m=float(np.linalg.norm(q-v,axis=1).max()),minimum_jacobian=float(np.linalg.det(J).min()))
 rows.append(report);print(strength,report['volume_ml'],[(x['part'],x['crossings'])for x in checks if x['crossings']],flush=True)
 if strength and all(x['crossings']<=y['crossings']and not x['truncated']for x,y in zip(checks,rows[0]['checks'])):chosen=report
atomic_json(out/'fit.json',dict(passed=chosen is not None,chosen=chosen,trials=rows,computational_mouth_closures_only=True,rendered_lumen_disks=False,elapsed_seconds=time.monotonic()-start))
if chosen is None:raise ValueError('No size increase passed current contact bounds')
