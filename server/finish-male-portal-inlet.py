"""Connect the actual foregut portal rim to its finite intrahepatic mouth."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.annular_bridge import find_annular_patch,attach_bridge
b=ROOT/'.local-models/geometry-workbench/liver-completion';out=b/'male-portal-inlet-v263';out.mkdir(exist_ok=False);start=time.monotonic();path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());r=AtlasReader();parts={p['id']:p for p in a['parts']};frame='model:local-ta98:assembled-metres';source=r.read_part(a,parts['TA98HEPATIC:male:portal-distal'],frame=frame,catalogue_path=path);v,ix=np.unique(source.vertices,axis=0,return_inverse=True);f=ix[source.faces]
from geometry_tools.algorithms.shared_closed_cell_clip import SharedCutBank
bank=SharedCutBank(v);upper=[];lower=[];limit=1.232
for t in f:
 for sign,target in [(1,upper),(-1,lower)]:
  poly=[]
  for i in range(3):
   x,y=int(t[i]),int(t[(i+1)%3]);dx,dy=v[x,1]-limit,v[y,1]-limit
   if sign*dx>=0:poly.append(x)
   if dx*dy<0:poly.append(bank.intersection(x,y,np.array([0.,1.,0.]),limit,'portal-neck'))
  for i in range(1,len(poly)-1):target.append([poly[0],poly[i],poly[i+1]])
v=np.asarray(bank.positions);f=np.asarray(upper,int);ls=boundary_loops(f);root=next(loop for loop in ls if abs(v[loop,1].mean()-limit)<1e-9);Mesh(v,np.asarray(lower,int),frame).save(out/'confluence-context.npz');outer=v[root];center=outer.mean(0);pv=np.r_[outer,center+.93*(outer-center)];n=len(root);pf=np.array([[i,(i+1)%n,n+(i+1)%n]for i in range(n)]+[[i,n+(i+1)%n,n+i]for i in range(n)]);tri=pv[pf];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]).sum(0)
if normal[1]<0:pf=pf[:,[0,2,1]]
child=r.read_part(a,parts['TA98HEPATIC:male:portal'],frame=frame,catalogue_path=path);cc=find_annular_patch(child.vertices,child.faces,[-.0058649,1.242,.032422],[0,1,0]);q=attach_bridge(pv,pf,np.arange(len(pf)),child.vertices,child.faces,cc,stations=16,controls=32,straight_shear=True);q['faces']=np.r_[q['faces'],pf[:,[0,2,1]]];origin=np.array([0,1.25,0]);m=Mesh((q['vertices']-origin).astype('<f4').astype(float)+origin,q['faces'],frame);m.save(out/'parent.npz');check=audit(m,require_closed=True,max_pairs=500000,max_tests=8000000);atomic_json(out/'parent-check.json',check);atomic_json(out/'junction.json',dict(passed=check['passed'],role='portal',parent_id='TA98HEPATIC:male:portal-distal',parent_interface_faces=q['interface_faces'].tolist(),child_interface_faces=cc.tolist(),native_foregut_outer_rim_literal=True,lumen_disks_added=0,elapsed_seconds=time.monotonic()-start));print(check['passed'],check['failures'],flush=True)
if not check['passed']:raise ValueError('Portal inlet failed; not installed')
