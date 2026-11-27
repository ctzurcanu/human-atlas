"""Give one remaining near-collinear carrier facet a finite altitude."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.common import atomic_json,digest_file
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'whole-wall-v24/whole-wall.npz';out=b/'whole-wall-v25';out.mkdir(exist_ok=False);m=Mesh.load(source);v=m.vertices.copy();f=m.faces;o=np.array([0.,1.375,.03125]);q=v[f];oldnorm=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);edges=q-np.roll(q,1,axis=1);height=np.linalg.norm(oldnorm,axis=1)/np.linalg.norm(edges,axis=2).max(1);bad=np.flatnonzero(height<.00001);assert len(bad)==1
i=int(bad[0]);tri=f[i];pairs=[(0,1,2),(1,2,0),(2,0,1)];a,c,j=max(pairs,key=lambda row:np.linalg.norm(v[tri[row[0]]]-v[tri[row[1]]]));pa,pc,pj=v[tri[[a,c,j]]];edge=pc-pa;projection=pa+edge*((pj-pa)@edge)/(edge@edge);direction=pj-projection;length=np.linalg.norm(direction);assert length>0;direction/=length;delta=direction*(.00002-length);candidates=[]
for pattern in [{int(tri[j]):1.},{int(tri[a]):-1.},{int(tri[c]):-1.},{int(tri[a]):-.5,int(tri[c]):-.5}]:
 ids=np.array(list(pattern));trial=v.copy()
 for k,weight in pattern.items():trial[k]=(trial[k]-o+delta*weight).astype('float32').astype(float)+o
 affected=np.flatnonzero(np.isin(f,ids).any(1));t=trial[f[affected]];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);movement=np.linalg.norm(trial[ids]-v[ids],axis=1)
 if movement.max()>.00003 or np.any(np.einsum('ij,ij->i',normal,oldnorm[affected])<=0):continue
 t=trial[f[i]];newheight=np.linalg.norm(np.cross(t[1]-t[0],t[2]-t[0]))/np.linalg.norm(t-np.roll(t,1,axis=0),axis=1).max()
 if newheight<.00001:continue
 candidates.append((float(movement.max()),trial,ids,newheight))
assert candidates,'No orientation-preserving bounded sliver conditioning';maximum,v,ids,newheight=min(candidates,key=lambda p:p[0]);mesh=Mesh(v,f,m.frame);mesh.save(out/'whole-wall.npz');proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof);result=dict(passed=proof['passed'],source_sha256=digest_file(source),changed_controls=ids.tolist(),maximum_control_displacement_m=maximum,source_face=i,old_altitude_m=float(height[i]),new_altitude_m=float(newheight),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
