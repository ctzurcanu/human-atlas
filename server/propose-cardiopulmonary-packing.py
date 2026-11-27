"""Bounded common-lung growth and whole-heart pose proposal; never installs.

The skeleton and recipient muscles remain literal. Lung territories and their
connected airways share one positive similarity, without anisotropic stretch.
Actual triangle and occupied-volume checks remain mandatory after sampling.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
from scipy.optimize import least_squares
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.open3d_backend import library

ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args()
start=time.monotonic();path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame='model:local-male-ta98:assembled-metres'
out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def read(pid):return reader.read_part(a,parts[pid],frame=frame,catalogue_path=path)
def union(ids):
 ms=[read(i)for i in ids];off=np.r_[0,np.cumsum([len(m.vertices)for m in ms])]
 v,iv=np.unique(np.concatenate([m.vertices for m in ms]),axis=0,return_inverse=True);f=np.concatenate([iv[m.faces+off[i]]for i,m in enumerate(ms)])
 _,ix,count=np.unique(np.sort(f,axis=1),axis=0,return_index=True,return_counts=True);assert count.max()<=2
 return Mesh(v,f[ix[count==1]],frame)
lung=union(a['pulmonarySegmentIntegration']['partIds']);heart=read('TA98HEART:male:epicardium')
lp=(lung.bounds[0]+lung.bounds[1])*.5
hp=np.mean([read('TA98HEART:male:fibrous:'+s+'-annulus').vertices.mean(0)for s in ['left','right']],0)
hidden=set(a.get('defaultHidden',[]));lo,hi=lung.bounds+np.array([[-.04]*3,[.04]*3])
fixed=[p['id']for p in a['parts']if p['system']in ['skeletal','muscular']and p['id']not in hidden and not p.get('suppressed')and not p.get('defaultHidden')and np.all(np.array(p['bounds'])[1]>=lo)and np.all(np.array(p['bounds'])[0]<=hi)]
wall=union(fixed);diaphragm=union(a['finiteDiaphragmIntegration']['partIds']);o3d=library()
def index(m):
 origin=m.vertices.mean(0);scale=float(np.ptp(m.vertices,axis=0).max());s=o3d.t.geometry.RaycastingScene(nthreads=2)
 s.add_triangles(o3d.core.Tensor(((m.vertices-origin)/scale).astype(np.float32)),o3d.core.Tensor(m.faces.astype(np.uint32)))
 def d(q):return s.compute_signed_distance(o3d.core.Tensor(((q-origin)/scale).astype(np.float32)),nthreads=2,nsamples=3).numpy().astype(float)*scale
 return d,s,origin,scale
ld,_,_,_=index(lung);wd,_,_,_=index(wall);_,ds,do,dscale=index(diaphragm)
def above(q):
 origins=q.copy();origins[:,1]=2.
 rays=np.c_[(origins-do)/dscale,np.broadcast_to([0.,-1.,0.],q.shape)].astype(np.float32)
 travel=ds.cast_rays(o3d.core.Tensor(rays),nthreads=2)['t_hit'].numpy().astype(float)*dscale
 return np.where(np.isfinite(travel),q[:,1]-(2.-travel),1.)
rng=np.random.default_rng(467)
def samples(m,n):
 tri=m.triangles;area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1);i=rng.choice(len(tri),min(n,len(tri)),replace=False,p=area/area.sum())
 return np.r_[tri[i].mean(1),m.vertices[np.linspace(0,len(m.vertices)-1,n//3,dtype=int)]]
hq=samples(heart,4500);lq=samples(lung,4500)
def transformed(x):
 hr=Rotation.from_rotvec(np.deg2rad(x[3:6])).as_matrix();hs=x[6]/100;ls=x[10]/100
 return hp+(hq-hp)@hr.T*hs+x[:3]/1000,lp+(lq-lp)*ls+x[7:10]/1000
calls=0
def residual(x):
 global calls
 calls+=1;h,l=transformed(x);ls=x[10]/100
 # Query the original lung field in inverse candidate coordinates.
 hd=ld(lp+(h-lp-x[7:10]/1000)/ls)*ls
 return np.r_[np.maximum(0,.0007-hd)*1000,np.maximum(0,.0007-wd(h))*3000,np.maximum(0,.0005-above(h))*3000,np.maximum(0,.0004-wd(l))*3000,np.maximum(0,.0004-above(l))*3000,x[:3]*.002,x[3:6]*.002,(x[6]-100)*.004,x[7:10]*.003,(x[10]-110)*.008]
def jac(x):
 base=residual(x);step=.05
 return np.column_stack([(residual(x+np.eye(11)[i]*step)-base)/step for i in range(11)])
x0=np.array([5,10,0,5,-10,-5,94,0,0,0,100.001],float)
solve=least_squares(residual,x0,jac=jac,bounds=([-30,-30,-30,-35,-35,-35,90,-10,-10,-10,100],[30,30,30,35,35,35,100.0001,10,10,10,115]),max_nfev=70,ftol=1e-6,xtol=1e-6,gtol=1e-6)
x=solve.x;h,l=transformed(x);ls=x[10]/100;hl=ld(lp+(h-lp-x[7:10]/1000)/ls)*ls
def matrix(p,r,s,t):
 m=np.eye(4);m[:3,:3]=s*r;m[:3,3]=p-s*r@p+t;return m.tolist()
report=dict(source_catalogue_sha256=digest_file(path),heart_matrix=matrix(hp,Rotation.from_rotvec(np.deg2rad(x[3:6])).as_matrix(),x[6]/100,x[:3]/1000),lung_matrix=matrix(lp,np.eye(3),ls,x[7:10]/1000),parameters=x.tolist(),heart_relative_scale=float(x[6]/100),lung_relative_scale=float(ls),minimum_heart_lung_sample_m=float(hl.min()),minimum_heart_wall_sample_m=float(wd(h).min()),minimum_lung_wall_sample_m=float(wd(l).min()),minimum_heart_above_diaphragm_m=float(above(h).min()),minimum_lung_above_diaphragm_m=float(above(l).min()),heart_lung_penetrating_samples=int((hl<0).sum()),lung_wall_penetrating_samples=int((wd(l)<0).sum()),optimizer_status=int(solve.status),evaluations=calls,elapsed_seconds=time.monotonic()-start,proposal_only=True,installed=False)
atomic_json(out/'proposal.json',report);atomic_npz(out/'samples.npz',heart_before=hq,heart_after=h,lung_before=lq,lung_after=l);print(json.dumps(report,indent=2),flush=True)
