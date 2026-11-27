"""Bounded pose proposal; native distance samples are not an acceptance gate.

All recipient tissues remain fixed. A later literal-triangle and occupied-volume
audit is required before any pose may be installed.
"""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
from scipy.optimize import least_squares,differential_evolution
from scipy.spatial.transform import Rotation
from scipy.ndimage import map_coordinates
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.common import atomic_json,digest_file,atomic_npz
from geometry_tools.open3d_backend import library

ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--stage',required=True);ap.add_argument('--catalogue');ap.add_argument('--seed',type=float,nargs=7);ap.add_argument('--global-search',action='store_true');ap.add_argument('--surface-wall-only',action='store_true');ap.add_argument('--rotation-limit-deg',type=float,default=15);ap.add_argument('--minimum-scale',type=float,default=.9);args=ap.parse_args()
assert 15<=args.rotation_limit_deg<=35 and .8<=args.minimum_scale<=1
start=time.monotonic();path=ROOT/(args.catalogue or f'.local-models/ta98-review/bowel-wall-local-atlas/{args.sex}/atlas.json');a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres'
out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
def read(pid):return reader.read_part(a,parts[pid],frame=frame,catalogue_path=path)
heart=read(f'TA98HEART:{args.sex}:epicardium')
annuli=[read(f'TA98HEART:{args.sex}:fibrous:{s}-annulus')for s in ['right','left']];pivot=np.mean([m.vertices.mean(0)for m in annuli],0)
def union(ids):
 ms=[read(i)for i in ids];offsets=np.r_[0,np.cumsum([len(m.vertices)for m in ms])]
 v,iv=np.unique(np.concatenate([m.vertices for m in ms]),axis=0,return_inverse=True)
 f=np.concatenate([iv[m.faces+offsets[i]]for i,m in enumerate(ms)])
 _,ix,count=np.unique(np.sort(f,axis=1),axis=0,return_index=True,return_counts=True)
 assert count.max()<=2
 return Mesh(v,f[ix[count==1]],frame)
lung=union(a['pulmonarySegmentIntegration']['partIds'])
edges=np.concatenate([lung.faces[:,[0,1]],lung.faces[:,[1,2]],lung.faces[:,[2,0]]]);_,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
assert np.all(counts==2),'Pulmonary envelope must be closed before signed-distance proposal'
hidden=set(a.get('defaultHidden',[]));lo,hi=heart.bounds+np.array([[-.025]*3,[.025]*3])
wall_ids=[p['id']for p in a['parts']if p['id']not in hidden and not p.get('suppressed')and not p.get('defaultHidden')and p['system']in ['muscular','skeletal']and np.all(np.array(p['bounds'])[1]>=lo)and np.all(np.array(p['bounds'])[0]<=hi)]
wall=union(wall_ids)
o3d=library()
def index(mesh,signed=True):
 origin=mesh.vertices.mean(0);scale=float(np.ptp(mesh.vertices,axis=0).max());v=((mesh.vertices-origin)/scale).astype(np.float32);f=mesh.faces.astype(np.uint32)
 scene=o3d.t.geometry.RaycastingScene(nthreads=2);scene.add_triangles(o3d.core.Tensor(v),o3d.core.Tensor(f))
 def distance(q):
  tensor=o3d.core.Tensor(((q-origin)/scale).astype(np.float32))
  result=scene.compute_signed_distance(tensor,nthreads=2,nsamples=3)if signed else scene.compute_distance(tensor,nthreads=2)
  return result.numpy().astype(float)*scale
 return distance
lung_distance=index(lung);wall_distance=index(wall,signed=not args.surface_wall_only)
diaphragm=union(a['finiteDiaphragmIntegration']['partIds'])
diaphragm_origin=diaphragm.vertices.mean(0);diaphragm_scale=float(np.ptp(diaphragm.vertices,axis=0).max())
diaphragm_scene=o3d.t.geometry.RaycastingScene(nthreads=2)
diaphragm_scene.add_triangles(o3d.core.Tensor(((diaphragm.vertices-diaphragm_origin)/diaphragm_scale).astype(np.float32)),o3d.core.Tensor(diaphragm.faces.astype(np.uint32)))
def above_diaphragm(q):
 if args.surface_wall_only:return np.ones(len(q))
 origins=q.copy();origins[:,1]=2.
 rays=np.c_[(origins-diaphragm_origin)/diaphragm_scale,np.broadcast_to([0.,-1.,0.],q.shape)].astype(np.float32)
 travel=diaphragm_scene.cast_rays(o3d.core.Tensor(rays),nthreads=2)['t_hit'].numpy().astype(float)*diaphragm_scale
 return np.where(np.isfinite(travel),q[:,1]-(2.-travel),1.)
rng=np.random.default_rng(437);tri=heart.triangles;area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
face_ids=rng.choice(len(tri),size=5000,replace=False,p=area/area.sum());q=np.concatenate([tri[face_ids].mean(1),heart.vertices[np.linspace(0,len(heart.vertices)-1,2000,dtype=int)]])
def pose(x):return pivot+(q-pivot)@Rotation.from_rotvec(np.deg2rad(x[3:6])).as_matrix().T*(x[6]/100)+x[:3]/1000
calls=0
def residual(x):
 global calls
 calls+=1;p=pose(x);ld=lung_distance(p);wd=wall_distance(p);dd=above_diaphragm(p)
 return np.r_[np.maximum(0,.001-ld)*1000,np.maximum(0,.001-wd)*3000,np.maximum(0,.0005-dd)*3000,x[:3]*.003,x[3:6]*.003,(x[6]-100)*.008]
def jac(x):
 base=residual(x);step=.05
 return np.column_stack([(residual(x+np.eye(7)[i]*step)-base)/step for i in range(7)])
x0=np.array(args.seed if args.seed else [0.,0.,0.,0.,0.,0.,100.]);initial=residual(np.array([0.,0.,0.,0.,0.,0.,100.]))
r=args.rotation_limit_deg
lower=np.array([-30.,-5. if args.surface_wall_only else -30.,-30.,-r,-r,-r,args.minimum_scale*100]);upper=np.array([30.,15. if args.surface_wall_only else 30.,30.,r,r,r,100.0001])
if args.global_search:
 full_q=q.copy();q=full_q[np.linspace(0,len(full_q)-1,1000,dtype=int)]
 # Build distance fields once for the broad search. Querying the native
 # triangle scene for every individual trial wasted most of the old budget.
 # Trilinear fields only propose poses; the final refinement and independent
 # emitted-triangle gates still use the actual meshes.
 radius=float(np.linalg.norm(full_q-pivot,axis=1).max())+.035
 step=.002;grid_origin=pivot-radius;grid_shape=np.ceil(np.full(3,2*radius)/step).astype(int)+1
 coords=np.stack(np.meshgrid(*[grid_origin[i]+np.arange(grid_shape[i])*step for i in range(3)],indexing='ij'),axis=-1).reshape(-1,3)
 def field(distance):
  values=np.concatenate([distance(coords[i:i+100000])for i in range(0,len(coords),100000)])
  return values.reshape(tuple(grid_shape)).astype(np.float32)
 lung_field=field(lung_distance);wall_field=field(wall_distance)
 xz=np.stack(np.meshgrid(grid_origin[0]+np.arange(grid_shape[0])*step,grid_origin[2]+np.arange(grid_shape[2])*step,indexing='ij'),axis=-1).reshape(-1,2)
 query=np.c_[xz[:,0],np.full(len(xz),2.),xz[:,1]];heights=(2.-above_diaphragm(query)).reshape(grid_shape[0],grid_shape[2]);del coords
 def coarse_residual(x):
  p=pose(x);c=((p-grid_origin)/step).T;ld=map_coordinates(lung_field,c,order=1,mode='constant',cval=-1.,prefilter=False);wd=map_coordinates(wall_field,c,order=1,mode='constant',cval=-1.,prefilter=False);dh=map_coordinates(heights,c[[0,2]],order=1,mode='constant',cval=0.,prefilter=False)
  return np.r_[np.maximum(0,.0015-ld)*1000,np.maximum(0,.001-wd)*3000,np.maximum(0,.0005-(p[:,1]-dh))*3000,x[:3]*.003,x[3:6]*.003,(x[6]-100)*.008]
 def objective(x):
  v=coarse_residual(x);return float(v@v)
 broad=differential_evolution(objective,list(zip(lower,upper)),maxiter=65,popsize=6,polish=False,seed=472,tol=.001,workers=1,updating='immediate')
 atomic_json(out/'broad-search.json',dict(grid_spacing_m=step,grid_shape=grid_shape.tolist(),parameters=broad.x.tolist(),objective=float(broad.fun),evaluations=int(broad.nfev),elapsed_seconds=time.monotonic()-start,proposal_only=True));print('Saved bounded grid search; refining against actual triangles',flush=True)
 x0=broad.x;q=full_q
solve=least_squares(residual,x0,jac=jac,bounds=(lower,upper),max_nfev=60,ftol=1e-6,xtol=1e-6,gtol=1e-6)
x=solve.x;p=pose(x);ld=lung_distance(p);wd=wall_distance(p)
atomic_npz(out/'sample-points.npz',before=q,after=p,lung_signed_distance=ld,chest_wall_signed_distance=wd)
rotation=Rotation.from_rotvec(np.deg2rad(x[3:6])).as_matrix();scale=float(x[6]/100);translation=x[:3]/1000
matrix=np.eye(4);matrix[:3,:3]=scale*rotation;matrix[:3,3]=pivot-scale*rotation@pivot+translation
report=dict(source_catalogue_sha256=digest_file(path),sex=args.sex,pivot_m=pivot.tolist(),matrix=matrix.tolist(),translation_m=translation.tolist(),rotation_degrees=x[3:6].tolist(),relative_scale=scale,minimum_lung_sample_distance_m=float(ld.min()),minimum_wall_sample_distance_m=float(wd.min()),wall_distance_method='unsigned-to-actual-surfaces'if args.surface_wall_only else'signed-union',diaphragm_vertical_above_constraint=not args.surface_wall_only,lung_penetrating_samples=int((ld<0).sum()),wall_penetrating_samples=int((wd<0).sum()),initial_maximum_penetration_mm=float(initial[:2*len(q)].max()),optimizer_status=int(solve.status),evaluations=calls,elapsed_seconds=time.monotonic()-start,proposal_only=True,actual_triangle_and_volume_checks_required=True,installed=False)
atomic_json(out/'proposal.json',report);print(json.dumps(report,indent=2),flush=True)
