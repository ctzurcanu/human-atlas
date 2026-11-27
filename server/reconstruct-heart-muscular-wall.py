"""Bounded whole-wall prototype; no valve, root or tissue-layer claims."""
from pathlib import Path
import argparse,ast,json,re,struct,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--spacing-mm',type=float,default=.55);ap.add_argument('--stage',required=True);ap.add_argument('--sign-method',choices=['winding','rays'],default='winding');ap.add_argument('--unsigned-from');ap.add_argument('--field-only',action='store_true');ap.add_argument('--resume-field',action='store_true');args=ap.parse_args()
import open3d as o3d
from geometry_tools.algorithms.regular_grid_isosurface import extract
base=ROOT/'.local-models/geometry-workbench/heart';out=base/args.stage;out.mkdir(exist_ok=args.resume_field);source=base/'wall-seams-v2/wall.npz';z=np.load(source);v=z['vertices'];f=z['faces'];owners=z['face_owner'];origin=np.array([0.,1.38,.03]);start=time.monotonic();h=args.spacing_mm/1000
if not .00025<=h<=.001:raise ValueError('Explicit bounded grid spacing required')
# Cancel only literal opposed cross-piece tissue interface faces for the
# numerical whole-wall classifier. Source pieces remain unchanged on disk.
keys=np.sort(f,axis=1);u,inv,count=np.unique(keys,axis=0,return_inverse=True,return_counts=True);keep=np.ones(len(f),bool)
for key in np.flatnonzero(count==2):
 ids=np.flatnonzero(inv==key);q=v[f[ids]];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
 if owners[ids[0]]!=owners[ids[1]]and normal[0]@normal[1]<0:keep[ids]=False
ff=f[keep];face_owner=owners[keep];scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((v-origin).astype('float32')),o3d.core.Tensor(ff.astype('uint32')))
lo=np.floor((v.min(0)-.004)/h)*h;hi=np.ceil((v.max(0)+.004)/h)*h;shape=np.ceil((hi-lo)/h).astype(int)+1;shape_tuple=tuple(map(int,shape));total=int(np.prod(shape));assert total<=65000000,('Grid size exceeds explicit memory budget',total)
def read_completed_cache(path):
 # Older NumPy emitted np.int64(...) in headers when given NumPy shape
 # scalars. Parse just this literal integer wrapper; never evaluate code.
 with path.open('rb')as stream:
  version=np.lib.format.read_magic(stream);assert version==(1,0)
  length=struct.unpack('<H',stream.read(2))[0];header=stream.read(length).decode('latin1');offset=stream.tell()
 metadata=ast.literal_eval(re.sub(r'np\.int64\((\d+)\)',r'\1',header))
 assert metadata==dict(descr='<f4',fortran_order=False,shape=shape_tuple)
 assert path.stat().st_size==offset+total*4
 return np.memmap(path,mode='r',dtype='<f4',offset=offset,shape=shape_tuple)
def signed(points):
 if args.sign_method=='rays':return scene.compute_signed_distance(o3d.core.Tensor((points-origin).astype('float32')),nthreads=1,nsamples=3).numpy()
 from geometry_tools.winding_backend import winding_numbers
 unsigned=scene.compute_distance(o3d.core.Tensor((points-origin).astype('float32')),nthreads=1).numpy();w=winding_numbers(v,ff,points);return np.where(w>.5,-unsigned,unsigned)
if args.resume_field:
 completed=json.loads((out/'field-complete.json').read_text());assert completed['source_sha256']==digest_file(source)and completed['field_sha256']==digest_file(out/'distance.npy');grid=np.load(out/'distance.npy',mmap_mode='r');assert grid.shape==tuple(shape)and completed['sign_method']==args.sign_method
else:
 atomic_json(out/'plan.json',dict(source=str(source),source_sha256=digest_file(source),spacing_m=h,shape=shape.tolist(),samples=total,origin=lo.tolist(),source_faces=len(ff),cancelled_opposed_faces=int((~keep).sum()),maximum_samples=65000000,sign_method=args.sign_method,anatomical_acceptance=False))
 cached=None
 if args.unsigned_from:
  cache=base/args.unsigned_from;plan=json.loads((cache/'plan.json').read_text());done=json.loads((cache/'construction.json').read_text());assert done['constructed']and plan['source_sha256']==digest_file(source)and plan['shape']==shape.tolist()and plan['spacing_m']==h;cached=read_completed_cache(cache/'distance.npy').reshape(-1)
 grid=np.lib.format.open_memmap(out/'distance.npy',mode='w+',dtype='<f4',shape=shape_tuple);flat=grid.reshape(-1);step=180000
 for first in range(0,total,step):
  ids=np.arange(first,min(first+step,total));ijk=np.column_stack(np.unravel_index(ids,tuple(shape)));points=lo+ijk*h
  if cached is not None:
   from geometry_tools.winding_backend import winding_numbers
   distance=np.where(winding_numbers(v,ff,points)>.5,-np.abs(cached[ids]),np.abs(cached[ids]))
  else:distance=signed(points)
  flat[first:first+len(ids)]=distance
  if first//step%20==0:print('Signed field',round(100*(first+len(ids))/total,1),'percent;',round(time.monotonic()-start,1),'seconds',flush=True)
 grid.flush();atomic_json(out/'field-complete.json',dict(passed=True,source_sha256=digest_file(source),field_sha256=digest_file(out/'distance.npy'),sign_method=args.sign_method,elapsed_seconds=time.monotonic()-start,unsigned_distance_reuse=args.unsigned_from,anatomical_acceptance=False))
if args.field_only:print('Complete source-bound field saved',flush=True);raise SystemExit(0)
verts,faces,extraction=extract(grid,lo,h);atomic_npz(out/'dense-wall.npz',vertices=verts,faces=faces,frame=np.array('model:local-ta98:assembled-metres'));atomic_json(out/'extraction.json',extraction)
# Bound the local viewer cost; source deviation and emitted topology are
# measured again after decimation, rather than assuming that it is safe.
if len(faces)>100000:
 dense=o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(verts),o3d.utility.Vector3iVector(faces.astype('int32')));simple=dense.simplify_quadric_decimation(100000,maximum_error=2.5e-8);verts=np.asarray(simple.vertices).copy();faces=np.asarray(simple.triangles).copy()
# Marching cubes winding is checked against the source distance field itself.
t=verts[faces];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-30);sample=np.arange(0,len(faces),max(1,len(faces)//3000));centres=t[sample].mean(1);probe=np.r_[centres+normal[sample]*h*.2,centres-normal[sample]*h*.2];dist=signed(probe);score=float(np.median(dist[:len(sample)]-dist[len(sample):]));assert abs(score)>1e-8,'Cannot establish muscular surface orientation'
if score<0:faces=faces[:,[0,2,1]]
closest=scene.compute_closest_points(o3d.core.Tensor((verts-origin).astype('float32')),nthreads=1);distance=np.linalg.norm((closest['points'].numpy()+origin)-verts,axis=1);primitive=closest['primitive_ids'].numpy();mesh=Mesh(verts,faces,'model:local-ta98:assembled-metres');mesh.save(out/'whole-wall.npz');atomic_npz(out/'source-ownership.npz',source_vertex_owner=face_owner[primitive],source_face=primitive)
result=dict(constructed=True,extraction=extraction,source_envelope_deviation_quantiles_m=np.quantile(distance,[0,.5,.9,.99,1]).tolist(),vertices=len(verts),triangles=len(faces),candidate_geometry_sha256=mesh.geometry_sha,elapsed_seconds=time.monotonic()-start,source_wall_modified=False,installed=False,anatomical_acceptance=False,remaining=['Strict emitted topology/self-intersection audit','Whole-wall chamber and opening identity','Separate finite wall layers','Valve and great-vessel interfaces','Rendered anatomical review'])
atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
