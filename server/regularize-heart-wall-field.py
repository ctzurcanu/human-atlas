"""Close only sub-millimetre source seam defects, then contour natively.

The physical tolerance is explicit. This is a candidate, never an accepted
heart: roots, chamber continuity and the rendered sections require review.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from scipy.ndimage import grey_opening,label
import point_cloud_utils as pcu
import manifold3d as manifold
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);ap.add_argument('--source',default='whole-wall-v6');args=ap.parse_args()
base=ROOT/'.local-models/geometry-workbench/heart';source=base/args.source;out=base/args.stage;out.mkdir(exist_ok=False);start=time.monotonic();plan=json.loads((source/'plan.json').read_text());done=json.loads((source/'field-complete.json').read_text());assert done['field_sha256']==digest_file(source/'distance.npy')and done['sign_method']=='winding'
d=np.load(source/'distance.npy');h=plan['spacing_m'];origin=np.array(plan['origin']);shape=np.array(d.shape);radius=1.5*h
xyz=np.stack(np.meshgrid(*[np.arange(-1,2)]*3,indexing='ij'),axis=-1);footprint=np.linalg.norm(xyz,axis=-1)<=1.5
# With negative tissue SDF, min then max closes positive cracks. The
# footprint radius is .825 mm on the current .55 mm source grid.
closed=grey_opening(d,footprint=footprint,mode='nearest');lab,n=label(closed<0);counts=np.bincount(lab.ravel());counts[0]=0;main=int(counts.argmax());discarded=np.flatnonzero((counts>0)&(np.arange(len(counts))!=main));assert all(counts[q]<=8 for q in discarded),'Non-numerical source tissue component requires review'
closed[(lab!=0)&(lab!=main)]=np.abs(closed[(lab!=0)&(lab!=main)]);del lab
atomic_npz(out/'field.npz',distance=closed,origin=origin,spacing=np.array(h))
lo=np.full(tuple(shape-1),np.inf,dtype='<f4');hi=-lo.copy();offsets=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
for off in offsets:
 q=closed[tuple(slice(int(k),int(k+s-1))for k,s in zip(off,shape))];np.minimum(lo,q,out=lo);np.maximum(hi,q,out=hi)
cells=np.column_stack(np.nonzero((lo<0)&(hi>=0)));del lo,hi
assert len(cells)<1000000
strides=np.array([shape[1]*shape[2],shape[2],1]);cube=(cells@strides)[:,None]+(offsets@strides)[None,:];ids,inv=np.unique(cube,return_inverse=True);ijk=np.column_stack(np.unravel_index(ids,tuple(shape)));positions=ijk*h
v,f=pcu.marching_cubes_sparse_voxel_grid(closed.reshape(-1)[ids].astype(float),positions,inv.reshape(-1,8).astype('int64'),0.)
# Orient by the signed integral before the native manifold validator.
t=v[f];volume=float(np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()/6)
if volume<0:f=f[:,[0,2,1]]
solid=manifold.Manifold(manifold.Mesh64(v.astype(float),f.astype('uint64')))
status=str(solid.status());atomic_json(out/'contour.json',dict(status=status,vertices=len(v),triangles=len(f),active_cells=len(cells),elapsed_seconds=time.monotonic()-start));assert status=='Error.NoError',status
# Native tolerance simplification retains manifoldness and limits movement.
simple=solid.simplify(.00012);encoded=simple.to_mesh64();v=encoded.vert_properties[:,:3]+origin;f=encoded.tri_verts.astype('int64');mesh=Mesh(v,f,'model:local-ta98:assembled-metres');mesh.save(out/'whole-wall.npz')
result=dict(constructed=True,source_field_sha256=done['field_sha256'],source_wall_sha256=done['source_sha256'],closing_radius_m=radius,simplification_tolerance_m=.00012,removed_numerical_samples=int(counts[discarded].sum()),tissue_volume_mL=float(simple.volume()*1e6),genus=int(simple.genus()),components=len(simple.decompose()),vertices=len(v),triangles=len(f),added_tissue_mL=float(((closed<0)&(d>=0)).sum()*h**3*1e6),removed_tissue_mL=float(((closed>=0)&(d<0)).sum()*h**3*1e6),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False)
atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
