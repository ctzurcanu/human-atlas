"""Screened Poisson candidate from literal native oriented tissue samples."""
from pathlib import Path
import json,sys,time
import numpy as np
import open3d as o3d
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';p=b/'wall-seams-v2/wall.npz';out=b/'whole-wall-v13';out.mkdir(exist_ok=False);z=np.load(p);v=z['vertices'];f=z['faces'];owners=z['face_owner'];o=np.array([0,1.38,.03]);keys=np.sort(f,axis=1);_,inv,count=np.unique(keys,axis=0,return_inverse=True,return_counts=True);keep=np.ones(len(f),bool)
for key in np.flatnonzero(count==2):
 ids=np.flatnonzero(inv==key);t=v[f[ids]];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
 if owners[ids[0]]!=owners[ids[1]]and n[0]@n[1]<0:keep[ids]=False
source=o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v-o),o3d.utility.Vector3iVector(f[keep].astype('int32')));source.compute_triangle_normals();o3d.utility.random.seed(98);cloud=source.sample_points_uniformly(150000,use_triangle_normal=True);surface,density=o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(cloud,depth=9,scale=1.05,linear_fit=True,n_threads=1);vv=np.asarray(surface.vertices).copy();ff=np.asarray(surface.triangles).copy();Mesh(vv+o,ff,'model:local-ta98:assembled-metres').save(out/'dense-wall.npz');print('Poisson contour saved',len(ff),'faces',round(time.monotonic()-start,1),'s',flush=True)
solid=md.Manifold(md.Mesh64(vv,ff.astype('uint64')));status=str(solid.status());atomic_json(out/'contour.json',dict(status=status,elapsed_seconds=time.monotonic()-start));assert status=='Error.NoError',status
pieces=solid.decompose();vol=np.array([q.volume()for q in pieces]);main=int(vol.argmax());component_rows=[dict(volume_mL=float(q.volume()*1e6),genus=int(q.genus()),triangles=q.num_tri())for q in pieces];atomic_json(out/'components.json',component_rows)
assert all(abs(x)*1e6<.1 for i,x in enumerate(vol)if i!=main),'Large unreviewed Poisson tissue component'
simple=pieces[main].set_tolerance(.000001).simplify(.00012);packed=simple.to_mesh64();relative=packed.vert_properties[:,:3];faces=packed.tri_verts.astype('int64');mesh=Mesh(relative.astype('float32').astype(float)+o,faces,'model:local-ta98:assembled-metres');mesh.save(out/'whole-wall.npz');scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((v-o).astype('float32')),o3d.core.Tensor(f[keep].astype('uint32')));distance=scene.compute_distance(o3d.core.Tensor((mesh.vertices-o).astype('float32')),nthreads=1).numpy();print('Encoded contour ready',len(faces),'faces',round(time.monotonic()-start,1),'s',flush=True);proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof)
report=dict(passed=proof['passed'],source_sha256=digest_file(p),point_samples=150000,poisson_depth=9,vertices=len(relative),triangles=len(faces),genus=int(simple.genus()),volume_mL=simple.volume()*1e6,source_deviation_quantiles_m=np.quantile(distance,[0,.5,.9,.99,1]).tolist(),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
