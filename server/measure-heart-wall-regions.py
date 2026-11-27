"""Source-bound chamber ownership and external/internal accessibility charts.

Outward hemisphere visibility distinguishes outer heart covering from blood
cavity lining. It is an explicit surface classifier, not histological proof.
"""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
ap=argparse.ArgumentParser();ap.add_argument('--source',default='whole-wall-v18');ap.add_argument('--stage',default='wall-regions-v20');args=ap.parse_args()
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/args.source;assert json.loads((source/'construction.json').read_text())['passed'];out=b/args.stage;out.mkdir(exist_ok=False);m=Mesh.load(source/'whole-wall.npz');q=m.triangles;c=q.mean(1);n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];o=np.array([0.,1.375,.03125]);scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces.astype('uint32')))
axis=np.eye(3)[np.argmin(abs(n),axis=1)];u=np.cross(n,axis);u/=np.linalg.norm(u,axis=1)[:,None];v=np.cross(n,u);escaped=np.zeros(len(q),int);rays_per_face=64
for i in range(rays_per_face):
 r=np.sqrt((i+.5)/rays_per_face);angle=i*np.pi*(3-np.sqrt(5));direction=r*np.cos(angle)*u+r*np.sin(angle)*v+np.sqrt(1-r*r)*n;rays=np.c_[c-o+n*.0000002,direction].astype('float32');hit=scene.cast_rays(o3d.core.Tensor(rays),nthreads=1)['t_hit'].numpy();escaped+=~np.isfinite(hit)
visibility=escaped/rays_per_face;inner=visibility<.5;native=np.load(b/'wall-seams-v2/wall.npz');src=o3d.t.geometry.RaycastingScene(nthreads=1);src.add_triangles(o3d.core.Tensor((native['vertices']-o).astype('float32')),o3d.core.Tensor(native['faces'].astype('uint32')));near=src.compute_closest_points(o3d.core.Tensor((c-o).astype('float32')),nthreads=1);owner=native['face_owner'][near['primitive_ids'].numpy()];atomic_npz(out/'charts.npz',face_owner=owner,endocardial_faces=inner,external_visibility=visibility,face_centres=c,face_normals=n)
areas=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)/2;rows={}
for i,name in enumerate(['right-atrium','right-ventricle','left-atrium','left-ventricle']):
 rows[name]=dict(endocardial_faces=int(np.sum((owner==i)&inner)),epicardial_faces=int(np.sum((owner==i)&~inner)),endocardial_area_m2=float(areas[(owner==i)&inner].sum()),epicardial_area_m2=float(areas[(owner==i)&~inner].sum()))
report=dict(computed=True,source_geometry_sha256=m.geometry_sha,source_archive_sha256=digest_file(source/'whole-wall.npz'),source_ownership_sha256=digest_file(b/'wall-seams-v2/wall.npz'),rays_per_face=rays_per_face,internal_visibility_cutoff=.5,ambiguous_faces=int(np.sum((visibility>.3)&(visibility<.7))),chambers=rows,visibility_quantiles=np.quantile(visibility,[0,.1,.25,.5,.75,.9,1]).tolist(),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False);atomic_json(out/'classification.json',report);print(json.dumps(report),flush=True)
