"""Audit a true eroded core and its complementary finite surface covering."""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';out=base/'surface-covering-v44';out.mkdir(exist_ok=False);source=base/'surface-covering-v43/core.npz';m=Mesh.load(source);a=Mesh.load(base/'whole-wall-v25/whole-wall.npz');o=np.array([0.,1.375,.03125]);s=md.Manifold(md.Mesh64(m.vertices-o,m.faces.astype('uint64')));assert str(s.status())=='Error.NoError';simple=s.set_tolerance(.0000001).simplify(.000001);z=simple.to_mesh64();conditioned=normalize_shared_boundary(z.vert_properties[:,:3]+o,{'core':z.tri_verts},o);core=Mesh(conditioned['positions'],conditioned['cells']['core'],m.frame);core.save(out/'core.npz');atomic_json(out/'encoding.json',conditioned['evidence']);print('Conditioned eroded core',len(core.faces),'faces',round(time.monotonic()-start,1),'s',flush=True)
cov=Mesh(np.r_[a.vertices,core.vertices],np.r_[a.faces,core.faces[:,[0,2,1]]+len(a.vertices)],m.frame);cov.save(out/'covering.npz');checks={}
for name,mesh in [('core',core),('covering',cov)]:
 proof=audit(mesh,require_closed=True,max_pairs=5000000,max_tests=30000000);atomic_json(out/(name+'-audit.json'),proof);checks[name]=dict(passed=proof['passed'],failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs']);print(name,json.dumps(checks[name]),flush=True)
scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((a.vertices-o).astype('float32')),o3d.core.Tensor(a.faces.astype('uint32')));distance=scene.compute_signed_distance(o3d.core.Tensor((core.vertices-o).astype('float32')),nthreads=1,nsamples=3).numpy();inside=bool(np.all(distance<0));report=dict(passed=all(x['passed']for x in checks.values())and inside,source_sha256=digest_file(source),outer_geometry_sha256=a.geometry_sha,inner_geometry_sha256=core.geometry_sha,simplification_tolerance_m=.000001,core_vertices=len(core.vertices),core_triangles=len(core.faces),covering_triangles=len(cov.faces),core_vertices_inside_outer=inside,signed_distance_to_outer_quantiles_m=np.quantile(distance,[0,.1,.5,.9,1]).tolist(),checks=checks,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
