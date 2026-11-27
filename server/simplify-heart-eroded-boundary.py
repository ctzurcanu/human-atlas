"""Remove sub-micrometre offset slivers before actual Float32 validation."""
from pathlib import Path
import sys,time,json
import numpy as np
import manifold3d as md,open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'surface-covering-v49';out.mkdir(exist_ok=False);m=Mesh.load(b/'surface-covering-v43/core.npz');a=Mesh.load(b/'whole-wall-v25/whole-wall.npz');o=np.array([0.,1.375,.03125]);solid=md.Manifold(md.Mesh64(m.vertices-o,m.faces.astype('uint64')));simple=solid.set_tolerance(.00001).simplify(.00001);z=simple.to_mesh64();r=normalize_shared_boundary(z.vert_properties[:,:3]+o,{'core':z.tri_verts},o);core=Mesh(r['positions'],r['cells']['core'],m.frame);core.save(out/'core.npz');atomic_json(out/'encoding.json',r['evidence']);print('Offset sliver conditioning',len(core.faces),'triangles',round(time.monotonic()-start,2),'s',flush=True);proof=audit(core,require_closed=True,max_pairs=4500000,max_tests=28000000);atomic_json(out/'core-audit.json',proof)
scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((a.vertices-o).astype('float32')),o3d.core.Tensor(a.faces.astype('uint32')));d=scene.compute_signed_distance(o3d.core.Tensor((core.vertices-o).astype('float32')),nthreads=1,nsamples=3).numpy();volume=proof['topology']['oriented_triangle_integral_m3']*1e6;report=dict(passed=proof['passed']and bool(np.all(d<0))and abs(volume-310.5767051825339)<.05,source_sha256=digest_file(b/'surface-covering-v43/core.npz'),simplification_tolerance_m=.00001,source_volume_mL=310.5767051825339,volume_mL=volume,triangles=len(core.faces),core_inside_outer=bool(np.all(d<0)),closest_outer_m=float(-d.max()),failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
