"""Local numerical cleanup of closed venous Boolean output, with bounds."""
from pathlib import Path
import sys,time,json
import numpy as np
import pymeshfix,open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'native-pv-mouths-v67';out.mkdir(exist_ok=False);source=b/'native-pv-mouths-v66/whole-wall.npz';m=Mesh.load(source);o=np.array([0.,1.375,.03125]);tool=pymeshfix.PyTMesh();tool.set_quiet(True);tool.load_array((m.vertices-o)*1000,m.faces.astype('int32'));tool.clean(max_iters=3,inner_loops=2);v,f=tool.return_arrays();v=(v/1000).astype('float32').astype(float)+o;wall=Mesh(v,f.astype('int64'),m.frame);wall.save(out/'whole-wall.npz');proof=audit(wall,require_closed=True,max_pairs=3500000,max_tests=26000000);atomic_json(out/'emitted-audit.json',proof);distances=[]
for first,second in [(m,wall),(wall,m)]:
 scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((second.vertices-o).astype('float32')),o3d.core.Tensor(second.faces.astype('uint32')));p=np.r_[first.vertices,first.triangles.mean(1)];d=scene.compute_distance(o3d.core.Tensor((p-o).astype('float32')),nthreads=1).numpy();distances.append(float(d.max()))
before=json.loads((b/'native-pv-mouths-v66/emitted-audit.json').read_text())['topology']['oriented_triangle_integral_m3'];after=proof['topology']['oriented_triangle_integral_m3'];error=abs(after-before)*1e6;passed=proof['passed']and max(distances)<.00002 and error<.01;report=dict(passed=passed,source_sha256=digest_file(source),geometry_sha256=wall.geometry_sha,maximum_two_way_surface_sample_error_m=distances,maximum_allowed_error_m=.00002,volume_change_mL=error,maximum_allowed_volume_change_mL=.01,triangles=len(wall.faces),failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
