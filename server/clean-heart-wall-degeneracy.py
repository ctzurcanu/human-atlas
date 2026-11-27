"""MeshFix cleanup of a complete closed wall, never open anatomical mouths.

The input has no boundary edges: holes created while removing degenerate
triangles are numerical patch boundaries, not the vascular passages. Actual
vessel sections and both blood circuits remain separate acceptance checks.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import pymeshfix
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'whole-wall-v15/whole-wall.npz';out=b/'whole-wall-v16';out.mkdir(exist_ok=False);m=Mesh.load(source);o=np.array([0.,1.375,.03125]);before=json.loads((source.parent/'emitted-audit.json').read_text());assert before['topology']['open_edges']==0
tool=pymeshfix.PyTMesh();tool.set_quiet(True);tool.load_array(m.vertices-o,m.faces.astype('int32'));print('Repair connected wall',tool.n_points,tool.n_faces,flush=True);complete=tool.clean(max_iters=3,inner_loops=2);v,f=tool.return_arrays();v=(v.astype('float32').astype(float)+o);mesh=Mesh(v,f.astype('int64'),m.frame);mesh.save(out/'whole-wall.npz');print('Repaired contour',len(f),'faces',round(time.monotonic()-start,1),'s',flush=True)
scene=o3d.t.geometry.RaycastingScene(nthreads=1);scene.add_triangles(o3d.core.Tensor((m.vertices-o).astype('float32')),o3d.core.Tensor(m.faces.astype('uint32')));distance=scene.compute_distance(o3d.core.Tensor((v-o).astype('float32')),nthreads=1).numpy();proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof);report=dict(passed=proof['passed']and float(distance.max())<.001,meshfix_completed=bool(complete),source_sha256=digest_file(source),source_boundary_edges=0,components_joined=False,source_vertices=len(m.vertices),vertices=len(v),triangles=len(f),source_deviation_quantiles_m=np.quantile(distance,[0,.5,.9,.99,1]).tolist(),maximum_cleanup_displacement_limit_m=.001,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
