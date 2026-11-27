"""Bounded Taubin fairing of the repaired wall, with full emitted checks."""
from pathlib import Path
import json,sys,time
import numpy as np
import open3d as o3d
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';p=b/'whole-wall-v18/whole-wall.npz';assert json.loads((p.parent/'construction.json').read_text())['passed'];out=b/'whole-wall-v22';out.mkdir(exist_ok=False);m=Mesh.load(p);o=np.array([0,1.375,.03125]);source=o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(m.vertices-o),o3d.utility.Vector3iVector(m.faces.astype('int32')));smooth=source.filter_smooth_taubin(number_of_iterations=5,lambda_filter=.5,mu=-.53);delta=np.asarray(smooth.vertices)-(m.vertices-o);length=np.linalg.norm(delta,axis=1);delta*=np.minimum(1.,.0006/np.maximum(length,1e-30))[:,None];v=(m.vertices-o+delta).astype('float32').astype(float)+o;displacement=np.linalg.norm(v-m.vertices,axis=1);assert displacement.max()<.00075,'Fairing exceeded explicit source control bound';mesh=Mesh(v,m.faces,m.frame);mesh.save(out/'whole-wall.npz');proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof);result=dict(passed=proof['passed'],source_sha256=digest_file(p),method='Open3D Taubin, per-control bounded displacement',iterations=5,lambda_filter=.5,mu=-.53,maximum_displacement_limit_m=.00075,displacement_quantiles_m=np.quantile(displacement,[0,.5,.9,.99,1]).tolist(),triangles=len(m.faces),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
