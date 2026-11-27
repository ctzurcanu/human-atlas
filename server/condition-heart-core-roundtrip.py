"""Condition the exact eroded core at its declared storage precision."""
from pathlib import Path
import sys,time,json
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';out=b/'surface-covering-v46';out.mkdir(exist_ok=False)
m=Mesh.load(b/'surface-covering-v43/core.npz');o=np.array([0.,1.375,.03125]);v=(m.vertices-o).astype('float32').astype(float)
raw=md.Mesh64(v,m.faces.astype('uint64'),tolerance=2e-7);raw.merge();solid=md.Manifold(raw)
assert str(solid.status())=='Error.NoError',str(solid.status())
solid=solid.set_tolerance(2e-7).simplify(2e-7);z=solid.to_mesh64();r=normalize_shared_boundary(z.vert_properties[:,:3]+o,{'core':z.tri_verts},o)
core=Mesh(r['positions'],r['cells']['core'],m.frame);core.save(out/'core.npz');atomic_json(out/'encoding.json',r['evidence'])
print('Conditioned core',len(core.faces),'triangles',round(time.monotonic()-start,2),'s',flush=True)
proof=audit(core,require_closed=True,max_pairs=5000000,max_tests=30000000);atomic_json(out/'core-audit.json',proof)
report=dict(passed=proof['passed'],source_sha256=digest_file(b/'surface-covering-v43/core.npz'),source_volume_mL=310.5767051825339,volume_mL=proof['topology']['oriented_triangle_integral_m3']*1e6,storage_conditioning_tolerance_m=2e-7,triangles=len(core.faces),failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False)
atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
