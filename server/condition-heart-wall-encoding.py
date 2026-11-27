"""Remove sub-micrometre contour slivers before actual Float32 audits."""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'whole-wall-v9/whole-wall.npz';out=b/'whole-wall-v10';out.mkdir(exist_ok=False);m=Mesh.load(source);o=np.array([0.,1.38,.03]);solid=md.Manifold(md.Mesh64(m.vertices-o,m.faces.astype('uint64')));assert str(solid.status())=='Error.NoError';result=solid.set_tolerance(.000001);z=result.to_mesh64();v=z.vert_properties[:,:3];f=z.tri_verts.astype('int64');encoded=Mesh(v.astype('float32').astype(float)+o,f,m.frame);encoded.save(out/'whole-wall.npz');proof=audit(encoded,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof)
report=dict(passed=proof['passed'],source_sha256=digest_file(source),numerical_tolerance_m=.000001,vertices=len(v),triangles=len(f),volume_mL=result.volume()*1e6,genus=int(result.genus()),maximum_rounding_error_m=float(np.linalg.norm(v.astype('float32').astype(float)-v,axis=1).max()),elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
