"""New bounded candidate from a fully saved, dense contour checkpoint.

The previous audit timed out. Its candidate is not installed or accepted.
This explicitly changes the geometry under a .4 mm contour tolerance, then
runs all geometry gates again on the new, complete browser representation.
"""
from pathlib import Path
import json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';source=base/'whole-wall-v14/whole-wall.npz';out=base/'whole-wall-v15';out.mkdir(exist_ok=False);m=Mesh.load(source);o=np.array([0.,1.375,.03125]);solid=md.Manifold(md.Mesh64(m.vertices-o,m.faces.astype('uint64')));assert str(solid.status())=='Error.NoError';simple=solid.simplify(.0004);z=simple.to_mesh64();v=z.vert_properties[:,:3];f=z.tri_verts.astype('int64');assert np.array_equal(v,v.astype('float32').astype(float));mesh=Mesh(v+o,f,m.frame);mesh.save(out/'whole-wall.npz');print('Compact candidate',len(f),'faces',round(time.monotonic()-start,1),'s',flush=True);proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof);report=dict(passed=proof['passed'],source_checkpoint_sha256=digest_file(source),source_checkpoint_audit='timed-out; not accepted or installed',new_candidate_geometry_sha256=mesh.geometry_sha,simplification_tolerance_m=.0004,vertices=len(v),triangles=len(f),volume_mL=simple.volume()*1e6,genus=int(simple.genus()),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
