"""Conditioned intersection cleanup; preserve both double and GPU evidence."""
from pathlib import Path
import json,sys,time
import numpy as np
import pymeshfix
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
start=time.monotonic();b=ROOT/'.local-models/geometry-workbench/heart';source=b/'whole-wall-v16/whole-wall.npz';out=b/'whole-wall-v17';out.mkdir(exist_ok=False);m=Mesh.load(source);o=np.array([0.,1.375,.03125]);tool=pymeshfix.PyTMesh();tool.set_quiet(True);tool.load_array((m.vertices-o)*1000,m.faces.astype('int32'));changed=tool.clean(max_iters=10,inner_loops=3);v,f=tool.return_arrays();v=v/1000;f=f.astype('int64');double=Mesh(v+o,f,m.frame);double.save(out/'double-wall.npz');encoded=Mesh(v.astype('float32').astype(float)+o,f,m.frame);encoded.save(out/'whole-wall.npz');print('Cleaned wall',len(f),'faces',round(time.monotonic()-start,1),'s',flush=True)
checks={}
for name,mesh in [('double',double),('emitted',encoded)]:
 proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);checks[name]=dict(passed=proof['passed'],failures=proof['failures'],self_crossings=proof['self_intersections']['defect_pairs']);atomic_json(out/(name+'-audit.json'),proof);print(name,json.dumps(checks[name]),flush=True)
report=dict(passed=all(x['passed']for x in checks.values()),meshfix_changed=bool(changed),checks=checks,source_sha256=digest_file(source),coordinate_conditioning_scale=1000,vertices=len(v),triangles=len(f),elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
