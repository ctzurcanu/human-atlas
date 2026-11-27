"""Contract only diagnosed sub-voxel numerical fins on the current wall."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.audit import audit
from geometry_tools.algorithms.shared_edge_collapse import collapse_short_edges
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/heart';source=base/'whole-wall-v17';out=base/'whole-wall-v18';out.mkdir(exist_ok=True);assert not any(out.iterdir());m=Mesh.load(source/'whole-wall.npz');d=Mesh.load(source/'double-wall.npz');a=json.loads((source/'emitted-audit.json').read_text());assert a['self_intersections']['defect_pairs']==len(a['self_intersections']['examples']);q=d.triangles;normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];words=np.broadcast_to(normal[:,None,:],(len(q),3,3)).copy()
cells={'wall':[dict(surface='wall',side=1,neighbor=None)]};repair=collapse_short_edges(m.vertices,{'wall':m.faces},cells,{'wall':words},desired_positions=d.vertices,maximum_displacement_m=.0006,maximum_steps=12,wall_seconds=12.,candidate_face_rings=1,defect_face_pairs=a['self_intersections']['examples']);surfaces=repair['surfaces'];report=repair['evidence'];used,inv=np.unique(surfaces['wall'],return_inverse=True);mesh=Mesh(m.vertices[used],inv.reshape(-1,3),m.frame);mesh.save(out/'whole-wall.npz');atomic_json(out/'contractions.json',report);proof=audit(mesh,require_closed=True,expected_components=1,max_pairs=2000000,max_tests=16000000);atomic_json(out/'emitted-audit.json',proof);result=dict(passed=proof['passed'],source_sha256=digest_file(source/'whole-wall.npz'),maximum_contraction_bound_m=.0006,triangles=len(mesh.faces),vertices=len(mesh.vertices),remaining_intersections=proof['self_intersections']['defect_pairs'],failures=proof['failures'],elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,installed=False);atomic_json(out/'construction.json',result);print(json.dumps(result),flush=True)
