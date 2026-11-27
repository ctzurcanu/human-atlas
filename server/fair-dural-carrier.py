"""Fair only failed coat patches within an explicit 0.12 mm source bound."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops

ap=argparse.ArgumentParser();ap.add_argument('carrier');ap.add_argument('--fields',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();source=Path(args.carrier).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
z=np.load(source/'surface.npz');v=z['vertices'].copy();f=z['faces'];original=v.copy();used,inv=np.unique(f,return_inverse=True)
proof=json.loads((Path(args.fields)/'field-check.json').read_text());bad=np.unique(np.concatenate([np.asarray(x,int)for x in proof['prisms']['invalid_faces_by_interval']]))
affected=set(map(int,f[bad].ravel()));adj=[set()for _ in v]
for tri in f:
    for a,b in zip(tri,np.roll(tri,-1)):adj[int(a)].add(int(b));adj[int(b)].add(int(a))
for _ in range(2):affected.update(j for i in list(affected)for j in adj[i])
fixed=set(map(int,np.concatenate(boundary_loops(f))));fixed.add(int(np.argmin(v[:,1])))
active=np.array(sorted(affected-fixed),int);limit=.00012
atomic_json(out/'plan.json',dict(source_sha256=digest_file(source/'surface.npz'),
    reason='Local source cuffs and the caudal connecting strip cannot support two finite coats without surface fairing',
    maximum_displacement_m=limit,active_vertices=active.tolist(),fixed_cranial_rim=sorted(fixed),
    every_skeleton_and_nerve_coordinate_fixed=True,anatomical_acceptance=False,installed=False))
for iteration in range(8):
    for weight in [.35,-.36]:
        delta=np.array([v[list(adj[i])].mean(0)-v[i]for i in active]);candidate=v[active]+weight*delta
        movement=candidate-original[active];length=np.linalg.norm(movement,axis=1)
        movement*=np.minimum(1,limit/np.maximum(length,1e-100))[:,None];v[active]=original[active]+movement
mesh=Mesh(v,f,'model:local-female-ta98:assembled-metres');check=audit(mesh,require_closed=False,max_pairs=250000,max_tests=2000000)
atomic_npz(out/'surface.npz',vertices=v,faces=f,frame=np.array(mesh.frame))
metadata=json.loads((source/'surface.json').read_text());metadata.update(source_carrier_sha256=digest_file(source/'surface.npz'),
    original_upper_positions_literal=False,unaffected_controls_literal=True,modified_vertices=len(active),
    maximum_surface_control_displacement_m=float(np.linalg.norm(v-original,axis=1).max()),
    audit=check,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start)
atomic_json(out/'surface.json',metadata);assert check['passed'],check['failures']
print(json.dumps(dict(vertices_faired=len(active),maximum_displacement_mm=float(np.linalg.norm(v-original,axis=1).max()*1000),
    self_crossings=check['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start)),flush=True)
