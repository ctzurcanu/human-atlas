"""Remove four explicitly diagnosed source slivers before finite coating.

Only the two 43 um join edges and two four-face source cavities change.
This is a bounded geometry proposal, not anatomical acceptance.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();source=Path(args.source).resolve();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
z=np.load(source/'surface.npz');raw=z['vertices'];used,inv=np.unique(z['faces'],return_inverse=True);v=raw[used];f=inv.reshape(-1,3);original=f.copy();edits=[]
def quality(ff):
    q=v[ff];area=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)
    return 2*np.sqrt(3)*area/np.sum((q-np.roll(q,1,axis=1))**2,axis=(1,2))
for remove,keep in [(4860,4878),(4890,4903)]:
    assert np.linalg.norm(v[remove]-v[keep])<.00005
    selected=np.any(f==remove,axis=1);before=f[selected].copy();f[f==remove]=keep
    retain=np.all(np.diff(np.sort(f,axis=1),axis=1)!=0,axis=1);f=f[retain]
    edits.append(dict(kind='contract-diagnosed-join-edge',removed_control=remove,kept_control=keep,
        displacement_bound_m=float(np.linalg.norm(v[remove]-v[keep])),source_faces=before.tolist()))
for vertex in [466,1212]:
    selected=np.any(f==vertex,axis=1);old=f[selected].copy();loops=boundary_loops(old)
    assert len(old)==4 and len(loops)==1 and len(loops[0])==4
    ring=loops[0];candidates=[ring[np.array([[0,1,2],[0,2,3]])],ring[np.array([[1,2,3],[1,3,0]])]]
    q=v[old];normal=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]).sum(0)
    candidates=[q for q in candidates if np.all(np.cross(v[q[:,1]]-v[q[:,0]],v[q[:,2]]-v[q[:,0]])@normal>0)]
    assert candidates,'No orientation-preserving cavity triangulation'
    chosen=max(candidates,key=lambda q:float(quality(q).min()))
    # A cavity is replaced only within its existing four-control hull. The
    # coarse hull diameter is a rigorous bound; sampled distances below
    # provide the finer diagnostic without claiming a Hausdorff certificate.
    bound=float(np.linalg.norm(v[ring][:,None]-v[ring][None,:],axis=2).max())
    assert bound<.005
    f=np.concatenate([f[~selected],chosen])
    edits.append(dict(kind='remove-four-face-source-sliver',removed_control=vertex,rim=ring.tolist(),
        before=old.tolist(),after=chosen.tolist(),patch_hull_diameter_m=bound))
mesh=Mesh(v,f,'model:local-female-ta98:assembled-metres');proof=audit(mesh,require_closed=False,max_pairs=250000,max_tests=2000000)
atomic_npz(out/'surface.npz',vertices=v,faces=f,frame=np.array(mesh.frame))
record=json.loads((source/'surface.json').read_text());record.update(source_sha256=digest_file(source/'surface.npz'),
    local_carrier_repairs=edits,retained_control_positions_literal=True,audit=proof,anatomical_acceptance=False,installed=False,
    elapsed_seconds=time.monotonic()-start)
atomic_json(out/'surface.json',record)
assert proof['passed'],proof['failures']
assert len(boundary_loops(f))==1 and len(boundary_loops(f)[0])==28
print(json.dumps(dict(passed=True,faces_before=len(original),faces_after=len(f),edits=len(edits),
    self_crossings=proof['self_intersections']['defect_pairs'],elapsed_seconds=time.monotonic()-start)),flush=True)
