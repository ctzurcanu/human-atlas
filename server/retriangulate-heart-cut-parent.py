"""Re-triangulate one diagnosed source chart without moving its boundary."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
import manifold3d as md
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH,box_distance_squared
from geometry_tools.distance import triangle_pairs
from geometry_tools.algorithms.tube_continuity import boundary_loops
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';src=h/args.source;out=h/args.stage;out.mkdir(exist_ok=False)
meta=json.loads((src/'arrangement.json').read_text());z=np.load(src/'arrangement.npz');p=z['positions'];arrays={k:z[k]for k in z.files}
rows=sorted(meta['surfaces'],key=lambda r:r['id']);f=np.concatenate([z[r['array']]for r in rows]);diag=json.loads((src/'diagnostics.json').read_text());assert diag['geometry_sha256']==Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha
sid='interface:8:-1';parent=4002;r=next(r for r in rows if r['id']==sid);local=np.flatnonzero(z[r['parents']]==parent);ff=z[r['array']];old_patch=ff[local];loops=boundary_loops(old_patch);assert len(loops)==1
loop=loops[0];assert len(np.unique(old_patch))==len(loop),'This repair accepts boundary-only simple charts'
old=h/'regional-bank-v102';e=json.loads((old/'emission.json').read_text());oz=np.load(old/'arrangement.npz');origin=np.array(meta['origin']);op=(oz['positions']-origin).astype('<f4').astype(float)+origin;q=op[oz[e['surfaces'][sid]][parent]];n=np.cross(q[1]-q[0],q[2]-q[0]);drop=np.argmax(abs(n));axes=np.delete(np.arange(3),drop);uv=p[loop][:,axes];uv-=uv.mean(0)
tri=np.asarray(md.triangulate([uv],epsilon=0.,allow_convex=False),int);assert len(tri)==len(old_patch)
replacement=loop[tri];nn=np.cross(p[replacement[:,1]]-p[replacement[:,0]],p[replacement[:,2]]-p[replacement[:,0]])
if np.all(nn@n<0):replacement=replacement[:,[0,2,1]];nn=-nn
print('Retriangulated source-parent signs',np.sign(nn@n).astype(int).tolist(),flush=True)
assert np.all(nn@n>0),'No reversed or zero source-parent facets permitted'
updated=ff.copy();updated[local]=replacement;arrays[r['array']]=updated
offset=sum(len(z[row['array']])for row in rows[:rows.index(r)]);ids=local+offset;newf=f.copy();newf[ids]=replacement;tree=BVH(Mesh(p,f,'source'))
pairs=set()
for face in ids:
    points=p[newf[face]];near=tree.query_bounds(points.min(0),points.max(0),distance_m=0.)
    pairs.update((min(int(face),int(j)),max(int(face),int(j)))for j in near if j!=face and j not in ids)
pairs.update((int(a),int(b))for a in ids for b in ids if a<b)
pairs=np.array(sorted(pairs),int).reshape(-1,2);tt=p[newf[pairs]];overlap=box_distance_squared(tt[:,0].min(1),tt[:,0].max(1),tt[:,1].min(1),tt[:,1].max(1))<=0;bad=np.zeros(len(pairs),bool);bad[overlap]=np.isin(triangle_pairs(tt[overlap,0],tt[overlap,1])['relation'],['crossing','coplanar-overlap'])
remaining={tuple(v)for v in diag['remaining_global_defect_pairs']if not np.isin(v,ids).any()}|set(map(tuple,pairs[bad]))
atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',meta)
atomic_json(out/'diagnostics.json',dict(geometry_sha256=Mesh(p,newf,'model:local-ta98:assembled-metres').geometry_sha,remaining_global_defect_pairs=sorted(remaining),complete_changed_triangle_cover=True))
report=dict(passed=not remaining,source_parent=parent,surface=sid,old_faces=len(local),new_faces=len(replacement),boundary_controls_literal=True,source_parent_orientations_positive=True,tested_pairs=int(overlap.sum()),remaining_global_crossings=len(remaining),elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',report);print(json.dumps(report),flush=True)
