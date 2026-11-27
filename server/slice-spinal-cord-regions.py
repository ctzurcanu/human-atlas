"""Apply root-registered regional planes to the validated shared cord cells."""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file
from geometry_tools.algorithms.plane_tissue_partition import partition
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('--plan',required=True);ap.add_argument('--output',required=True);args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion'
plan_path=Path(args.plan).resolve();plan=json.loads(plan_path.read_text());source=base/'female-regions-v2';z=np.load(source/'cells.npz');d=json.loads((source/'cells.json').read_text())
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
result=partition(z['positions'],{k:z[k]for k in d['surfaces']},d['cells'],plan['bands'],origin=[0.,1.35,-.05])
proof={}
for name,uu in result['cells'].items():
    faces=np.concatenate([result['surfaces'][u['surface']]if u['side']==1 else result['surfaces'][u['surface']][:,[0,2,1]]for u in uu])
    t=topology(Mesh(result['positions'],faces,'model:local-female-ta98:assembled-metres'))
    proof[name]={k:t[k]for k in ['closed_oriented_manifold','open_edges','nonmanifold_edges','invalid_vertex_links','zero_area_faces']}
bad=[n for n,t in proof.items()if not t['closed_oriented_manifold']]
atomic_npz(out/'cells.npz',positions=result['positions'],**result['surfaces'])
atomic_json(out/'cells.json',dict(cells=result['cells'],surfaces={k:len(f)for k,f in result['surfaces'].items()},
    attributes=plan['attributes'],bands=plan['bands'],construction=result['evidence'],topology=proof,
    source_sha256=digest_file(source/'cells.npz'),anatomical_acceptance=False,installed=False,
    passed=not bad,elapsed_seconds=time.monotonic()-start))
assert not bad,bad
print(json.dumps(dict(cells=len(proof),all_closed=True,elapsed_seconds=time.monotonic()-start)),flush=True)
