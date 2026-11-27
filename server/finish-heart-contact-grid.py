"""Bounded local contact repair, with a smaller bound for retained source points."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.bounded_float32_rounding import round_diagnosed_crossings
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);source=h/'pv-root-conditioned-complex-v170';bank=h/'pv-root-bank-v159'
meta=json.loads((source/'arrangement.json').read_text());z=np.load(source/'arrangement.npz');p=z['positions'];o=np.array(meta['origin']);anchor=np.load(bank/'arrangement.npz')['positions']
rows={r['id']:r for r in meta['surfaces']};f=np.concatenate([z[rows[sid]['array']]for sid in sorted(rows)]);proof=json.loads((source/'diagnostics.json').read_text());assert proof['geometry_sha256']==Mesh(p,f,'model:local-ta98:assembled-metres').geometry_sha;pairs=proof['remaining_global_defect_pairs']
base=np.load(h/'regional-bank-v102/arrangement.npz')['positions'];roots=[Mesh.load(h/f'pv-root-conditioned-source-v156/root-{i}.npz').vertices for i in range(4)];typ=np.dtype((np.void,12));source_mask=np.isin((anchor-o).astype('<f4').view(typ).ravel(),(np.r_[base,*roots]-o).astype('<f4').view(typ).ravel())
guard=LocalCrossingGuard(p,f,pairs,maximum_displacement_m=2e-8,maximum_pair_tests=4000000)
selected=np.unique(f[np.array(pairs,int)]);incident=np.unique(np.concatenate([guard.incident[guard.starts[v]:guard.starts[v+1]]for v in selected]));local=np.unique(f[incident]);allowed=local[source_mask[local]];assert len(allowed)<=512
fixed=np.flatnonzero(source_mask&~np.isin(np.arange(len(p)),allowed));limits=np.full(len(p),2e-8);limits[allowed]=1e-8
print('Final contacts',len(pairs),'; retained source controls eligible',len(allowed),flush=True)
r=round_diagnosed_crossings(p,f,pairs,origin=o,protected=fixed,seconds=20.,maximum_moves=96,
    maximum_pair_tests=4000000,maximum_displacement_m=2e-8,anchor_positions=anchor,
    maximum_displacement_by_vertex=limits)
atomic_json(out/'rounding.json',r['evidence']);arrays={k:z[k]for k in z.files};arrays['positions']=r['positions'];atomic_npz(out/'arrangement.npz',**arrays);atomic_json(out/'arrangement.json',meta)
changed=np.flatnonzero(np.any(r['positions']!=anchor,axis=1)&source_mask)
proof=dict(passed=r['evidence']['passed'],remaining_global_crossings=r['evidence']['remaining_global_defect_pair_count'],
    retained_source_controls_changed=len(changed),changed_source_ids=changed.tolist(),
    maximum_source_adjustment_m=float(np.linalg.norm(r['positions'][changed]-anchor[changed],axis=1).max(initial=0)),
    maximum_total_adjustment_m=float(np.linalg.norm(r['positions']-anchor,axis=1).max()),
    stop=r['evidence']['stop_reason'],elapsed_seconds=time.monotonic()-start,installed=False)
atomic_json(out/'construction.json',proof);atomic_json(out/'diagnostics.json',dict(geometry_sha256=Mesh(r['positions'],f,'model:local-ta98:assembled-metres').geometry_sha,remaining_global_defect_pairs=r['evidence']['remaining_global_defect_pairs'],complete_guarded_worklist=True));print(json.dumps(proof),flush=True)
