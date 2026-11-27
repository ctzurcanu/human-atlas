"""Bounded browser-storage rounding of the complete finite PV source walls."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json
from geometry_tools.audit import audit
from geometry_tools.algorithms.bounded_float32_rounding import round_diagnosed_crossings
ap=argparse.ArgumentParser();ap.add_argument('--source',required=True);ap.add_argument('--stage',required=True)
args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart'
source=h/args.source;out=h/args.stage;out.mkdir(exist_ok=False)
base=Mesh.load(h/'blood-spaces-v126/installed-boundary.npz');o=np.array([0.,1.375,.03125])
dtype=np.dtype((np.void,12));basekeys=(base.vertices-o).astype('<f4').view(dtype).ravel();rows=[]
for i in range(4):
    m=Mesh.load(source/f'root-{i}.npz');proof=json.loads((source/f'root-{i}-audit.json').read_text())
    defects=proof['self_intersections'];assert defects['complete_within_pair_budget']
    assert len(defects['examples'])==defects['defect_pairs']
    fixed=np.flatnonzero(np.isin((m.vertices-o).astype('<f4').view(dtype).ravel(),basekeys))
    r=round_diagnosed_crossings(m.vertices,m.faces,defects['examples'],origin=o,
        protected=fixed,seconds=10.,maximum_moves=64,maximum_pair_tests=3000000,
        maximum_displacement_m=2e-8)
    atomic_json(out/f'root-{i}-rounding.json',r['evidence'])
    z=Mesh(r['positions'],m.faces,m.frame);z.save(out/f'root-{i}.npz')
    a=audit(z,require_closed=True,max_pairs=1000000,max_tests=7000000)
    atomic_json(out/f'root-{i}-audit.json',a)
    row=dict(i=i,passed=a['passed'],failures=a['failures'],
        crossings=a['self_intersections']['defect_pairs'],rounding=r['evidence'])
    rows.append(row);print(json.dumps(row),flush=True)
atomic_json(out/'construction.json',dict(passed=all(r['passed']for r in rows),
    roots=rows,elapsed_seconds=time.monotonic()-start,installed=False))
