"""Audit saved source cells in bounded, source-bound batches, then finalize.

Construction and validation have separate deadlines. An incomplete source never
gets passed=True; every final record must match the literal geometry and checker.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

ap=argparse.ArgumentParser();ap.add_argument('source');ap.add_argument('--start',type=int,default=0)
ap.add_argument('--count',type=int,default=0);ap.add_argument('--finalize',action='store_true')
args=ap.parse_args();started=time.monotonic();source=Path(args.source).resolve()
path=source/'cells.json';m=json.loads(path.read_text())
assert m.get('completed_construction_stage')and m.get('source_topology_pending')
assert m['construction']['all_positive_reference_jacobians']
binding=dict(source_sha256=digest_file(source/'cells.npz'),checker_sha256=digest_file(Path(__file__)),
    topology_sha256=digest_file(Path(__file__).parents[1]/'geometry_tools/audit.py'))
folder=source/'topology-checks';folder.mkdir(exist_ok=True);names=sorted(m['cells'])
if args.finalize:
    checks={}
    for p in sorted(folder.glob('batch-*.json')):
        row=json.loads(p.read_text());assert row['binding']==binding and row['passed']
        for name,result in row['checks'].items():
            if name in checks:assert checks[name]==result
            checks[name]=result
    assert set(checks)==set(names),'Every source cell must be independently checked'
    assert all(c['closed_oriented_manifold']and c['oriented_triangle_integral_m3']>0 for c in checks.values())
    m.update(passed=True,source_topology_pending=False,topology=checks,topology_stage_binding=binding)
    atomic_json(path,m,overwrite=True)
    print(json.dumps(dict(passed=True,complete_cells=len(checks),elapsed_seconds=time.monotonic()-started)),flush=True)
else:
    assert args.count>0 and args.start>=0
    with np.load(source/'cells.npz')as z:arrays={k:z[k]for k in z.files}
    points=arrays['positions'];checks={}
    for name in names[args.start:args.start+args.count]:
        faces=np.concatenate([arrays[u['surface']]if u['side']==1 else arrays[u['surface']][:,[0,2,1]]for u in m['cells'][name]])
        used,inv=np.unique(faces,return_inverse=True)
        check=topology(Mesh(points[used],inv.reshape(-1,3),'female-shared-regional-cord'))
        checks[name]={k:check[k]for k in ['closed_oriented_manifold','open_edges','nonmanifold_edges','zero_area_faces','components','oriented_triangle_integral_m3']}
    assert checks
    passed=all(c['closed_oriented_manifold']and c['oriented_triangle_integral_m3']>0 for c in checks.values())
    atomic_json(folder/f'batch-{args.start}.json',dict(binding=binding,checks=checks,passed=passed,
        elapsed_seconds=time.monotonic()-started))
    print(json.dumps(dict(passed=passed,checked=len(checks),total=len(names),elapsed_seconds=time.monotonic()-started)),flush=True)
    assert passed,[n for n,c in checks.items()if not c['closed_oriented_manifold']or c['oriented_triangle_integral_m3']<=0]
