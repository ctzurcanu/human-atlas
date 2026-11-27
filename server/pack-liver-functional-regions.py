"""Audit the actual Float32 eight-cell partition and reciprocal interfaces."""
from pathlib import Path
from collections import defaultdict
import sys,json,time,argparse
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.tissue import pack
from geometry_tools.algorithms.quantized_surface_topology import conform_quantized_surfaces

ap=argparse.ArgumentParser();ap.add_argument('--version',default='v144');args=ap.parse_args()
start=time.monotonic();base=ROOT/'.local-models/geometry-workbench/liver-completion';out=base/('functional-territories-'+args.version)
z=np.load(out/'partition.npz');labels=[k for k in z.files if k!='positions'];origin=np.array([0.,1.25,0.])
mapped,ancestry,conform=conform_quantized_surfaces(z['positions'],{label:z[label]for label in labels},origin)
owners=defaultdict(list)
for label in labels:
    for f in mapped[label]:owners[tuple(sorted(map(int,f)))].append((label,f))
patches=defaultdict(list)
for uses in owners.values():
    assert len(uses)<=2
    if len(uses)==2:
        (a,f),(b,g)=uses;assert any(np.array_equal(np.roll(f,k),g[[0,2,1]])for k in range(3))
        patches[(a,b)].append(f)
    else:
        a,f=uses[0];patches[(a,None)].append(f)
surfaces={};cells={label:[]for label in labels}
for i,((a,b),ff)in enumerate(patches.items()):
    sid='functional:'+str(i);surfaces[sid]=np.array(ff,int)
    cells[a].append(dict(surface=sid,side=1,neighbor=b))
    if b:cells[b].append(dict(surface=sid,side=-1,neighbor=a))
ids={label:'TA98LIVER:male:segment:'+label for label in labels}
result=pack(z['positions'],surfaces,cells,ids,'model:local-ta98:assembled-metres',out/'packed',origin=origin,
            max_pairs=2000000,max_tests=16000000,checkpoint_directory=out/'checkpoints')
atomic_json(out/'encoding-evidence.json',dict(passed=result['passed'],conform=conform,
    bytes=result['total_bytes'],elapsed_seconds=time.monotonic()-start,
    method='Actual native eight-segment guide fields partition the repaired whole; the eight territories share reciprocal faces and one position/normal bank.',
    remaining_limits=['Functional boundaries remain anatomical approximations of the supplied native segment guides.',
      'Hilar bulk closures are calculation boundaries, not membranes. Conduit envelopes are subtracted during section rendering.']))
print(json.dumps(dict(passed=result['passed'],failures=result['checks']['failures'],bytes=result['total_bytes'],elapsed_seconds=time.monotonic()-start)),flush=True)
