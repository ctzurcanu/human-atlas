"""Quick actual-cell topology checkpoint before an expensive emitted audit."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json, digest_file

ap = argparse.ArgumentParser()
ap.add_argument('--stage', required=True)
ap.add_argument('--encoded', action='store_true')
args = ap.parse_args()
start = time.monotonic()
d = ROOT / '.local-models/geometry-workbench/heart' / args.stage
z = np.load(d / 'cells.npz')
v = z['positions']
if args.encoded:
    v = v.astype('<f4').astype(float)
rows = {}
for name in ['serosa', 'subserosa']:
    f = z[name+'_faces']
    used, inv = np.unique(f, return_inverse=True)
    row = topology(Mesh(v[used], inv.reshape(-1, 3), 'cardiac-canonical-bank-metres'))
    rows[name] = row
    print(name, 'closed', row['closed_oriented_manifold'], 'nonmanifold', row['nonmanifold_edges'],
          'zero', row['zero_area_faces'], 'components', len(row['components']), flush=True)
atomic_json(d / ('encoded-topology.json' if args.encoded else 'source-topology.json'),
            dict(source_sha256=digest_file(d / 'cells.npz'), encoded=args.encoded, cells=rows,
                 elapsed_seconds=time.monotonic()-start, anatomy_accepted=False))
