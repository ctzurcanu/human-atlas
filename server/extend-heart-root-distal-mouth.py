"""Continue an already checked paired root to a clear distal review mouth.

Existing controls, cusp hinges and source wall correspondences remain literal.
The added end rings translate together, preserving their section and thickness.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops, annular_rim, finite_wall_cells

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True)
ap.add_argument('--millimeters', required=True, type=float)
ap.add_argument('--stage', required=True)
ap.add_argument('--direction', nargs=3, type=float,
                help='Explicit anatomical continuation direction in the source body frame')
args = ap.parse_args()
assert 0 < args.millimeters <= 10
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
source, output = h / args.input, h / args.stage
output.mkdir(exist_ok=False)
meta = json.loads((source / 'construction.json').read_text())
assert meta['passed'] and meta.get('root_continuation')
z = np.load(source / 'source-walls.npz')
outer, inner, faces, levels = [z[k].copy() for k in
    ['outer_positions', 'inner_positions', 'source_faces', 'station_indices']]
ring = max(boundary_loops(faces), key=lambda r: outer[r].mean(0)[1])
level = int(levels[ring[0]])
neighbours = np.unique(faces[np.any(np.isin(faces, ring), axis=1)])
previous = np.setdiff1d(neighbours, ring)
direction = outer[ring].mean(0) - outer[previous].mean(0)
if args.direction:
    direction = np.asarray(args.direction, dtype=float)
direction /= np.linalg.norm(direction)
delta = direction * args.millimeters / 1000
old_count = len(outer)
old_outer, old_inner = outer[ring].copy(), inner[ring].copy()
for i, t in enumerate([.5, 1.]):
    next_ring = np.arange(len(outer), len(outer) + len(ring))
    faces = np.r_[faces, annular_rim(ring, next_ring)]
    outer = np.r_[outer, old_outer + t * delta]
    inner = np.r_[inner, old_inner + t * delta]
    levels = np.r_[levels, np.full(len(ring), int(levels.max()) + 1)]
    ring = next_ring
full = [faces, faces[:, [0, 2, 1]] + len(outer)]
for ring in boundary_loops(faces):
    full.append(annular_rim(ring, ring + len(outer)))
mesh = Mesh(np.r_[outer, inner], np.concatenate(full), 'model:local-ta98:assembled-metres')
wall = finite_wall_cells(np.stack([outer, inner]), faces)
check = audit(mesh, require_closed=True, max_pairs=500000, max_tests=4000000)
passed = wall['finite_wall_cells_pass'] and check['passed']
mesh.save(output / ('finite-root.npz' if passed else 'diagnostic-root.npz'))
atomic_npz(output / 'source-walls.npz', outer_positions=outer, inner_positions=inner,
    source_faces=faces, station_indices=levels)
meta.update(passed=bool(passed), source_stage=args.input, distal_extension_m=delta.tolist(),
    prior_controls_literal=bool(np.array_equal(outer[:old_count], z['outer_positions']) and
        np.array_equal(inner[:old_count], z['inner_positions'])),
    wall=wall, audit=check, elapsed_seconds=time.monotonic()-start, installed=False)
atomic_json(output / 'construction.json', meta)
print(json.dumps(dict(passed=bool(passed), failures=check['failures'],
    finite_wall=wall['finite_wall_cells_pass'], extension_mm=args.millimeters,
    elapsed_seconds=time.monotonic()-start)), flush=True)
raise SystemExit(0 if passed else 2)
