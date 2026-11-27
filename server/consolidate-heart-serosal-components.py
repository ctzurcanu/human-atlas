"""Keep the continuous outer stratum; reassign isolated micro-islands inward.

No tissue is deleted. Opposed interfaces cancel literally and the union stays
the complete original coat. Unsupported small attachment allocations are
reported rather than exposed as floating fragments.
"""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh

ap = argparse.ArgumentParser()
ap.add_argument('--input', required=True)
ap.add_argument('--stage', required=True)
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
src, out = h / args.input, h / args.stage
out.mkdir(exist_ok=False)
z = np.load(src / 'cells.npz')
p, a, b = z['positions'], z['serosa_faces'], z['subserosa_faces']
edges = np.r_[a[:, [0, 1]], a[:, [1, 2]], a[:, [2, 0]]]
graph = coo_matrix((np.ones(len(edges)), (edges[:, 0], edges[:, 1])), shape=(len(p), len(p)))
_, labels = connected_components(graph, directed=False)
face_labels = labels[a[:, 0]]
tri = p[a]
area = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)/2
areas = np.bincount(face_labels, weights=area)
main = int(areas.argmax())
retained = face_labels == main
removed = a[~retained]
fraction = float(area[~retained].sum()/area.sum())
assert fraction < .01, 'More than 1% of the serosal candidate is disconnected'
combined = np.r_[b, removed]
canonical = np.sort(combined, axis=1)
_, first, inverse, count = np.unique(canonical, axis=0, return_index=True,
                                     return_inverse=True, return_counts=True)
assert count.max() <= 2
for i in np.flatnonzero(count == 2):
    pair = combined[np.flatnonzero(inverse == i)]
    assert any(np.array_equal(np.roll(pair[0], k), pair[1][[0, 2, 1]]) for k in range(3))
inner = combined[first[count == 1]]
atomic_npz(out / 'cells.npz', positions=p, serosa_faces=a[retained], subserosa_faces=inner)
for name, f in [('serosa', a[retained]), ('subserosa', inner)]:
    used, inv = np.unique(f, return_inverse=True)
    Mesh(p[used], inv.reshape(-1, 3), 'cardiac-canonical-bank-metres').save(out / (name+'.npz'))
meta = json.loads((src / 'construction.json').read_text())
meta.update(continuous_serosa_only=True, source_partition_sha256=digest_file(src / 'cells.npz'),
            reassigned_disconnected_component_count=int(len(np.unique(face_labels))-1),
            reassigned_surface_fraction=fraction,
            removed_tissue_volume_m3=0, no_tissue_deleted=True, source_union_literal=True,
            cells=[dict(name=name, faces=len(f)) for name, f in [('serosa', a[retained]), ('subserosa', inner)]],
            elapsed_seconds=time.monotonic()-start)
atomic_json(out / 'construction.json', meta)
print('Continuous serosa retained;', len(np.unique(face_labels))-1,
      'isolated allocations returned to underlying tissue; surface fraction', fraction,
      'elapsed', time.monotonic()-start, flush=True)
