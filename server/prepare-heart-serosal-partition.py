"""Prepare a real serosal/subserosal split inside the installed visceral coat.

The calculation mask retreats only the exposed epicardial surface. It never
expands the heart, fills a blood cavity or copies reference-model geometry.
This is a candidate construction; anatomy and emitted geometry require review.
"""
from pathlib import Path
import argparse, json, sys, time
from collections import defaultdict
import numpy as np
from scipy.optimize import linprog

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

ap = argparse.ArgumentParser()
ap.add_argument('--stage', required=True)
ap.add_argument('--bank', default='annular-finished-bank-v378')
ap.add_argument('--serosal-depth-m', type=float, default=.000010)
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
e = json.loads((h / args.bank / 'emission.json').read_text())
assert digest_file(h / args.bank / 'cells.bin') == e['bank_sha256']
src = h / e['source_stage']
meta = json.loads((src / 'ownership.json').read_text())
z = np.load(src / 'arrangement.npz')
p = z['positions'].astype('<f4').astype(float)
surfaces = {r['id']: z[r['faces']] for r in meta['surfaces']}
owners = defaultdict(list)
for pid, uses in meta['cells'].items():
    for use in uses:
        owners[use['surface']].append((pid, use))
epi = 'TA98HEART:male:epicardium'
free = []
refs = []
epifree = []
for sid, fs in surfaces.items():
    if len(owners[sid]) != 1:
        continue
    pid, use = owners[sid][0]
    face = fs if use['side'] == 1 else fs[:, [0, 2, 1]]
    free.append(face)
    refs.extend(dict(surface=sid, source_face=i) for i in range(len(fs)))
    if pid == epi:
        epifree.append(face)
assert epifree
f = np.concatenate(free)
ef = np.concatenate(epifree)
tri = p[ef]
cross = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
normal = cross / np.linalg.norm(cross, axis=1)[:, None]
sums = np.zeros_like(p)
incident = defaultdict(list)
for k in range(3):
    np.add.at(sums, ef[:, k], cross)
    for i, v in enumerate(ef[:, k]):
        incident[int(v)].append(i)
used = np.unique(ef)
directions = np.zeros_like(p)
directions[used] = -sums[used]/np.linalg.norm(sums[used], axis=1)[:, None]
fixed, solved = [], []
for v in used:
    normals = normal[incident[int(v)]]
    if np.all(normals @ directions[v] < -.01):
        continue
    result = linprog([0., 0., 0., -1.], A_ub=np.c_[normals, np.ones(len(normals))],
                     b_ub=np.zeros(len(normals)), bounds=[(-1, 1)]*3+[(0, 1)],
                     method='highs', options={'time_limit': .2})
    if result.success and np.linalg.norm(result.x[:3]) > 0:
        direction = result.x[:3]/np.linalg.norm(result.x[:3])
        if np.all(normals @ direction < -.01):
            directions[v] = direction
            solved.append(int(v))
            continue
    fixed.append(int(v))
    directions[v] = 0
# At a tissue-to-vessel junction, the calculation mask may move its shared
# edge, but the installed neighbouring tissue is never moved by this script.
depth = args.serosal_depth_m
assert 0 < depth <= .000015
q = p + directions*depth
q = q.astype('<f4').astype(float)
ids, inv = np.unique(f, return_inverse=True)
mesh = Mesh(q[ids], inv.reshape(-1, 3), 'cardiac-canonical-bank-metres')
mesh.save(out / 'mask.npz')
atomic_npz(out / 'mask-source-map.npz', source_control_ids=ids,
           source_surface_ids=np.array([r['surface'] for r in refs]),
           source_face_ids=np.array([r['source_face'] for r in refs]))
atomic_json(out / 'construction.json', dict(
    source_bank_sha256=e['bank_sha256'], source_stage=e['source_stage'],
    source_catalogue_sha256=digest_file(ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'),
    source_epi_part=epi, region_id='TA98HEART:male:pericardium:subserosa',
    serosa_term='A12.1.08.008', subserosa_term='A12.1.08.009',
    calculation_retreat_m=depth, actual_histological_thickness_verified=False,
    supported_outer_controls=len(used)-len(fixed), singular_fixed_controls=fixed,
    constrained_direction_controls=solved, no_installed_context_moved=True,
    geometry_kind='finite-tissue-mask', anatomy_accepted=False,
    elapsed_seconds=time.monotonic()-start))
print('Serosal partition mask:', len(mesh.vertices), 'controls,', len(mesh.faces),
      'faces;', len(fixed), 'fixed singular controls', flush=True)
proof = audit(mesh, require_closed=True, max_pairs=2000000, max_tests=16000000)
atomic_json(out / 'mask-audit.json', proof)
print('Calculation-mask quality:', proof['passed'], proof['failures'],
      'elapsed', time.monotonic()-start, flush=True)
assert proof['passed'], 'Calculation mask failed; no tissue partition may be installed'
