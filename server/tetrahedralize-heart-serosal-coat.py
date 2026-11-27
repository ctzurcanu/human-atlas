"""Calculation tetrahedra for the literal current epicardial tissue only."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
import meshpy.tet as tet

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file

ap = argparse.ArgumentParser()
ap.add_argument('--stage', required=True)
ap.add_argument('--bank', default='annular-finished-bank-v378')
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
epi = 'TA98HEART:male:epicardium'
uses = meta['cells'][epi]
f = np.concatenate([surfaces[u['surface']] if u['side'] == 1 else
                    surfaces[u['surface']][:, [0, 2, 1]] for u in uses])
ids, inv = np.unique(f, return_inverse=True)
v, fs = p[ids], inv.reshape(-1, 3)
atomic_npz(out / 'source.npz', positions=v, faces=fs, source_controls=ids)
info = tet.MeshInfo()
info.set_points((v*1000).tolist())
info.set_facets(fs.tolist())
print('Tetrahedralizing existing coat:', len(v), 'controls', len(fs), 'faces', flush=True)
result = tet.build(info, options=tet.Options('pQY'))
q = np.asarray(result.points)/1000
t = np.asarray(result.elements, np.int64)
assert np.max(np.linalg.norm(q[:len(v)]-v, axis=1)) < 1e-12
q[:len(v)] = v
det = np.linalg.det(np.transpose(q[t[:, 1:]]-q[t[:, :1]], (0, 2, 1)))
assert np.all(det != 0)
negative = det < 0
t[negative] = t[negative][:, [0, 2, 1, 3]]
volume = abs(det).sum()/6
tri = v[fs]
reference = np.einsum('ij,ij->i', tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum()/6
passed = reference > 0 and abs(volume-reference) < reference*1e-8
atomic_npz(out / 'tetrahedra.npz', positions=q, tetrahedra=t, source_faces=fs,
           source_controls=ids, original_positions=np.array(len(v)))
atomic_json(out / 'construction.json', dict(
    passed=bool(passed), source_bank_sha256=e['bank_sha256'],
    source_stage=e['source_stage'], original_positions=len(v), positions=len(q),
    tetrahedra=len(t), tissue_volume_m3=float(volume), source_volume_m3=float(reference),
    original_surface_coordinates_literal=True, original_surface_faces_literal=True,
    calculation_only=True, anatomy_accepted=False, installed=False,
    elapsed_seconds=time.monotonic()-start))
print('Existing coat calculation:', passed, len(t), 'tetrahedra;',
      'volume replay', float(volume-reference), 'm3;', time.monotonic()-start, 's', flush=True)
assert passed, 'Calculation filled an unintended space or changed tissue volume'
