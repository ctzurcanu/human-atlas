"""Complementary finite serosa/subserosa inside the literal existing coat.

Distances to the actual exposed and myocardial interfaces set an illustrative
1:3 allocation. This is not a measured patient-specific histological ratio.
"""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.native_labelled_tetrahedra import partition_constraints_native
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

ap = argparse.ArgumentParser()
ap.add_argument('--calculation', required=True)
ap.add_argument('--stage', required=True)
ap.add_argument('--attachment-taper-m', type=float, default=0.)
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
calc = h / args.calculation
proof = json.loads((calc / 'construction.json').read_text())
assert proof['passed']
meta = json.loads((h / proof['source_stage'] / 'ownership.json').read_text())
source = np.load(h / proof['source_stage'] / 'arrangement.npz')
p = source['positions'].astype('<f4').astype(float)
surfaces = {r['id']: source[r['faces']] for r in meta['surfaces']}
epi = 'TA98HEART:male:epicardium'
uses = meta['cells'][epi]
outer = np.concatenate([surfaces[u['surface']] for u in uses if not u.get('neighbor')])
inner = np.concatenate([surfaces[u['surface']] for u in uses
                        if u.get('neighbor', '').startswith('TA98HEART:male:myocardium:')
                        or u.get('neighbor') == 'TA98HEART:male:endocardium'])
z = np.load(calc / 'tetrahedra.npz')
v, t = z['positions'], z['tetrahedra']
do = np.sqrt(igl.point_mesh_squared_distance(v, p, outer)[0])
di = np.sqrt(igl.point_mesh_squared_distance(v, p, inner)[0])
field = .25*di-.75*do
# These are literal carrier incidences, not proximity-based tissue identities.
original_ids = z['source_controls']
on_outer = np.isin(original_ids, np.unique(outer))
on_inner = np.isin(original_ids, np.unique(inner))
field[:len(original_ids)][on_outer & ~on_inner] = np.maximum(field[:len(original_ids)][on_outer & ~on_inner], 1e-10)
field[:len(original_ids)][on_inner & ~on_outer] = np.minimum(field[:len(original_ids)][on_inner & ~on_outer], -1e-10)
field[:len(original_ids)][on_outer & on_inner] = 0
if args.attachment_taper_m:
    # The old illustrative coat pinches to its neighbouring tissues at some
    # attachment stars. Terminate the new outer stratum before those singular
    # contacts; the complete original tissue stays in the inner allocation.
    assert 0 < args.attachment_taper_m <= .000002
    unsupported = do+di < args.attachment_taper_m
    field[unsupported] = -args.attachment_taper_m
else:
    unsupported = np.zeros(len(v), bool)
atomic_npz(out / 'fields.npz', positions=v, tetrahedra=t, field=field)
result = partition_constraints_native(v, t, field[:, None],
                                     [[(0, 1)], [(0, -1)]], np.array([-1]),
                                     out / 'native', seconds=60., maximum_tetrahedra=200000)
q = result['positions']
arrays = dict(positions=q)
records = []
for i, name in enumerate(['serosa', 'subserosa']):
    f = result['cells'][i]
    arrays[name+'_faces'] = f
    ids, inv = np.unique(f, return_inverse=True)
    m = Mesh(q[ids], inv.reshape(-1, 3), 'cardiac-canonical-bank-metres')
    m.save(out / (name+'.npz'))
    row = dict(name=name, vertices=len(ids), faces=len(f), geometry_sha256=m.geometry_sha)
    records.append(row)
atomic_npz(out / 'cells.npz', **arrays)
atomic_json(out / 'construction.json', dict(
    constructed=True, source_calculation_sha256=digest_file(calc / 'tetrahedra.npz'),
    source_bank_sha256=proof['source_bank_sha256'], source_stage=proof['source_stage'],
    carrier_literal=True, exterior_not_expanded=True, target_serosal_fraction=.25,
    anatomical_thickness_verified=False, cells=records,
    attachment_termination_m=args.attachment_taper_m,
    attachment_termination_control_count=int(unsupported.sum()),
    native_evidence=result['evidence'], installed=False, anatomy_accepted=False,
    elapsed_seconds=time.monotonic()-start))
print('Finite coat split saved:', [(r['name'], r['faces']) for r in records],
      'elapsed', time.monotonic()-start, flush=True)
