"""Finish and audit a finite, constraint-feasible course construction stage.

Optimizer convergence is recorded separately. A timed-out or incomplete
optimizer stage is never eligible; actual triangle and context review follows.
"""
from pathlib import Path
import hashlib
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit, winding_contains
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.tube_sections import section_metrics, rotation_minimizing_frames

WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1'
OUTPUT = WORK / 'completed-feasible-course-v6'
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name, value):
    (OUTPUT / name).write_text(json.dumps(value, indent=2, default=lambda x:
        x.tolist() if hasattr(x, 'tolist') else str(x)) + '\n')

started = time.monotonic()
OUTPUT.mkdir(exist_ok=False)
info = json.loads((WORK / 'source-inputs.json').read_text())
assert sha(CAT) == info['source_sha256']
template = np.load(WORK / 'superior-course-candidate-v2-clamped/source-and-candidate-fields.npz')
fields_out, checks, isometry, meshes, inputs = {}, {}, {}, {}, []
for side, suffix in [('left', 'l'), ('right', 'r')]:
    folder = WORK / 'graph-trust-v6' / side
    report = json.loads((folder / 'solve-proof.json').read_text())
    assert report['completed'] and report['elapsed_seconds'] < 32, 'Incomplete/deadline-expired input is rejected'
    assert report['source_catalogue_sha256'] == sha(CAT)
    assert report['minimum_sampled_clearance_mm'] >= .2
    assert report['maximum_middle_curvature_support'] <= .85
    z = np.load(folder / 'diagnostic-checkpoint.npz')
    original = template[side + '_source_fields']
    rings = template[side + '_rings']
    fixed = np.asarray(report['fixed_source_profiles'])
    assert np.array_equal(z['fields'][:, fixed], original[:, fixed])
    assert z['distances'].min() >= .2 and z['support'].max() <= .85
    assert np.linalg.det(z['rotations']).min() > 1 - 1e-12
    before, after = section_metrics(original), section_metrics(z['fields'])
    area_error = float(np.max(np.abs(after['area'] / before['area'] - 1)))
    shape_error = float(np.max(np.abs(after['shape_distances'] - before['shape_distances'])))
    assert area_error < 1e-11 and shape_error < 1e-12
    record = next(r for r in info['sources'] if r['part_id'] == f'ZA:Anterior cerebral artery.{suffix}')
    source = Mesh.load(WORK / record['mesh_file'])
    vertices = source.vertices.copy()
    vertices[rings] = z['fields'][0]
    attrs = {k: a.copy() for k, a in source.attributes.items()}
    n = np.einsum('sij,svj->svi', z['rotations'], source.attributes['normals'][rings].astype(float))
    n /= np.linalg.norm(n, axis=2)[..., None]
    n = np.rint(n * 32767).astype(np.int16)
    n[fixed] = source.attributes['normals'][rings[fixed]]
    attrs['normals'][rings] = n
    m = Mesh(vertices, source.faces, source.frame, attrs, dict(source_id=record['part_id'],
        source_catalogue_sha256=sha(CAT), finite_wall=False))
    check = audit(m, require_closed=False, expected_components=1, max_pairs=250000, max_tests=1000000)
    checks[side] = check
    save('actual-geometry-gates.json', dict(complete=False, checks=checks, passed=False, installed=False))
    assert check['passed'], (side, check['failures'])
    m.save(OUTPUT / f'aca-{side}.npz')
    meshes[side] = m
    isometry[side] = dict(source_id=record['part_id'], maximum_relative_section_area_error=area_error,
        maximum_section_pairwise_shape_error_m=shape_error, fixed_source_end_profiles_normals_UVs_literal=True,
        UVs_all_literal=True, source_calibre_held=True,
        optimizer_energy_converged=report['passed'],
        geometry_constraint_stage_complete=True, constraint_stage_is_not_anatomical_acceptance=True)
    for key, value in dict(rings=rings, source_fields=original, fields=z['fields'], centres=z['centres'],
        rotations=z['rotations'], source_normals=source.attributes['normals'], new_normals=attrs['normals']).items():
        fields_out[side + '_' + key] = value
    inputs.append(dict(side=side, report=str(folder / 'solve-proof.json'), report_sha256=sha(folder / 'solve-proof.json'),
        checkpoint_sha256=sha(folder / 'diagnostic-checkpoint.npz'), optimizer_energy_converged=False,
        geometry_constraints_met=True, no_deadline_expiration=True))

radius = .0006871863787458552
ends = np.array([fields_out['left_centres'][13], fields_out['right_centres'][13]])
path = np.linspace(*ends, 9)
frames = rotation_minimizing_frames(path)
angle = np.arange(16) * np.pi / 8
field = path[:, None] + radius * (np.cos(angle)[None, :, None] * frames[:, None, :, 0]
    + np.sin(angle)[None, :, None] * frames[:, None, :, 1])
v = field.reshape(-1, 3)
faces = []
for k in range(8):
    for j in range(16):
        a, b = k * 16 + j, k * 16 + (j + 1) % 16
        faces.extend([[a, b, a + 16], [b, b + 16, a + 16]])
faces = np.asarray(faces)
radial = (field - path[:, None]).reshape(-1, 3) / radius
alignment = np.einsum('ij,ij->i', np.cross(v[faces[:, 1]] - v[faces[:, 0]], v[faces[:, 2]] - v[faces[:, 0]]), radial[faces].mean(1))
if alignment.sum() < 0:
    faces = faces[:, ::-1].copy()
uv = np.stack(np.meshgrid(np.arange(16) / 16, np.linspace(0, np.linalg.norm(ends[1] - ends[0]), 9), indexing='xy'), -1).reshape(-1, 2).astype(np.float32)
bridge = Mesh(v, faces, meshes['left'].frame,
    dict(normals=np.rint(radial * 32767).astype(np.int16), uvs=uv),
    dict(independently_constructed=True, source_primitive_not_copied=True, finite_wall=False))
checks['bridge'] = audit(bridge, require_closed=False, expected_components=1, max_pairs=250000, max_tests=1000000)
assert checks['bridge']['passed']
bridge.save(OUTPUT / 'acom-open-bridge.npz')
fields_out.update(bridge_fields=field, bridge_centres=path, bridge_radius_m=np.asarray(radius))
np.savez_compressed(OUTPUT / 'source-and-candidate-fields.npz', **fields_out)
save('section-isometry-proof.json', isometry)
save('construction-input-certification.json', dict(complete=True, inputs=inputs, passed=True,
    scope='A completed physical-constraint certificate and actual source-surface replay; not convergence to an energy minimum or anatomical acceptance.',
    elapsed_seconds=time.monotonic() - started, installed=False))
save('candidate-manifest.json', dict(schema='ACA-superior-communication-course-candidate/v1',
    source_catalogue_sha256=sha(CAT), frame=bridge.frame, bridge_radius_m=radius,
    intended_carotid_terminal_union_required=True, finite_wall=False, installed=False,
    meshes=[dict(part_id=f'ZA:Anterior cerebral artery.{suffix}', candidate=f'aca-{side}.npz')
        for side, suffix in [('left', 'l'), ('right', 'r')]] + [dict(part_id='ZA:Anterior communicating artery', candidate='acom-open-bridge.npz')]))
save('actual-geometry-gates.json', dict(complete=True, checks=checks, passed=True,
    source_catalogue_unchanged=True, installed=False, anatomical_acceptance=False))
print(json.dumps(dict(output=str(OUTPUT), passed=True, left_self=checks['left']['self_intersections']['defect_pairs'],
    right_self=checks['right']['self_intersections']['defect_pairs'], elapsed_seconds=time.monotonic() - started)))
