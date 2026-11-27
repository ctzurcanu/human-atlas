"""Staged reuse of the frozen basal close_branch/outer-minus-blood constructor.

Preparation only binds a new full404 profile field. Construction additionally
requires an explicit passing, SHA-bound retained-sweep/port gate. No live atlas,
anatomical assignment, packing or installation is performed by this driver.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit, winding_contains
from geometry_tools.algorithms.branch_junction_envelopes import close_branch
from geometry_tools.algorithms.finite_tissue_boolean import make_solid, provenance, barycentric
from geometry_tools.algorithms.tube_continuity import boundary_loops

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
FROZEN = BASE / 'collector-left-common-basal-v1'
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
CAT_SHA = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
INPUT_SHA = 'fd97fcea6a46624aa621226c172a63c7a425d22cac29f8eea58f4c582087bc98'
FRAME = 'model:local-ta98:assembled-metres'
CURVES = [60, 61, 62, 63, 131, 132, 133, 134]
GRAPH = {60: (131, 2), 61: (60, 84), 62: (60, 72), 63: (62, 24),
         132: (131, 36), 133: (132, 24), 134: (131, 60)}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=lambda x: x.tolist() if hasattr(x, 'tolist') else x.item()) + '\n')


def archive(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k].copy() for k in z.files}


def bindings(contract):
    assert sha(CAT) == contract['source_catalogue_sha256'] == CAT_SHA
    for key in ['fields', 'frozen_inputs', 'fixed_parent', 'field_proof']:
        row = contract[key]
        assert sha(row['path']) == row['sha256'], 'Changed input: ' + key


def prepare(fields_path, field_proof, output):
    start = time.monotonic()
    assert sha(CAT) == CAT_SHA and sha(FROZEN / 'construction-inputs.npz') == INPUT_SHA
    assert not output.exists()
    fields = archive(fields_path)
    frozen = archive(FROZEN / 'construction-inputs.npz')
    slices = {int(c): (int(a), int(b)) for c, a, b in fields['curve_slices']}
    assert set(slices) == set(CURVES)
    assert all(fields[k].shape == (404, 12, 3) for k in ['all_outer', 'all_inner', 'all_source_outer', 'all_source_inner'])
    assert fields['all_centres'].shape == (404, 3)
    assert sum(b-a for a, b in slices.values()) == 404
    arrays = {k: v for k, v in frozen.items() if k.startswith('parent_')}
    rows = []
    for curve in CURVES:
        a, b = slices[curve]
        prefix = 'curve_' + str(curve)
        trim = 8 if curve == 60 else 0
        old_outer = frozen[prefix + '_outer'].reshape(-1, 12, 3)
        old_inner = frozen[prefix + '_inner'].reshape(-1, 12, 3)
        expected_native_outer = frozen[prefix + '_retained_native_sections'] if curve == 131 else old_outer
        assert np.array_equal(fields['all_source_outer'][a+trim:b], expected_native_outer)
        assert np.array_equal(fields['all_source_inner'][a+trim:b], old_inner)
        outer = fields['all_outer'][a+trim:b].copy()
        inner = fields['all_inner'][a+trim:b].copy()
        correspondence = None
        if curve == 131:
            # The current encoded parent is canonical; the literal native root
            # polygon has a disclosed 0.173275mm phase mismatch. This is the
            # same declared strip0 bridge, not a new course or radius change.
            _, correspondence = cKDTree(outer[0]).query(arrays['parent_outer'])
            assert len(set(correspondence)) == 12
            outer[0, correspondence] = arrays['parent_outer']
            inner[0, correspondence] = arrays['parent_inner']
        arrays[prefix + '_outer'] = outer.reshape(-1, 3)
        arrays[prefix + '_inner'] = inner.reshape(-1, 3)
        for suffix in ['_faces', '_source_face_ids', '_retained_native_sections']:
            arrays[prefix + suffix] = frozen[prefix + suffix]
        arrays[prefix + '_native_station_ids'] = np.arange(trim, b-a)
        arrays[prefix + '_posed_centres'] = fields['all_centres'][a+trim:b]
        if curve in GRAPH:
            parent, station = GRAPH[curve]
            first, _ = slices[parent]
            arrays[prefix + '_internal_anchor'] = fields['all_centres'][first+station]
        if 'all_station_rotations' in fields:
            arrays[prefix + '_station_rotations'] = fields['all_station_rotations'][a+trim:b]
        # Original UV and normal charts remain source-bound. Proper transport
        # rotations are exported separately; no original UV is regenerated.
        source_file = FROZEN / ('source-inferior-basal.npz' if curve < 100 else 'source-superior-basal.npz')
        source = Mesh.load(source_file)
        error, ids = cKDTree(source.vertices).query(fields['all_source_outer'][a+trim:b].reshape(-1, 3))
        # Canonical curve131 profile0 replaced the native polygon, so its
        # source attributes are separately held by the exact parent packet.
        literal = error == 0
        arrays[prefix + '_literal_source_control_mask'] = literal
        arrays[prefix + '_nearest_source_control_ids'] = ids
        for name, values in source.attributes.items():
            arrays[prefix + '_original_source_' + name] = values[ids]
        rows.append(dict(curve=curve, source_profile_count=b-a, retained_profile_count=len(outer), trim_native_stations=list(range(trim)),
                         constructed_parent=GRAPH.get(curve), native_root_identity='60:0 -> 131:0' if curve == 60 else None,
                         root_loop=1 if curve == 63 else 0, root_polygon_correspondence=None if correspondence is None else correspondence.tolist(),
                         source_npz=str(source_file), source_npz_sha256=sha(source_file)))
    output.mkdir(parents=True)
    np.savez_compressed(output / 'construction-inputs.npz', **arrays)
    sweep_checks = []
    for curve in CURVES:
        prefix = 'curve_' + str(curve)
        faces = arrays[prefix + '_faces']
        for side in ['outer', 'inner']:
            points = arrays[prefix + '_' + side]
            centres = points.reshape(-1, 12, 3).mean(1)
            radial = points - np.repeat(centres, 12, axis=0)
            t = points[faces]
            dots = np.einsum('ij,ij->i', np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), radial[faces].mean(1))
            check = audit(Mesh(points, faces, FRAME), require_closed=False, expected_components=1)
            sweep_checks.append(dict(curve=curve, side=side, nonpositive_radial_faces=np.flatnonzero(dots <= 0).tolist(),
                                     minimum_radial_dot=float(dots.min()), audit=check))
    write(output / 'prepared-retained-sweep-gate.json', dict(passed=all(not p['nonpositive_radial_faces'] and p['audit']['passed'] for p in sweep_checks), checks=sweep_checks))
    assert all(not p['nonpositive_radial_faces'] and p['audit']['passed'] for p in sweep_checks), 'Stop before Boolean: actual prepared sweep failed.'
    contract = dict(schema='left-basal-new-field-construction-contract/v1', source_catalogue_sha256=CAT_SHA,
                    fields=dict(path=str(fields_path), sha256=sha(fields_path)),
                    field_proof=dict(path=str(field_proof), sha256=sha(field_proof)),
                    frozen_inputs=dict(path=str(FROZEN / 'construction-inputs.npz'), sha256=INPUT_SHA),
                    fixed_parent=dict(path=str(FROZEN / 'fixed-current-lipv.npz'), sha256=sha(FROZEN / 'fixed-current-lipv.npz')),
                    construction_inputs_sha256=sha(output / 'construction-inputs.npz'), rows=rows, frame=FRAME,
                    prepared_retained_gate_sha256=sha(output / 'prepared-retained-sweep-gate.json'),
                    incidence=GRAPH, inner_provenance='396 existing illustrative inner profiles + 8 explicitly extended consumed curve60 root profiles; all404 transported, never regenerated from posed outer.',
                    source_preservation='Full404 transported first; consume ONLY60stations0..7. Construct60retained8 -> 131station2; native incidence0 -> 0 remains distinct. Curve131strip0 is exact parent-rim bridge; all other prepared fields literal.',
                    fixed_interfaces='All5960 current LIPV faces/coordinates/normalwords/UV/signs, atrial seams and V6 annulus literal. Only new collector geometry is emitted.',
                    algorithm='Existing close_branch junction sleeves + one outer union minus one connected fluid union. Fresh per-primitive face provenance; no inherited output face IDs.',
                    target='Whole communicating8-course basal wall before any .405 root-only ownership partition.',
                    required_gate_schema=dict(passed=True, field_sha256=sha(fields_path), retained_sweeps_passed=True, profile8_halfspace_passed=True),
                    numerical_ancestry_only=True, anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-start)
    write(output / 'contract-before-construction.json', contract)
    print(json.dumps(dict(prepared=str(output), contract_sha256=sha(output / 'contract-before-construction.json'), elapsed_seconds=contract['elapsed_seconds'])))


def construct(output, gate_path):
    start = time.monotonic()
    contract_path = output / 'contract-before-construction.json'
    contract = json.loads(contract_path.read_text())
    bindings(contract)
    assert sha(output / 'construction-inputs.npz') == contract['construction_inputs_sha256']
    assert sha(output / 'prepared-retained-sweep-gate.json') == contract['prepared_retained_gate_sha256']
    assert json.loads((output / 'prepared-retained-sweep-gate.json').read_text())['passed']
    gate = json.loads(gate_path.read_text())
    assert all(gate.get(k) is True for k in ['passed', 'retained_sweeps_passed', 'profile8_halfspace_passed'])
    assert gate['field_sha256'] == contract['fields']['sha256']
    candidate = output / 'candidate-stage-v1'
    assert not candidate.exists(), 'One construction only; preserve failed stages.'
    candidate.mkdir()
    write(candidate / 'state.json', dict(complete=False, passed=False, installed=False, gate_sha256=sha(gate_path)))
    x = archive(output / 'construction-inputs.npz')
    primitives = []
    unions = {'outer': [], 'fluid': []}
    component_rows = []
    for curve in CURVES:
        prefix = 'curve_' + str(curve)
        anchors = {1 if curve == 63 else 0: x[prefix + '_internal_anchor']} if curve in GRAPH else {}
        outer, fluid = close_branch(x[prefix + '_outer'], x[prefix + '_inner'], x[prefix + '_faces'], anchors,
                                    origin=x['parent_center'], fluid_through_external=True)
        component_rows.append(dict(curve=curve, outer_mouths=outer['mouths'], fluid_mouths=fluid['mouths']))
        for role, data in [('outer', outer), ('fluid', fluid)]:
            solid, oid = make_solid(data['positions'], data['faces'])
            unions[role].append(solid)
            primitives.append(dict(curve=curve, role=role, oid=oid, **data))
    # Every declared internal anchor must actually be inside the posed parent
    # fluid volume. Classification is computational, not anatomical acceptance.
    anchor_checks = []
    for curve, (parent, station) in GRAPH.items():
        primitive = next(p for p in primitives if p['curve'] == parent and p['role'] == 'fluid')
        check = winding_contains(Mesh(primitive['positions'], primitive['faces'], FRAME), x['curve_' + str(curve) + '_internal_anchor'][None])
        anchor_checks.append(dict(curve=curve, parent=parent, parent_station=station, certificate=check))
        assert check['inside'] == [True] and check['boundary_ambiguous'] == [False]
    import manifold3d as md
    outer_union = md.Manifold.batch_boolean(unions['outer'], md.OpType.Add)
    blood_union = md.Manifold.batch_boolean(unions['fluid'], md.OpType.Add)
    wall = outer_union - blood_union
    lookup = {p['oid']: p for p in primitives}
    primitive_arrays = {}
    for i, p in enumerate(primitives):
        for name in ['positions', 'faces', 'disk_faces', 'junction_faces']:
            primitive_arrays['primitive_' + str(i) + '_' + name] = p[name]
    np.savez_compressed(candidate / 'source-primitives.npz', **primitive_arrays)
    exported = {}
    for name, solid in [('outer-classifier', outer_union), ('blood-classifier', blood_union), ('finite-basal-wall', wall)]:
        assert solid.status() == md.Error.NoError
        data = solid.to_mesh64()
        positions = np.asarray(data.vert_properties)[:, :3]
        faces = np.asarray(data.tri_verts, int)
        origins, source_faces = provenance(data)
        mesh = Mesh(positions, faces, FRAME)
        mesh.save(candidate / (name + '.npz'))
        parent_triangles = np.array([lookup[int(o)]['positions'][lookup[int(o)]['faces'][int(f)]] for o, f in zip(origins, source_faces)])
        weights = np.array([barycentric(t, p) for t, p in zip(mesh.triangles, parent_triangles)])
        np.savez_compressed(candidate / (name + '-ancestry.npz'), primitive_ids=origins, source_face_ids=source_faces,
                            parent_triangles=parent_triangles, corner_barycentric=weights)
        exported[name] = mesh
    wall_check = audit(exported['finite-basal-wall'], require_closed=True, expected_components=1, expected_genus=[8])
    blood_check = audit(exported['blood-classifier'], require_closed=True, expected_components=1, expected_genus=[0])
    proof = dict(schema='left-basal-new-field-finite-construction/v1', source_catalogue_sha256=CAT_SHA,
                 contract_sha256=sha(contract_path), inputs_sha256=sha(output / 'construction-inputs.npz'), field_gate_sha256=sha(gate_path),
                 field_sha256=contract['fields']['sha256'], primitives=[{k: v for k, v in p.items() if k not in ['positions', 'faces', 'disk_faces', 'junction_faces']} for p in primitives],
                 component_construction=component_rows, anchor_checks=anchor_checks,
                 wall_audit=wall_check, blood_audit=blood_check, wall_npz_sha256=sha(candidate / 'finite-basal-wall.npz'),
                 blood_npz_sha256=sha(candidate / 'blood-classifier.npz'), source_primitives_sha256=sha(candidate / 'source-primitives.npz'),
                 canonical_annulus_pending=True, actual_nine_port_replay_pending=True, fixed_context_pending=True,
                 anatomical_acceptance=False, installed=False, complete=True,
                 passed=wall_check['passed'] and blood_check['passed'], elapsed_seconds=time.monotonic()-start)
    write(candidate / 'proof.json', proof)
    write(candidate / 'state.json', proof)
    print(json.dumps(dict(output=str(candidate), passed=proof['passed'], elapsed_seconds=proof['elapsed_seconds'])))
    assert proof['passed'], 'Stop on positive geometry gate; no automatic cleanup or retry.'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare', 'construct'])
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--fields', type=Path)
    parser.add_argument('--field-proof', type=Path)
    parser.add_argument('--gate', type=Path)
    args = parser.parse_args()
    if args.stage == 'prepare':
        assert args.fields and args.field_proof
        prepare(args.fields.resolve(), args.field_proof.resolve(), args.output.resolve())
    else:
        assert args.gate
        construct(args.output.resolve(), args.gate.resolve())


if __name__ == '__main__':
    main()
