"""One bounded two-interval C2 guide feasibility attempt. No tissue wall emission."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.tube_c2_course import C2Course

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
SOURCE = BASE / 'root-masked-course-flow-v2'
OUTPUT = BASE / 'curvature-constrained-root-transition-v1'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def unit(vectors):
    lengths = np.linalg.norm(vectors, axis=-1, keepdims=True)
    assert np.all(lengths > 1e-12)
    return vectors / lengths


def arclength(points):
    return np.r_[0., np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))] * 1000


def curvature_jet(parameters, tangents, station, incoming):
    ids = np.arange(station - 2, station + 1) if incoming else np.arange(station, station + 3)
    local = parameters[ids] - parameters[station]
    derivative = np.polyfit(local, tangents[ids], 2)[1]
    tangent = tangents[station]
    return derivative - np.dot(derivative, tangent) * tangent


def angle_degrees(first, second):
    return float(np.degrees(np.arccos(np.clip(np.dot(first, second), -1, 1))))


def main():
    began = time.monotonic()
    assert digest(CATALOGUE) == EXPECTED
    OUTPUT.mkdir(exist_ok=False)
    source_fields = SOURCE / 'early-course-fields.npz'
    assert digest(source_fields) == '89ed9551da4ffbb491fed66c6880ce8143e361fe5524797486de9f8c5dcd93e0'
    with np.load(source_fields, allow_pickle=False) as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    with np.load(BASE / 'collector-left-common-basal-v1/construction-inputs.npz', allow_pickle=False) as archive:
        plane_origin = archive['parent_center'].copy()
        plane_normal = archive['parent_distal_normal'].copy()
    fixed_parent = BASE / 'collector-left-common-basal-v1/fixed-current-lipv.npz'
    fixed_parent_sha = digest(fixed_parent)
    origin = plane_origin.copy()
    slices = {int(curve): (int(start), int(end)) for curve, start, end in fields['curve_slices']}
    specs = []
    courses = {}
    boundary_arrays = {}

    for curve, first, last in [(131, 7, 20), (60, 8, 32)]:
        begin, end = slices[curve]
        original = fields['all_source_centres'][begin:end]
        held = fields['all_centres'][begin:end]
        rotations = fields['all_station_rotations'][begin:end]
        outer = fields['all_source_outer'][begin:end]
        source_tangents = unit(np.gradient(original, axis=0, edge_order=2))
        frame_tangents = np.einsum('nij,nj->ni', rotations, source_tangents)
        source_parameters = arclength(original)
        held_parameters = arclength(held)
        length = float(held_parameters[last] - held_parameters[first])
        knot_parameters = np.linspace(held_parameters[first], held_parameters[last], 8)
        reference = np.column_stack([
            np.interp(knot_parameters, held_parameters, (held[:, axis] - origin[axis]) * 1000)
            for axis in range(3)
        ])
        radius = np.linalg.norm(outer[first:last + 1] - original[first:last + 1, None], axis=2).max(1)
        minimum_radius = 1.25 * float(radius.max()) * 1000
        diameter = max(float(np.linalg.norm(profile[:, None] - profile[None, :], axis=2).max()) for profile in outer[first:last + 1]) * 1000
        course = C2Course(reference, length, 8, minimum_radius,
                          speed_lower=.85, trust_radius_mm=diameter,
                          boundary_condition='not-a-knot')
        start_tangent = frame_tangents[first].copy()
        end_tangent = frame_tangents[last].copy()
        start_acceleration = curvature_jet(source_parameters, source_tangents, first, True)
        end_acceleration = curvature_jet(held_parameters, frame_tangents, last, False)
        equalities = [
            (0., 1, start_tangent, f'curve{curve}:literal-start-profile-frame'),
            (0., 2, start_acceleration, f'curve{curve}:incoming-profile-frame-curvature-estimate'),
            (length, 1, end_tangent, f'curve{curve}:literal-held-end-profile-frame'),
            (length, 2, end_acceleration, f'curve{curve}:outgoing-profile-frame-curvature-estimate'),
        ]
        discrete = unit(np.gradient(held, held_parameters, axis=0, edge_order=2))
        start_projection = float(start_tangent @ course.directions[0])
        end_projection = float(end_tangent @ course.directions[-1])
        halfspaces = [(float(s), plane_normal.tolist(), 0., {'kind': 'fixed-basal-plane-centre', 'curve': curve})
                      for s in np.linspace(0, length, 64)]
        spec = {
            'curve': curve, 'ownership_stations_inclusive': [first, last],
            'physical_parameter_length_mm': length, 'controls': 8,
            'minimum_radius_mm': minimum_radius, 'speed_lower': .85,
            'trust_radius_mm': diameter, 'reference_controls_mm': reference.tolist(),
            'reference_directions': course.directions.tolist(),
            'endpoint_unit_jet_projections': [start_projection, end_projection],
            'endpoint_tangents_from_literal_Q_times_original_source_tangent': [start_tangent.tolist(), end_tangent.tolist()],
            'tangent_vs_discrete_held_derivative_degrees': [angle_degrees(start_tangent, discrete[first]), angle_degrees(end_tangent, discrete[last])],
            'endpoint_second_derivative_jets_per_mm': [start_acceleration.tolist(), end_acceleration.tolist()],
            'second_jet_qualification': 'One-sided three-profile finite estimates of actual profile-frame tangent variation per physical arclength; normal projection; unit endpoint speed. These are declared jets, not certified native continuous derivatives.',
            'position_equalities': 'Pinned once by C2Course reference endpoints; no duplicate position equality rows.',
            'profile_boundary_policy': 'Complete endpoint profile positions and proper Q frames remain literal after the guide solve. Numerical endpoint residual never authorizes changing shared controls.',
            'roll_policy_if_guide_passes': 'Proper station frames along analytic guide tangent, rotation-minimizing transport with a C2 roll correction to the exact endpoint Q. Actual join strips and all retained12-control sweeps must pass.',
            'plane_policy': 'Conic guide-centre halfspaces only; actual full-profile distal-plane and triangle gates are still mandatory before wall emission.',
        }
        specs.append(spec)
        courses[curve] = (course, halfspaces, equalities)
        for label, station in [('start', first), ('end', last)]:
            index = begin + station
            for key in ['all_outer', 'all_inner', 'all_centres', 'all_station_rotations']:
                boundary_arrays[f'curve_{curve}_{label}_{key}'] = fields[key][index]

    incidence = json.loads((BASE / 'root-masked-course-flow-v1/proof.json').read_text())['declared_native_branch_incidence']
    contract = {
        'schema': 'one-calibrated-two-interval-C2-course-contract/v1',
        'source_catalogue_sha256': EXPECTED, 'input_fields_sha256': digest(source_fields),
        'source_curvature_diagnostic_sha256': digest(SOURCE / 'curvature-feasibility-v1/diagnostic.json'),
        'helper': 'geometry_tools/algorithms/tube_c2_course.py',
        'helper_sha256': digest(ROOT / 'geometry_tools/algorithms/tube_c2_course.py'),
        'solve_order': [131, 60], 'one_call_per_interval': True,
        'native_time_limit_seconds_per_call': 9, 'total_solve_budget_seconds': 20,
        'maximum_native_iterations': 90, 'no_restarts_or_parameter_sweeps': True,
        'specs': specs, 'all_seven_native_branch_incidence_obligations': incidence,
        'construction_confluence_distinguished': 'Native60root0->131station0 is protected. Consumed-root construction60retained8->131station2 remains a separate fixed confluence obligation; no extra consumption.',
        'fixed_fields': 'Source60prefix0..8, source131prefix0..7; held60stations32..end and13120..end, all other six courses, literal boundary profiles/Q, all radii/areas/ancestry.',
        'fixed_parent_path': str(fixed_parent), 'fixed_parent_sha256': fixed_parent_sha,
        'fixed_parent_heart_V6_skeleton_and_HRA_context': True,
        'all_retained_ring_triangle_and_adjacent_literal_strip_gates_required': True,
        'no_wall_or_Boolean_emission': True, 'anatomical_acceptance': False,
    }
    write(OUTPUT / 'before-construction-contract.json', contract)
    np.savez_compressed(OUTPUT / 'literal-boundary-profiles.npz', **boundary_arrays)
    result = {'schema': 'one-calibrated-two-interval-C2-course-result/v1',
              'contract_sha256': digest(OUTPUT / 'before-construction-contract.json'),
              'calls': [], 'passed': False, 'wall_emitted': False,
              'profile_sweep_gate_not_waived': True}

    for spec in specs:
        curve = spec['curve']
        if min(spec['endpoint_unit_jet_projections']) < .85:
            result.update(status='ReferenceConeAlgebraicallyIncompatible', stopped_curve=curve,
                          witness=spec['endpoint_unit_jet_projections'],
                          scope='Exact endpoint unit derivative cannot meet the fixed reference projection lower bound; no solve was called.')
            break
        remaining = 20. - sum(call['elapsed_seconds'] for call in result['calls'])
        if remaining <= 0:
            result.update(status='TotalSolveBudgetExhausted', stopped_curve=curve)
            break
        course, inequalities, equalities = courses[curve]
        controls, diagnostic = course.solve(inequalities, time_limit=min(9., remaining), equalities=equalities)
        result['calls'].append({'curve': curve, **diagnostic})
        if controls is None:
            result.update(status=diagnostic['status'], stopped_curve=curve,
                          scope='One declared conic guide problem only; not global anatomical path infeasibility. No restart or bound relaxation.')
            break
        np.savez_compressed(OUTPUT / f'curve-{curve}-guide.npz', controls_mm=controls,
                            knots_mm=course.knots, reference_mm=course.reference)
        if not diagnostic['certificate']['continuous_curvature_bound_pass']:
            result.update(status='ContinuousGuideCertificateFailed', stopped_curve=curve)
            break
    else:
        result.update(status='BothGuideCertificatesPassed', passed=True,
                      scope='Guide certificates only; actual proper-frame profile transport, full sweeps and join strips still pending. No wall may be emitted.')

    assert digest(CATALOGUE) == EXPECTED
    assert digest(fixed_parent) == fixed_parent_sha
    result['elapsed_seconds'] = time.monotonic() - began
    write(OUTPUT / 'result.json', result)
    print(json.dumps({'output': str(OUTPUT), 'result_sha256': digest(OUTPUT / 'result.json'),
                      'status': result['status'], 'stopped_curve': result.get('stopped_curve'),
                      'solve_calls': len(result['calls']), 'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
