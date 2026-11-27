"""One bounded biarc guide per graph-owned basal stem; no profile or wall edits."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.tube_biarc import biarc_path
from geometry_tools.algorithms.tube_sections import TubeGeometryError

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
OUTPUT = BASE / 'whole-unbranched-root-stem-biarc-v1'
SOURCE = BASE / 'root-masked-course-flow-v2/early-course-fields.npz'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
OLD_CONTRACT = BASE / 'whole-unbranched-root-stem-c2-v1/before-construction-contract.json'
PARENT = BASE / 'collector-left-common-basal-v1/fixed-current-lipv.npz'
HELPER = ROOT / 'geometry_tools/algorithms/tube_biarc.py'
EXPECTED_SOURCE = '89ed9551da4ffbb491fed66c6880ce8143e361fe5524797486de9f8c5dcd93e0'
EXPECTED_CATALOGUE = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def verify_sources(contract):
    for path, key in [(CATALOGUE, 'source_catalogue_sha256'),
                      (SOURCE, 'source_fields_sha256'),
                      (PARENT, 'fixed_parent_sha256'),
                      (HELPER, 'helper_sha256'),
                      (OLD_CONTRACT, 'previous_contract_sha256')]:
        assert digest(path) == contract[key], f'Changed binding: {path}'


def prepare():
    assert digest(CATALOGUE) == EXPECTED_CATALOGUE
    assert digest(SOURCE) == EXPECTED_SOURCE
    OUTPUT.mkdir(exist_ok=False)
    old = json.loads(OLD_CONTRACT.read_text())
    with np.load(SOURCE, allow_pickle=False) as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    slices = {int(curve): (int(first), int(last))
              for curve, first, last in fields['curve_slices']}
    specs = []
    boundaries = {}
    for previous in old['specs']:
        curve = previous['curve']
        first, last = previous['ownership_stations_inclusive']
        begin, _ = slices[curve]
        start, end = fields['all_centres'][[begin + first, begin + last]]
        tangent0, tangent1 = np.asarray(previous['endpoint_unit_tangents'])
        chord = end - start
        alpha, beta = float(chord @ tangent0), float(chord @ tangent1)
        coefficient = 1. - float(tangent0 @ tangent1)
        numerator_constant = float(chord @ chord) / 2.
        # Both actual cases have N=A-alpha*d0 and D=beta+c*d0.
        # alpha,beta>0 and c>=0 imply d1>0 exactly on 0<d0<A/alpha.
        assert alpha > 0. and beta > 0. and coefficient >= 0.
        upper = numerator_constant / alpha
        native = fields['all_source_outer'][begin + first:begin + last + 1]
        centers = fields['all_source_centres'][begin + first:begin + last + 1]
        radii = np.linalg.norm(native - centers[:, None], axis=2).max(axis=1)
        radius = float(radii.max())
        assert abs(radius * 1250. - previous['minimum_bend_radius_mm']) < 1e-10
        spec = {
            'curve': curve, 'ownership_stations_inclusive': [first, last],
            'first_downstream_branch_station_fixed': previous['first_downstream_branch_station_fixed'],
            'start_m': start.tolist(), 'end_m': end.tolist(),
            'endpoint_unit_tangents': [tangent0.tolist(), tangent1.tolist()],
            'endpoint_frame_vs_discrete_held_tangent_degrees': previous['endpoint_frame_tangent_vs_discrete_held_derivative_degrees'],
            'endpoint_jet_qualification': previous['jet_qualification'],
            'maximum_actual_outer_circumradius_m': radius,
            'maximum_radius_source_station': first + int(radii.argmax()),
            'minimum_bend_radius_m': 1.25 * radius,
            'd1_formula': '(A-alpha*d0)/(beta+c*d0)',
            'A_m2': numerator_constant, 'alpha_m': alpha,
            'beta_m': beta, 'c': coefficient,
            'mathematically_legal_d0_interval_m_open': [0., upper],
            'legal_interval_derivation': 'N>0 iff d0<A/alpha; D>0 for all d0>0 because beta>0 and c>=0. Hence d1>0 exactly on this open interval.',
        }
        specs.append(spec)
        for label, station in [('start', first), ('end', last)]:
            for key in ['all_outer', 'all_inner', 'all_centres', 'all_station_rotations']:
                boundaries[f'curve_{curve}_{label}_{key}'] = fields[key][begin + station]
    boundary_file = OUTPUT / 'literal-boundary-profiles.npz'
    np.savez_compressed(boundary_file, **boundaries)
    contract = {
        'schema': 'whole-unbranched-root-stem-biarc-contract/v1',
        'source_catalogue_sha256': digest(CATALOGUE),
        'source_fields_sha256': digest(SOURCE),
        'fixed_parent_sha256': digest(PARENT),
        'previous_contract_sha256': digest(OLD_CONTRACT),
        'helper': str(HELPER.relative_to(ROOT)), 'helper_sha256': digest(HELPER),
        'script_sha256': digest(__file__),
        'literal_boundary_profiles_sha256': digest(boundary_file),
        'specs': specs,
        'continuity': 'C1 tangent-continuous circular arcs; curvature may jump at their shared tangent. C2 is not claimed.',
        'method_reason': 'The failed held-reference cone schedule is algebraically incompatible; this method has exact whole-arc radius bounds without frozen reference directions.',
        'fixed_fields': old['fixed_fields'],
        'all_seven_native_incidences_fixed': old['all_seven_native_incidences_fixed'],
        'construction_confluence_obligation': old['construction_confluence_obligation'],
        'fixed_heart_parent_V6_skeleton_HRA_and_other_anatomy': True,
        'fixed_basal_plane_origin_m': old['fixed_basal_plane_origin_m'],
        'fixed_basal_plane_distal_normal': old['fixed_basal_plane_distal_normal'],
        'scalar_method': 'scipy.optimize.minimize_scalar(method=bounded); one call per stem,131 then60, stop at first failure; no parameter sweep or restart.',
        'objective': 'minimize max(1/R_first,1/R_second), straight arc curvature=0',
        'termination': {'maximum_objective_evaluations_per_stem': 64,
                        'maximum_seconds_per_stem': 10., 'maximum_total_seconds': 20.,
                        'absolute_d0_tolerance_m': 1e-10},
        'coarse_sampling_for_exact_radii': {'max_step_m': 1e6, 'max_turn_rad': float(np.pi / 2), 'chord_error_m': 1e6},
        'fine_path_only_after_radius_pass': {'max_step_m': .0005, 'max_turn_rad': float(np.deg2rad(2)), 'chord_error_m': 1e-6},
        'profile_transport_and_sweep_gates': 'Separate stage mandatory: exact endpoint Q/profiles, unscaled sections, adjacent literal strips, all retained faces and basal halfspace; no wall before pass.',
        'no_profile_field_or_wall_changes': True, 'anatomical_acceptance': False,
    }
    write(OUTPUT / 'before-construction-contract.json', contract)
    verify_sources(contract)
    print(json.dumps({'contract': str(OUTPUT / 'before-construction-contract.json'),
                      'sha256': digest(OUTPUT / 'before-construction-contract.json'),
                      'scalar_intervals_m': [spec['mathematically_legal_d0_interval_m_open'] for spec in specs]}))


def solve():
    contract_file = OUTPUT / 'before-construction-contract.json'
    contract = json.loads(contract_file.read_text())
    verify_sources(contract)
    assert digest(__file__) == contract['script_sha256']
    result_file = OUTPUT / 'guide-result.json'
    assert not result_file.exists(), 'One stage only; preserve previous evidence'
    started = time.monotonic()
    result = {'schema': 'whole-unbranched-root-stem-biarc-guide-result/v1',
              'contract_sha256': digest(contract_file), 'calls': [],
              'passed': False, 'guide_only': True, 'wall_emitted': False}
    for spec in contract['specs']:
        curve = spec['curve']
        call_started = time.monotonic()
        history = []
        best = None

        def objective(distance):
            nonlocal best
            if len(history) >= 64 or time.monotonic() - call_started >= 10. or time.monotonic() - started >= 20.:
                raise TimeoutError('Declared finite scalar budget reached')
            try:
                _, _, _, diagnostic = biarc_path(
                    spec['start_m'], spec['endpoint_unit_tangents'][0],
                    spec['end_m'], spec['endpoint_unit_tangents'][1],
                    first_tangent_distance_m=distance,
                    **contract['coarse_sampling_for_exact_radii'])
                radii = [arc['radius_m'] for arc in diagnostic['arcs']]
                curvature = max(0. if radius is None else 1. / radius for radius in radii)
                entry = {'d0_m': float(distance), 'maximum_curvature_per_m': curvature,
                         'arc_radii_m': radii, 'diagnostic': diagnostic}
                if best is None or curvature < best['maximum_curvature_per_m']:
                    best = entry
            except TubeGeometryError as error:
                curvature = float('inf')
                entry = {'d0_m': float(distance), 'geometry_failure': str(error)}
            history.append(entry)
            return curvature

        try:
            optimized = minimize_scalar(
                objective, method='bounded',
                bounds=tuple(spec['mathematically_legal_d0_interval_m_open']),
                options={'xatol': 1e-10, 'maxiter': 64})
            termination = {'success': bool(optimized.success), 'status': int(optimized.status),
                           'message': str(optimized.message), 'nfev': int(optimized.nfev),
                           'x_m': float(optimized.x), 'fun_per_m': float(optimized.fun)}
        except TimeoutError as error:
            termination = {'success': False, 'message': str(error)}
        call = {'curve': curve, 'termination': termination, 'evaluations': history,
                'best': best, 'elapsed_seconds': time.monotonic() - call_started}
        result['calls'].append(call)
        if not termination['success'] or best is None:
            result.update(status='ScalarTerminationFailed', stopped_curve=curve)
            break
        radii = best['arc_radii_m']
        call['whole_arc_radius_bound_pass'] = all(radius is None or radius >= spec['minimum_bend_radius_m'] for radius in radii)
        if not call['whole_arc_radius_bound_pass']:
            result.update(status='PhysicalRadiusGateFailed', stopped_curve=curve)
            break
        points, tangents, curvature, fine = biarc_path(
            spec['start_m'], spec['endpoint_unit_tangents'][0],
            spec['end_m'], spec['endpoint_unit_tangents'][1],
            first_tangent_distance_m=best['d0_m'],
            minimum_bend_radius_m=spec['minimum_bend_radius_m'],
            **contract['fine_path_only_after_radius_pass'])
        assert np.array_equal(points[[0, -1]], np.asarray([spec['start_m'], spec['end_m']]))
        path = OUTPUT / f'curve-{curve}-guide.npz'
        np.savez_compressed(path, points=points, tangents=tangents,
                            curvature=curvature, endpoint_tangents=np.asarray(spec['endpoint_unit_tangents']))
        call.update(fine_path=str(path), fine_path_sha256=digest(path),
                    fine_diagnostic=fine, stations=len(points),
                    endpoint_tangent_max_error=float(np.max(np.abs(tangents[[0, -1]] - spec['endpoint_unit_tangents']))))
    else:
        result.update(status='BothWholeArcRadiusBoundsPassed', passed=True)
    result['elapsed_seconds'] = time.monotonic() - started
    result['actual_profile_sweep_and_join_gates'] = 'Pending separate stage; no profile field, wall or CSG emitted.'
    verify_sources(contract)
    write(result_file, result)
    review = ROOT / 'reports/mesh-audit/REVIEW.md'
    with review.open('a') as stream:
        stream.write('\n\n### Left basal whole-stem biarc diagnostic,2026-10-25\n\n')
        stream.write(f"Source346bed; contract `{digest(contract_file)}`. One bounded scalar call per attempted stem; status `{result['status']}`. C1 guide only, exact endpoint positions/profile-frame tangent bindings, whole-arc radius criterion1.25×actual outer circumradius; no calibre reduction. Actual proper-profile transport, retained/join strips and port halfspace remain pending. No wall or installation and no anatomical acceptance. Evidence: `{result_file.relative_to(ROOT)}`.\n")
    print(json.dumps({'result': str(result_file), 'sha256': digest(result_file),
                      'status': result['status'], 'calls': [{key: call.get(key) for key in ['curve', 'whole_arc_radius_bound_pass', 'elapsed_seconds', 'stations']} for call in result['calls']],
                      'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['prepare', 'solve'], required=True)
    arguments = parser.parse_args()
    if arguments.stage == 'prepare':
        prepare()
    else:
        solve()
