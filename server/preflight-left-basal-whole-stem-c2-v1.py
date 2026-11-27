"""Read-only algebraic preflight for two graph-owned unbranched root stems."""
from pathlib import Path
import hashlib
import json
import runpy
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.tube_c2_course import C2Course

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
OUTPUT = BASE / 'whole-unbranched-root-stem-c2-v1'
SOURCE = BASE / 'root-masked-course-flow-v2/early-course-fields.npz'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
helpers = runpy.run_path(str(ROOT / 'server/construct-left-basal-c2-root-transition-v1.py'))
unit = helpers['unit']
arclength = helpers['arclength']
curvature_jet = helpers['curvature_jet']
angle_degrees = helpers['angle_degrees']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def main():
    began = time.monotonic()
    assert digest(CATALOGUE) == EXPECTED
    assert digest(SOURCE) == '89ed9551da4ffbb491fed66c6880ce8143e361fe5524797486de9f8c5dcd93e0'
    OUTPUT.mkdir(exist_ok=False)
    with np.load(SOURCE, allow_pickle=False) as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    with np.load(BASE / 'collector-left-common-basal-v1/construction-inputs.npz', allow_pickle=False) as archive:
        origin = archive['parent_center'].copy()
        plane_normal = archive['parent_distal_normal'].copy()
    slices = {int(curve): (int(start), int(end)) for curve, start, end in fields['curve_slices']}
    specs = []
    boundary_arrays = {}

    for curve, first, last, next_branch in [(131, 7, 35, 36), (60, 8, 71, 72)]:
        begin, end = slices[curve]
        original = fields['all_source_centres'][begin:end]
        held = fields['all_centres'][begin:end]
        rotations = fields['all_station_rotations'][begin:end]
        outer = fields['all_source_outer'][begin:end]
        source_tangents = unit(np.gradient(original, axis=0, edge_order=2))
        frame_tangents = np.einsum('nij,nj->ni', rotations, source_tangents)
        source_arc = arclength(original)
        held_arc = arclength(held)
        length = float(held_arc[last] - held_arc[first])
        parameters = np.linspace(held_arc[first], held_arc[last], 8)
        reference = np.column_stack([
            np.interp(parameters, held_arc, (held[:, axis] - origin[axis]) * 1000)
            for axis in range(3)
        ])
        radius = np.linalg.norm(outer[first:last + 1] - original[first:last + 1, None], axis=2).max(1)
        maximum_radius = float(radius.max()) * 1000
        minimum_bend_radius = 1.25 * maximum_radius
        diameter = max(float(np.linalg.norm(profile[:, None] - profile[None, :], axis=2).max()) for profile in outer[first:last + 1]) * 1000
        course = C2Course(reference, length, 8, minimum_bend_radius,
                          speed_lower=.85, trust_radius_mm=diameter,
                          boundary_condition='not-a-knot')
        jets = [frame_tangents[first].copy(), frame_tangents[last].copy()]
        accelerations = [
            curvature_jet(source_arc, source_tangents, first, True),
            curvature_jet(held_arc, frame_tangents, last, False),
        ]
        projections = [float(jets[0] @ course.directions[0]), float(jets[1] @ course.directions[-1])]
        acceleration_norms = [float(np.linalg.norm(vector)) for vector in accelerations]
        acceleration_bound = .85 ** 2 / minimum_bend_radius
        velocity_change = float(np.linalg.norm(jets[1] - jets[0]))
        velocity_capacity = acceleration_bound * length
        reference_compatible = min(projections) >= .85
        jet_compatible = max(acceleration_norms) <= acceleration_bound
        velocity_compatible = velocity_change <= velocity_capacity
        discrete_tangents = unit(np.gradient(held, held_arc, axis=0, edge_order=2))
        specs.append({
            'curve': curve, 'ownership_stations_inclusive': [first, last],
            'first_downstream_branch_station_fixed': next_branch,
            'physical_parameter_length_mm': length,
            'native_interval_arclength_mm': float(source_arc[last] - source_arc[first]),
            'maximum_actual_outer_circumradius_mm': maximum_radius,
            'minimum_bend_radius_mm': minimum_bend_radius,
            'speed_lower': .85, 'acceleration_bound_per_mm': acceleration_bound,
            'trust_radius_mm': diameter, 'controls': 8,
            'reference_controls_mm': reference.tolist(),
            'reference_directions': course.directions.tolist(),
            'endpoint_unit_tangents': [vector.tolist() for vector in jets],
            'endpoint_second_derivative_jets_per_mm': [vector.tolist() for vector in accelerations],
            'endpoint_jet_norms_per_mm': acceleration_norms,
            'endpoint_reference_projections': projections,
            'endpoint_frame_tangent_vs_discrete_held_derivative_degrees': [angle_degrees(jets[0], discrete_tangents[first]), angle_degrees(jets[1], discrete_tangents[last])],
            'jet_qualification': 'Exact endpoint directions are Q times original source tangent. Second jets are one-sided three-profile finite estimates of actual profile-frame tangent variation per physical arclength, normal-projected, with declared unit parameter speed.',
            'velocity_change_norm': velocity_change,
            'acceleration_integral_upper_bound': velocity_capacity,
            'velocity_change_integral_residual': velocity_capacity - velocity_change,
            'necessary_parameter_length_lower_bound_mm': velocity_change / acceleration_bound,
            'reference_cone_compatible': reference_compatible,
            'endpoint_acceleration_jets_compatible': jet_compatible,
            'velocity_integral_compatible': velocity_compatible,
            'algebraic_preflight_passed': reference_compatible and jet_compatible and velocity_compatible,
            'position_endpoint_policy': 'Auto-pinned once to reference endpoints; no duplicate position equality rows.',
        })
        for label, station in [('start', first), ('end', last), ('next_branch', next_branch)]:
            for key in ['all_outer', 'all_inner', 'all_centres', 'all_station_rotations']:
                boundary_arrays[f'curve_{curve}_{label}_{key}'] = fields[key][begin + station]

    incidence = json.loads((BASE / 'root-masked-course-flow-v1/proof.json').read_text())['declared_native_branch_incidence']
    contract = {
        'schema': 'whole-unbranched-root-stem-C2-contract-and-preflight/v1',
        'source_catalogue_sha256': EXPECTED, 'source_fields_sha256': digest(SOURCE),
        'ownership_basis': 'Whole unbranched anatomical graph stems; endpoints stop before first downstream source branch roots, not a minimally sized interval or additional collar consumption.',
        'specs': specs, 'all_seven_native_incidences_fixed': incidence,
        'construction_confluence_obligation': 'Native60root0->131station0 remains distinct from consumed-root bridge60station8->131station2. Neither adds collar consumption.',
        'fixed_fields': 'Source60prefix0..8 and131prefix0..7, exact endpoint profiles/Q, held60stations71..end and13135..end, all other six courses and their branch roots, all radii/areas and original ancestry.',
        'fixed_heart_parent_V6_skeleton_HRA_and_other_anatomy': True,
        'fixed_parent_sha256': digest(BASE / 'collector-left-common-basal-v1/fixed-current-lipv.npz'),
        'helper': 'geometry_tools/algorithms/tube_c2_course.py',
        'helper_sha256': digest(ROOT / 'geometry_tools/algorithms/tube_c2_course.py'),
        'boundary_condition': 'not-a-knot', 'solve_order_if_compatible': [131, 60],
        'solver_calls_per_interval_if_compatible': 1,
        'native_time_limit_seconds_per_call': 9, 'total_solve_budget_seconds': 20,
        'maximum_native_iterations': 90, 'no_restarts_or_parameter_sweeps': True,
        'fixed_basal_plane_origin_m': origin.tolist(),
        'fixed_basal_plane_distal_normal': plane_normal.tolist(),
        'guide_halfspace_scope': 'Guide centres only; all actual12-control plane checks, retained triangles and joins including adjacent literal strips remain mandatory before wall emission.',
        'preflight_only_no_solver_called': True,
        'both_intervals_algebraically_compatible': all(spec['algebraic_preflight_passed'] for spec in specs),
        'algebraic_scope': 'Necessary conditions of the stated C2Course continuous speed/acceleration certificate, not proof of a globally impossible anatomical route.',
        'no_new_guide_fields_or_wall': True, 'anatomical_acceptance': False,
        'elapsed_seconds': time.monotonic() - began,
    }
    np.savez_compressed(OUTPUT / 'literal-boundary-profiles.npz', **boundary_arrays)
    contract['boundary_profiles_sha256'] = digest(OUTPUT / 'literal-boundary-profiles.npz')
    write(OUTPUT / 'before-construction-contract.json', contract)
    assert digest(CATALOGUE) == EXPECTED
    print(json.dumps({'path': str(OUTPUT / 'before-construction-contract.json'),
                      'sha256': digest(OUTPUT / 'before-construction-contract.json'),
                      'both_compatible': contract['both_intervals_algebraically_compatible'],
                      'specs': specs, 'solver_calls': 0,
                      'elapsed_seconds': contract['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
