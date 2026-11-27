"""Rigid section transport onto the one saved biarc candidate, with literal joins."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
GUIDES = BASE / 'whole-unbranched-root-stem-biarc-v1'
OUTPUT = GUIDES / 'proper-profile-transport-v1'
SOURCE = BASE / 'root-masked-course-flow-v2/early-course-fields.npz'
INPUTS = BASE / 'collector-left-common-basal-v1/construction-inputs.npz'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def unit(vector):
    length = np.linalg.norm(vector, axis=-1, keepdims=True)
    assert np.all(length > 1e-14)
    return vector / length


def rotation(first, second):
    first, second = unit(first), unit(second)
    axis = np.cross(first, second)
    cosine = float(np.clip(first @ second, -1., 1.))
    assert cosine > -1. + 1e-10, 'No arbitrary antipodal frame choice'
    skew = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    return np.eye(3) + skew + skew @ skew / (1. + cosine)


def roll(axis, angle):
    axis = unit(axis)
    skew = np.array([[0., -axis[2], axis[1]], [axis[2], 0., -axis[0]], [-axis[1], axis[0], 0.]])
    return np.eye(3) + np.sin(angle) * skew + (1. - np.cos(angle)) * skew @ skew


def frame(tangent, axis):
    axis = unit(axis - float(axis @ tangent) * tangent)
    return np.column_stack([axis, unit(np.cross(tangent, axis)), tangent])


def circle(start, tangent, end):
    delta = end - start
    length = float(np.linalg.norm(delta))
    side = delta - float(delta @ tangent) * tangent
    squared = float(side @ side)
    if squared <= length * length * 1e-24:
        return {'start': start, 'tangent': tangent, 'end': end, 'length': length, 'radius': None}
    center = start + side * (length * length / (2. * squared))
    radius = float(np.linalg.norm(start - center))
    radial = (start - center) / radius
    last = (end - center) / radius
    angle = float(np.arctan2(last @ tangent, last @ radial) % (2. * np.pi))
    return {'center': center, 'radial': radial, 'tangent': tangent,
            'radius': radius, 'length': radius * angle, 'angle': angle}


def evaluate_arc(arc, distance):
    if arc['radius'] is None:
        fraction = distance / arc['length']
        return arc['start'] + fraction * (arc['end'] - arc['start']), arc['tangent']
    angle = distance / arc['radius']
    radial, tangent = arc['radial'], arc['tangent']
    direction = np.cos(angle) * radial + np.sin(angle) * tangent
    return arc['center'] + arc['radius'] * direction, -np.sin(angle) * radial + np.cos(angle) * tangent


def profile_area(profiles):
    radial = profiles - profiles.mean(axis=1)[:, None]
    return .5 * np.linalg.norm(np.cross(radial, np.roll(radial, -1, axis=1)).sum(axis=1), axis=1)


def main():
    started = time.monotonic()
    contract = json.loads((GUIDES / 'before-construction-contract.json').read_text())
    result = json.loads((GUIDES / 'guide-result.json').read_text())
    assert result['passed'] and result['contract_sha256'] == digest(GUIDES / 'before-construction-contract.json')
    assert digest(SOURCE) == contract['source_fields_sha256']
    assert digest(CATALOGUE) == contract['source_catalogue_sha256']
    independent = GUIDES / 'independent-guide-replay.json'
    assert independent.exists(), 'Independent guide replay required before transport'
    audit = json.loads(independent.read_text())
    assert audit['schema'] == 'independent-biarc-saved-guide-analytic-replay/v1'
    assert audit['passed_guide_only'], 'Independent guide replay did not pass'
    assert audit['source_guide_result_sha256'] == digest(GUIDES / 'guide-result.json')
    assert audit['source_contract_sha256'] == digest(GUIDES / 'before-construction-contract.json')
    OUTPUT.mkdir(exist_ok=False)
    before = {
        'schema': 'biarc-proper-profile-transport-contract/v1',
        'guide_contract_sha256': digest(GUIDES / 'before-construction-contract.json'),
        'guide_result_sha256': digest(GUIDES / 'guide-result.json'),
        'independent_guide_replay_sha256': digest(independent),
        'source_fields_sha256': digest(SOURCE), 'source_catalogue_sha256': digest(CATALOGUE),
        'construction_inputs_sha256': digest(INPUTS), 'script_sha256': digest(__file__),
        'ownership': {str(spec['curve']): spec['ownership_stations_inclusive'] for spec in contract['specs']},
        'parameterization': 'Original source centre polyline arclength normalized over each owned interval, mapped monotonically to exact biarc arclength; no extra rings or collar consumption.',
        'proper_transport': 'Rotation-minimizing projected source first-axis along analytic circle tangents, explicitly passing arcjoin; quintic smoothstep roll to exact distal Q. Each source section undergoes one proper rotation plus translation, no scale/shear.',
        'literal_policy': 'Both boundary positions/Q/outer/inner profiles copied literally; all noninterior sections and six other courses remain literal baseline fields.',
        'inner_provenance': '396 preexisting illustrative .8 inner sections plus exactly8 declared curve60 root extensions; none is claimed as native inner anatomy.',
        'jacobian_policy': 'Old all_jacobians renamed baseline_flow_jacobians; baseline_flow_jacobian_valid is false on every modified interior. New section normals use actual proper station Q, not old flow J.',
        'continuity': 'C1 guide with allowed one-sided curvature jump at arcjoin; exact boundary profile frames retained. Actual adjacent strip and full retained sweep gates mandatory.',
        'basal_plane_gate': 'All outer/inner controls on both owned intervals and literal curve60 profile8 must be strictly distal. Report all other retained controls separately; existing root annulus lies on the boundary and is not displaced to manufacture strict positivity.',
        'all_seven_native_incidences_fixed': contract['all_seven_native_incidences_fixed'],
        'fixed_parent_heart_V6_skeleton_all_other_anatomy': True,
        'no_wall_Boolean_install': True, 'anatomical_acceptance': False,
    }
    write(OUTPUT / 'before-transport-contract.json', before)
    with np.load(SOURCE, allow_pickle=False) as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    source_centres = fields['all_source_centres']
    source_outer, source_inner = fields['all_source_outer'], fields['all_source_inner']
    centres, outer, inner = (fields[key].copy() for key in ['all_centres', 'all_outer', 'all_inner'])
    rotations = fields['all_station_rotations'].copy()
    modified = np.zeros(len(centres), dtype=bool)
    slices = {int(curve): (int(first), int(last)) for curve, first, last in fields['curve_slices']}
    calls = {call['curve']: call for call in result['calls']}
    transport = []
    for spec in contract['specs']:
        curve = spec['curve']
        first, last = spec['ownership_stations_inclusive']
        begin, end = slices[curve]
        indices = np.arange(begin + first, begin + last + 1)
        native_centres = source_centres[begin:end]
        source_tangents = unit(np.gradient(native_centres, axis=0, edge_order=2))[first:last + 1]
        source_frames = np.array([frame(tangent, section[0] - center)
                                  for tangent, section, center in zip(source_tangents, source_outer[indices], source_centres[indices])])
        p0, p1 = np.asarray(spec['start_m']), np.asarray(spec['end_m'])
        t0, t1 = np.asarray(spec['endpoint_unit_tangents'])
        best = calls[curve]['best']
        d0 = best['d0_m']
        d1 = best['diagnostic']['second_tangent_distance_m']
        middle = ((p0 + d0 * t0) * d1 + (p1 - d1 * t1) * d0) / (d0 + d1)
        arc0, arc1 = circle(p0, t0, middle), circle(p1, -t1, middle)
        length = arc0['length'] + arc1['length']
        parameters = np.r_[0., np.cumsum(np.linalg.norm(np.diff(source_centres[indices], axis=0), axis=1))]
        fractions = parameters / parameters[-1]
        positions = fractions * length

        def evaluate(distance):
            if distance <= arc0['length']:
                return evaluate_arc(arc0, distance)
            point, tangent = evaluate_arc(arc1, length - distance)
            return point, -tangent

        target_centres = np.array([evaluate(distance)[0] for distance in positions])
        target_tangents = unit(np.array([evaluate(distance)[1] for distance in positions]))
        target_centres[[0, -1]] = [p0, p1]
        target_tangents[[0, -1]] = [t0, t1]
        base_axis = rotations[indices[0]] @ source_frames[0, :, 0]
        previous_tangent = t0
        axes = []
        previous_distance = 0.
        for distance, tangent in zip(positions, target_tangents):
            if previous_distance < arc0['length'] < distance:
                join_tangent = evaluate(arc0['length'])[1]
                base_axis = rotation(previous_tangent, join_tangent) @ base_axis
                previous_tangent = join_tangent
            base_axis = rotation(previous_tangent, tangent) @ base_axis
            base_axis = unit(base_axis - float(base_axis @ tangent) * tangent)
            axes.append(base_axis.copy())
            previous_tangent, previous_distance = tangent, distance
        axes = np.asarray(axes)
        desired_last = rotations[indices[-1]] @ source_frames[-1, :, 0]
        desired_last = unit(desired_last - float(desired_last @ t1) * t1)
        roll_angle = float(np.arctan2(t1 @ np.cross(axes[-1], desired_last), axes[-1] @ desired_last))
        blend = fractions ** 3 * (10. - 15. * fractions + 6. * fractions ** 2)
        new_rotations = np.array([frame(tangent, roll(tangent, weight * roll_angle) @ axis) @ source_frame.T
                                  for tangent, axis, weight, source_frame in zip(target_tangents, axes, blend, source_frames)])
        new_rotations[[0, -1]] = rotations[indices[[0, -1]]]
        centres[indices] = target_centres
        rotations[indices] = new_rotations
        outer[indices] = np.einsum('nij,nkj->nki', new_rotations, source_outer[indices] - source_centres[indices, None]) + target_centres[:, None]
        inner[indices] = np.einsum('nij,nkj->nki', new_rotations, source_inner[indices] - source_centres[indices, None]) + target_centres[:, None]
        for key, output in [('all_centres', centres), ('all_outer', outer), ('all_inner', inner), ('all_station_rotations', rotations)]:
            output[indices[[0, -1]]] = fields[key][indices[[0, -1]]]
        modified[indices[1:-1]] = True
        transport.append({'curve': curve, 'source_parameter_length_m': float(parameters[-1]),
                          'exact_biarc_length_m': length, 'quintic_endpoint_roll_rad': roll_angle,
                          'station_fractions': fractions.tolist(), 'analytic_tangents': target_tangents.tolist()})

    output_fields = {key: value for key, value in fields.items() if key != 'all_jacobians'}
    output_fields.update(all_centres=centres, all_outer=outer, all_inner=inner,
                         all_station_rotations=rotations,
                         baseline_flow_jacobians=fields['all_jacobians'],
                         baseline_flow_jacobian_valid=~modified,
                         biarc_modified_interior=modified)
    fields_file = OUTPUT / 'course-fields.npz'
    np.savez_compressed(fields_file, **output_fields)
    outer_ratios, inner_ratios = profile_area(outer) / profile_area(source_outer), profile_area(inner) / profile_area(source_inner)
    orthogonality = float(np.max(np.abs(np.swapaxes(rotations, 1, 2) @ rotations - np.eye(3))))
    determinant = np.linalg.det(rotations)
    fixed = ~modified
    literal = all(np.array_equal(output_fields[key][fixed], fields[key][fixed])
                  for key in ['all_centres', 'all_outer', 'all_inner', 'all_station_rotations'])
    sweep_rows, halfspaces, owned_halfspaces = [], [], []
    witnesses = {}
    with np.load(INPUTS, allow_pickle=False) as native:
        origin, plane_normal = native['parent_center'], native['parent_distal_normal']
        for curve, (begin, end) in slices.items():
            consumed = 8 if curve == 60 else 0
            faces = native[f'curve_{curve}_faces']
            original_ids = native[f'curve_{curve}_source_face_ids']
            row = {'curve': curve, 'retained_first_station': consumed, 'faces_tested': len(faces)}
            for label, profiles in [('outer', outer), ('inner', inner)]:
                vertices = profiles[begin + consumed:end].reshape(-1, 3)
                triangles = vertices[faces]
                normal = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
                radial = (triangles - centres[begin + consumed:end][faces // 12]).mean(axis=1)
                dot = np.einsum('ij,ij->i', normal, radial)
                bad = np.flatnonzero(dot <= 0.)
                join = np.zeros(len(faces), dtype=bool)
                if curve in [131, 60]:
                    spec = next(item for item in contract['specs'] if item['curve'] == curve)
                    for station in spec['ownership_stations_inclusive']:
                        global_stations = faces // 12 + consumed
                        join |= (global_stations.min(axis=1) <= station) & (global_stations.max(axis=1) >= station)
                row[label] = {'minimum_radial_dot': float(dot.min()), 'nonpositive_count': len(bad),
                              'zero_area_count': int((np.linalg.norm(normal, axis=1) == 0.).sum()),
                              'bad_native_face_ids': original_ids[bad].tolist(),
                              'bad_station_bands': sorted(set((faces[bad] // 12 + consumed).reshape(-1).tolist())),
                              'adjacent_literal_join_faces_tested': int(join.sum()),
                              'adjacent_literal_join_nonpositive_count': int((dot[join] <= 0.).sum())}
                witnesses[f'curve_{curve}_{label}_bad_triangles'] = triangles[bad]
                signed = (profiles[begin + consumed:end] - origin) @ plane_normal
                halfspaces.append({'curve': curve, 'side': label, 'minimum_signed_distance_m': float(signed.min()),
                                   'nonpositive_controls': int((signed <= 0.).sum())})
                if curve in [131, 60]:
                    spec = next(item for item in contract['specs'] if item['curve'] == curve)
                    first, last = spec['ownership_stations_inclusive']
                    owned_signed = (profiles[begin + first:begin + last + 1] - origin) @ plane_normal
                    owned_halfspaces.append({'curve': curve, 'side': label,
                                             'ownership_stations_inclusive': [first, last],
                                             'minimum_signed_distance_m': float(owned_signed.min()),
                                             'nonpositive_controls': int((owned_signed <= 0.).sum())})
            sweep_rows.append(row)
    begin60, _ = slices[60]
    profile8_plane = (outer[begin60 + 8] - origin) @ plane_normal
    incidence = []
    for edge in contract['all_seven_native_incidences_fixed']:
        child = slices[edge['child_curve']][0] + edge['child_station']
        parent = slices[edge['parent_curve']][0] + edge['parent_station']
        incidence.append({**edge, 'new_gap_m': float(np.linalg.norm(centres[child] - centres[parent])),
                          'both_endpoint_fields_literal': bool(not modified[child] and not modified[parent])})
    area_pass = bool(max(np.max(np.abs(outer_ratios - 1.)), np.max(np.abs(inner_ratios - 1.))) < 1e-10)
    rotation_pass = bool(orthogonality < 1e-10 and np.max(np.abs(determinant - 1.)) < 1e-10)
    sweep_pass = all(row[side]['nonpositive_count'] == 0 and row[side]['zero_area_count'] == 0 for row in sweep_rows for side in ['outer', 'inner'])
    # The exact fixed root collar uses curve60 profile8 as its distal mouth.
    # Report the stronger all-retained halfspace separately; do not hide crossings.
    halfspace_pass = all(row['nonpositive_controls'] == 0 for row in owned_halfspaces) and bool(np.all(profile8_plane > 0.))
    passed = literal and area_pass and rotation_pass and sweep_pass and halfspace_pass and all(edge['both_endpoint_fields_literal'] for edge in incidence)
    np.savez_compressed(OUTPUT / 'failure-witnesses.npz', **witnesses)
    proof = {'schema': 'biarc-proper-profile-actual-retained-gate/v1', 'passed': bool(passed),
             'contract_sha256': digest(OUTPUT / 'before-transport-contract.json'),
             'fields': str(fields_file), 'fields_sha256': digest(fields_file),
             'full_outer_sections': 404, 'full_inner_sections': 404,
             'inner_provenance': before['inner_provenance'], 'transport': transport,
             'modified_interior_sections': int(modified.sum()), 'all_other_profiles_Q_centres_bit_literal': literal,
             'all_outer_area_ratio_range': [float(outer_ratios.min()), float(outer_ratios.max())],
             'all_inner_area_ratio_range': [float(inner_ratios.min()), float(inner_ratios.max())],
             'proper_rotation_det_range': [float(determinant.min()), float(determinant.max())],
             'maximum_rotation_orthogonality_error': orthogonality,
             'baseline_J_not_used_for_modified_sections': True, 'baseline_J_validity_mask_false_count': int(modified.sum()),
             'retained_sweeps': sweep_rows, 'all_retained_profile_basal_halfspaces': halfspaces,
             'owned_interval_profile_basal_halfspaces': owned_halfspaces,
             'exact_curve60_profile8_halfspace_minimum_m': float(profile8_plane.min()),
             'all_seven_native_incidences': incidence,
             'rotation_and_area_passed': rotation_pass and area_pass,
             'all_actual_retained_triangles_and_adjacent_joins_passed': sweep_pass,
             'owned_profile_and_literal_root_port_halfspaces_passed': halfspace_pass,
             'no_wall_Boolean_install': True, 'anatomical_acceptance': False,
             'elapsed_seconds': time.monotonic() - started}
    assert digest(SOURCE) == before['source_fields_sha256'] and digest(CATALOGUE) == before['source_catalogue_sha256']
    write(OUTPUT / 'actual-retained-gate.json', proof)
    print(json.dumps({'proof': str(OUTPUT / 'actual-retained-gate.json'),
                      'sha256': digest(OUTPUT / 'actual-retained-gate.json'), 'passed': bool(passed),
                      'sweep_passed': sweep_pass, 'halfspace_passed': halfspace_pass,
                      'literal': literal, 'modified_sections': int(modified.sum()),
                      'bad_faces': [{ 'curve': row['curve'], 'outer': row['outer']['nonpositive_count'], 'inner': row['inner']['nonpositive_count']} for row in sweep_rows],
                      'elapsed_seconds': proof['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
