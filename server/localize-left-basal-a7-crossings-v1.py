"""Read-only exact ancestry and source-pose localization of the saved A7 contacts."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import exact_intersection_segment
from geometry_tools.common import atomic_json, json_bytes

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
CONSTRUCTION = BASE / 'biarc-field-finite-wall-v1'
V1, V2 = CONSTRUCTION / 'candidate-stage-v1', CONSTRUCTION / 'candidate-stage-v2'
OUTPUT = V2 / 'a7-crossing-localization-v1'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FIELDS = BASE / 'whole-unbranched-root-stem-biarc-v1/proper-profile-transport-v1/course-fields.npz'
HELD = BASE / 'root-masked-course-flow-v2/early-course-fields.npz'
NATIVE = ROOT / '.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1/native-curves-v3'
ARTERY = 'ZA:Medial basal segmental artery of left lung'
FRAME = 'model:local-ta98:assembled-metres'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def intervals(values):
    values = sorted(set(map(int, values)))
    groups = []
    for value in values:
        if not groups or value > groups[-1][1] + 1:
            groups.append([value, value])
        else:
            groups[-1][1] = value
    return groups


def bounds(points):
    return [points.min(axis=0).tolist(), points.max(axis=0).tolist()] if len(points) else None


def nearest_segments(points, curves):
    starts, ends, labels = [], [], []
    for curve, centers in curves.items():
        starts.extend(centers[:-1])
        ends.extend(centers[1:])
        labels.extend((curve, station) for station in range(len(centers) - 1))
    starts, ends = np.asarray(starts), np.asarray(ends)
    delta = ends - starts
    parameter = np.clip(np.einsum('nkj,kj->nk', points[:, None] - starts, delta) / np.sum(delta * delta, axis=1), 0., 1.)
    projected = starts + parameter[:, :, None] * delta
    squared = np.sum((points[:, None] - projected) ** 2, axis=2)
    chosen = squared.argmin(axis=1)
    return np.asarray([labels[index] for index in chosen]), parameter[np.arange(len(points)), chosen], np.sqrt(squared[np.arange(len(points)), chosen])


def main():
    started = time.monotonic()
    assert digest(CATALOGUE) == '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
    assert digest(V2 / 'finite-basal-wall.npz') == 'e5e65bc043e92de2a179a5ef65f8db4ea6b5ef5c7b0267cb019862cb2f8e810e'
    OUTPUT.mkdir(exist_ok=False)
    context_proof = json.loads((V2 / 'independent-fixed-context-v1/proof.json').read_text())
    context_row = next(row for row in context_proof['completed_owner_results'] if row['id'] == ARTERY)
    contact_file = Path(context_row['contact_archive'])
    assert digest(contact_file) == context_row['contact_archive_sha256']
    contact = load(contact_file)
    pairs = contact['pairs'][contact['relation'] == 'crossing']
    assert len(pairs) == 502
    wall = Mesh.load(V2 / 'finite-basal-wall.npz')
    atlas = json.loads(CATALOGUE.read_text())
    parts = {part['id']: part for part in atlas['parts']}
    reader = AtlasReader(ROOT)
    artery = reader.read_part(atlas, parts[ARTERY], frame=FRAME, catalogue_path=CATALOGUE, catalogue_sha=digest(CATALOGUE))
    assert artery.geometry_sha == context_row['context_geometry_sha256']
    assert artery.source['bank_sha256'] == context_row['bank_sha256']
    ancestry = load(V2 / 'patch-and-ancestry.npz')
    primitives = load(V1 / 'source-primitives.npz')
    specification = json.loads((V1 / 'proof.json').read_text())['primitives']
    by_oid = {row['oid']: (index, row) for index, row in enumerate(specification)}
    fields, held = load(FIELDS), load(HELD)
    original_inputs = load(BASE / 'collector-left-common-basal-v1/construction-inputs.npz')
    slices = {int(curve): (int(first), int(last)) for curve, first, last in fields['curve_slices']}
    affected = np.unique(pairs[:, 0])
    # All contacts must be outside the canonical root patch; do not guess new ancestry.
    assert affected.max() < len(ancestry['retained_corner_barycentric'])
    rows, source_patch, held_patch = [], [], []
    maximum_replay_error = 0.
    for face in affected:
        oid, source_face = int(ancestry['source_primitive_ids'][face]), int(ancestry['source_face_ids'][face])
        primitive_index, spec = by_oid[oid]
        curve, role = spec['curve'], spec['role']
        begin, end = slices[curve]
        consumed = 8 if curve == 60 else 0
        native_faces = original_inputs[f'curve_{curve}_faces']
        assert source_face < len(native_faces), 'A cap/anchor requires a distinct source comparison'
        primitive_faces = primitives[f'primitive_{primitive_index}_faces']
        face_controls = primitive_faces[source_face]
        assert np.array_equal(face_controls, native_faces[source_face])
        primitive_positions = primitives[f'primitive_{primitive_index}_positions']
        barycentric = ancestry['retained_corner_barycentric'][face]
        replay = barycentric @ primitive_positions[face_controls]
        maximum_replay_error = max(maximum_replay_error, float(np.max(np.abs(replay - wall.triangles[face]))))
        source_key = 'all_source_outer' if role == 'outer' else 'all_source_inner'
        posed_key = 'all_outer' if role == 'outer' else 'all_inner'
        native_vertices = fields[source_key][begin + consumed:end].reshape(-1, 3)
        held_vertices = held[posed_key][begin + consumed:end].reshape(-1, 3)
        source_patch.append(barycentric @ native_vertices[face_controls])
        held_patch.append(barycentric @ held_vertices[face_controls])
        rows.append({'wall_face': int(face), 'primitive_oid': oid, 'curve': curve, 'role': role,
                     'primitive_face': source_face,
                     'native_source_face': int(original_inputs[f'curve_{curve}_source_face_ids'][source_face]),
                     'native_stations': sorted(set((face_controls // 12 + consumed).tolist()))})
    assert maximum_replay_error < 1e-11, 'Saved literal ancestry did not replay the final contact triangles'
    source_patch, held_patch = np.asarray(source_patch), np.asarray(held_patch)
    states = {'source_original_ancestor_patch': source_patch,
              'held_physical_flow_v2_ancestor_patch': held_patch,
              'current_biarc_final_patch': wall.triangles[affected]}
    state_rows, contact_arrays = {}, {}
    for label, triangles in states.items():
        mesh = Mesh(triangles.reshape(-1, 3), np.arange(triangles.size // 3).reshape(-1, 3), FRAME)
        result, summary = pair_contacts(mesh, artery, 0., max_tests=200000, max_pairs=100000)
        assert not summary['truncated']
        crossing = result['pairs'][result['relation'] == 'crossing']
        state_rows[label] = {'geometry_sha256': mesh.geometry_sha, 'summary': summary,
                             'crossing_pairs': len(crossing),
                             'final_wall_face_ids_with_crossings': affected[np.unique(crossing[:, 0])].tolist()}
        contact_arrays[label + '_crossing_pairs'] = np.column_stack([affected[crossing[:, 0]], crossing[:, 1]])
    native_vein = reader.read_part(atlas, parts['ZA:Inferior basal vein of left lung'], frame=FRAME, catalogue_path=CATALOGUE, catalogue_sha=digest(CATALOGUE))
    native_contacts, native_summary = pair_contacts(native_vein, artery, 0., max_tests=200000, max_pairs=100000)
    assert not native_summary['truncated']
    contact_arrays['full_current_native_inferior_basal_source_crossing_pairs'] = native_contacts['pairs'][native_contacts['relation'] == 'crossing']

    exact_segments, exact_errors, exact_relation = [], [], []
    for first, second in pairs:
        result = exact_intersection_segment(wall.triangles[first], artery.triangles[second])
        exact_relation.append(result['relation'])
        exact_segments.append(result['endpoints'])
        exact_errors.append(result['coordinate_error_upper_m'])
    assert all(relation == 'crossing' for relation in exact_relation)
    exact_segments, exact_errors = np.asarray(exact_segments), np.asarray(exact_errors)
    evidence = json.loads((NATIVE / 'evidence.json').read_text())
    arterial_specs = [row for row in evidence['curves'] if row['part_id'] == ARTERY]
    sections = load(NATIVE / 'literal-section-centres.npz')
    arterial_curves = {row['index']: sections[row['key'] + '_centers'] for row in arterial_specs}
    labels, fractions, distances = nearest_segments(exact_segments.mean(axis=1), arterial_curves)
    arterial_face_rows = []
    coordinate_lookup = {}
    for spec in arterial_specs:
        profiles = sections[spec['key'] + '_literal_sections']
        for station, profile in enumerate(profiles):
            for point in profile:
                coordinate_lookup.setdefault(tuple(point), []).append((spec['index'], station))
    for face in np.unique(pairs[:, 1]):
        matches = [coordinate_lookup.get(tuple(point), []) for point in artery.triangles[face]]
        assert all(matches), 'Arterial contact triangle not bound to literal native profiles'
        curve_ids = set.intersection(*(set(curve for curve, _ in values) for values in matches))
        assert len(curve_ids) == 1
        curve = curve_ids.pop()
        stations = sorted(set(station for values in matches for which, station in values if which == curve))
        arterial_face_rows.append({'source_face': int(face), 'curve': curve, 'native_stations': stations})
    grouped = []
    for curve in sorted(set(row['curve'] for row in rows)):
        chosen = [row for row in rows if row['curve'] == curve]
        selected = np.array([row['wall_face'] for row in chosen])
        pairmask = np.isin(pairs[:, 0], selected)
        begin, _ = slices[curve]
        stations = sorted(set(station for row in chosen for station in row['native_stations']))
        current_motion = np.linalg.norm(fields['all_centres'][begin + np.array(stations)] - fields['all_source_centres'][begin + np.array(stations)], axis=1)
        biarc_motion = np.linalg.norm(fields['all_centres'][begin + np.array(stations)] - held['all_centres'][begin + np.array(stations)], axis=1)
        grouped.append({'curve': curve, 'source_part_id': 'ZA:Inferior basal vein of left lung',
                        'crossing_pairs': int(pairmask.sum()), 'final_wall_faces': len(chosen),
                        'roles': sorted(set(row['role'] for row in chosen)),
                        'native_station_intervals': intervals(stations),
                        'native_source_face_intervals': intervals(row['native_source_face'] for row in chosen),
                        'intersection_segment_bounds_m': bounds(exact_segments[pairmask].reshape(-1, 3)),
                        'source_to_current_center_displacement_minmax_m': [float(current_motion.min()), float(current_motion.max())],
                        'held_flow_to_biarc_center_displacement_minmax_m': [float(biarc_motion.min()), float(biarc_motion.max())]})
    nearest = [{'arterial_curve': curve, 'native_segment_intervals': intervals(labels[labels[:, 0] == curve, 1]),
                'crossing_witnesses_nearest': int((labels[:, 0] == curve).sum()),
                'centerline_distance_minmax_m': [float(distances[labels[:, 0] == curve].min()), float(distances[labels[:, 0] == curve].max())]}
               for curve in sorted(set(labels[:, 0]))]
    arterial_groups = [{'curve': curve, 'component': next(row['component'] for row in arterial_specs if row['index'] == curve),
                        'actual_source_contact_faces': sum(row['curve'] == curve for row in arterial_face_rows),
                        'native_station_intervals': intervals(station for row in arterial_face_rows if row['curve'] == curve for station in row['native_stations'])}
                       for curve in sorted(set(row['curve'] for row in arterial_face_rows))]
    arrays_file = OUTPUT / 'literal-localization-witnesses.npz'
    np.savez_compressed(arrays_file, exact_intersection_segments=exact_segments, exact_endpoint_coordinate_error_upper_m=exact_errors,
                        original_ancestor_triangles=source_patch, held_flow_ancestor_triangles=held_patch,
                        final_wall_face_ids=affected, nearest_arterial_curve_segment=labels,
                        nearest_arterial_segment_fraction=fractions, nearest_arterial_centerline_distance_m=distances,
                        **contact_arrays)
    report = {'schema': 'source-bound-final-basal-wall-A7-crossing-localization/v1',
              'source_catalogue_sha256': digest(CATALOGUE), 'candidate_wall_sha256': digest(V2 / 'finite-basal-wall.npz'),
              'contact_archive_sha256': digest(contact_file), 'final_ancestry_sha256': digest(V2 / 'patch-and-ancestry.npz'),
              'source_primitives_sha256': digest(V1 / 'source-primitives.npz'), 'field_sha256': digest(FIELDS),
              'held_field_sha256': digest(HELD), 'artery_geometry_sha256': artery.geometry_sha,
              'artery_bank_sha256': artery.source['bank_sha256'], 'script_sha256': digest(__file__),
              'crossing_pairs': len(pairs), 'distinct_wall_faces': len(affected), 'distinct_arterial_faces': len(arterial_face_rows),
              'exact_rational_crossings_replayed': len(exact_segments),
              'exact_intersection_segment_bounds_m': bounds(exact_segments.reshape(-1, 3)),
              'maximum_endpoint_Float64_conversion_error_m': float(exact_errors.max()),
              'affected_venous_courses': grouped, 'affected_arterial_components': arterial_groups,
              'nearest_arterial_centerline_segments': nearest,
              'all_final_contact_face_ancestry': rows, 'all_arterial_contact_face_stations': arterial_face_rows,
              'source_pose_comparison_same_descendant_face_patch': state_rows,
              'full_native_inferior_basal_source_vs_current_artery': native_summary,
              'ancestor_triangle_replay_max_coordinate_error_m': maximum_replay_error,
              'smallest_actual_dependency_set': 'Inferior-basal vein curve60 contact band plus its curve62 branch and their literal branch incidence; actual A7 artery components listed above. No superior-basal131 tree, atrial root, V6 or whole-lung refit inferred.',
              'interpretation': 'Same face-ancestry patch is replayed in native, held-flow and current poses. Pair counts compare surfaces only, not occupied volume or correct anatomy. A new crossing does not establish which vessel should move; source-relative sections and declared branch incidence constrain any later proposal.',
              'source_inner_qualification': 'Native outer is literal; source inner is previously declared illustrative .8 radial fraction, not native wall anatomy.',
              'witnesses_sha256': digest(arrays_file), 'no_fit_reroute_or_source_geometry_change': True,
              'anatomical_acceptance': False, 'elapsed_seconds': time.monotonic() - started}
    output = OUTPUT / 'proof.json'
    atomic_json(output, report)
    assert digest(CATALOGUE) == report['source_catalogue_sha256']
    print(json_bytes({'proof': str(output), 'sha256': digest(output), 'affected_vein_courses': grouped,
                      'affected_arterial_components': arterial_groups, 'nearest_arterial_segments': nearest,
                      'state_crossings': {label: row['crossing_pairs'] for label, row in state_rows.items()},
                      'native_whole_source_crossings': native_summary['crossing_face_pairs'],
                      'overlap_bounds_m': report['exact_intersection_segment_bounds_m'],
                      'elapsed_seconds': report['elapsed_seconds']}).decode(), flush=True)


if __name__ == '__main__':
    main()
