"""Current source-bound arterial cohort/interface inventory; no fitting."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from geometry_tools.spatial import BVH
from geometry_tools.winding_backend import winding_numbers

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
CURRENT = BASE / 'male-intrapulmonary-circulation-current-v1'
OUTPUT = CURRENT / 'left-lower-arterial-dependency-inventory-v1'
FRONTIER = BASE / 'pulmonary-circulation-frontier-v1'
NATIVE = FRONTIER / 'native-curves-v3'
INVENTORY = BASE / 'pulmonary-family-completion-v1/arterial-inventory-v1/inventory-with-geometry.json'
UPPER = BASE / 'pulmonary-family-completion-v1/left-upper-arterial-contract-v1/whole-left-upper-artery-contract.json'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
FRAME = 'model:local-ta98:assembled-metres'
PARTS = {
    'A6_source': 'ZA:Superior segmental artery of left lung',
    'A7_source': 'ZA:Medial basal segmental artery of left lung',
    'A8_source': 'ZA:Anterior basal segmental artery of left lung',
    'A9_source': 'ZA:Lateral basal segmental artery of left lung',
    'A10_source': 'ZA:Posterior basal segmental artery of left lung',
}
TERRITORIES = {'A6_source': 'superior', 'A7_source': 'medial-basal',
               'A8_source': 'anterior-basal', 'A9_source': 'lateral-basal', 'A10_source': 'posterior-basal'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path):
    with np.load(path, allow_pickle=False) as archive:
        return {key: archive[key].copy() for key in archive.files}


def ring_relation(first, second):
    distance = np.linalg.norm(first[:, None] - second[None, :], axis=2)
    phase = min(float(np.linalg.norm(first - np.roll(other, shift, axis=0), axis=1).max())
                for other in [second, second[::-1]] for shift in range(len(second)))
    return {'centre_gap_m': float(np.linalg.norm(first.mean(axis=0) - second.mean(axis=0))),
            'bidirectional_control_Hausdorff_m': float(max(distance.min(axis=0).max(), distance.min(axis=1).max())),
            'best_cyclic_or_reversed_max_control_gap_m': phase,
            'literal_shared_control_count': len(set(map(tuple, first)) & set(map(tuple, second))),
            'literal_shared_ring': phase == 0., 'communicating_lumen_established': False}


def main():
    began = time.monotonic()
    assert digest(CATALOGUE) == EXPECTED
    OUTPUT.mkdir(exist_ok=False)
    atlas = json.loads(CATALOGUE.read_text())
    parts = {part['id']: part for part in atlas['parts']}
    hidden = set(atlas.get('defaultHidden', []))
    inventory = json.loads(INVENTORY.read_text())
    assert inventory['source_catalogue_sha256'] == EXPECTED and inventory['inclusive_terms'] == 44
    upper = json.loads(UPPER.read_text())
    assert upper['source_catalogue_sha256'] == EXPECTED
    old_geometry = {row['id']: row for row in inventory['source_geometry']}
    native = json.loads((NATIVE / 'evidence.json').read_text())
    arrays_file = NATIVE / 'literal-section-centres.npz'
    assert digest(arrays_file) == native['arrays_sha256']
    arrays = load(arrays_file)
    curves = [curve for curve in native['curves'] if curve['part_id'] in PARTS.values()]
    assert len(curves) == 21
    by_curve = {curve['index']: curve for curve in curves}
    reader, meshes, source_rows = AtlasReader(ROOT), {}, []
    source_ids = list(PARTS.values()) + ['ZA:Left pulmonary artery', 'ZA:Pulmonary trunk', 'ZA:Bifurcation of pulmonary trunk', 'ZA:Right pulmonary artery']
    for identity in source_ids:
        part = parts[identity]
        mesh = reader.read_part(atlas, part, frame=FRAME, catalogue_path=CATALOGUE, catalogue_sha=EXPECTED)
        if identity in old_geometry:
            assert mesh.geometry_sha == old_geometry[identity]['geometry_sha256']
            assert mesh.source['bank_sha256'] == old_geometry[identity]['bank_sha256']
        meshes[identity] = mesh
        source_rows.append({'id': identity, 'source_id': part.get('sourceId'),
                            'concept_id': part.get('conceptId'), 'name': part['name'],
                            'geometry_sha256': mesh.geometry_sha, 'bank_sha256': mesh.source['bank_sha256'],
                            'record_sha256': hashlib.sha256(json.dumps(part, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                            'bank_path': mesh.source['bank_path'],
                            'actual_bounds_m': [mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)],
                            'suppressed': bool(part.get('suppressed')), 'catalogue_default_hidden': identity in hidden,
                            'record_default_hidden': bool(part.get('defaultHidden')),
                            'finite_wall_claimed': bool(part.get('pairedTissue') and part.get('sectionSolid')),
                            'native_courses': [curve['index'] for curve in curves if curve['part_id'] == identity]})
    for curve in curves:
        assert curve['geometry_sha256'] == meshes[curve['part_id']].geometry_sha

    # Reuse only source geometry ancestry after current literal hash replay.
    graph_file = FRONTIER / 'rooted-source-recipient-graph-v17/evidence.json'
    graph = json.loads(graph_file.read_text())['vascular_source_graph']
    incidences = []
    neighbors = {curve: set() for curve in by_curve}
    for node in graph['nodes']:
        members = [member for member in node['curve_station_members'] if member['curve'] in by_curve]
        if len(members) < 2:
            continue
        for first, second in combinations(members, 2):
            one, two = by_curve[first['curve']], by_curve[second['curve']]
            first_ring = arrays[one['key'] + '_literal_sections'][first['station']]
            second_ring = arrays[two['key'] + '_literal_sections'][second['station']]
            neighbors[first['curve']].add(second['curve'])
            neighbors[second['curve']].add(first['curve'])
            incidences.append({'source_graph_node': node['key'], 'first': first, 'second': second,
                               'cross_source_object': first['part_id'] != second['part_id'],
                               **ring_relation(first_ring, second_ring),
                               'geometry_only_ancestry_not_arterial_target': True})
    components, unseen = [], set(by_curve)
    while unseen:
        stack, group = [min(unseen)], []
        while stack:
            curve = stack.pop()
            if curve not in unseen:
                continue
            unseen.remove(curve)
            group.append(curve)
            stack.extend(neighbors[curve])
        components.append({'curves': sorted(group),
                           'source_part_ids': sorted(set(by_curve[curve]['part_id'] for curve in group)),
                           'source_graph_center_incidence_only': True})

    lpa = meshes['ZA:Left pulmonary artery']
    lpa_points = set(map(tuple, lpa.vertices))
    lpa_tree = BVH(lpa)
    lpa_ends = []
    for end in upper['fixed_heart_and_proximal_recipients']['full_current_LPA_ends']:
        ring = np.asarray(end['complete_bound_literal_loop']['positions_m'])
        assert all(tuple(point) in lpa_points for point in ring)
        pair_rows = []
        for curve in curves:
            profiles = arrays[curve['key'] + '_literal_sections']
            for endpoint in [0, len(profiles) - 1]:
                pair_rows.append({'curve': curve['index'], 'native_station': endpoint,
                                  **ring_relation(ring, profiles[endpoint])})
        pair_rows.sort(key=lambda row: row['centre_gap_m'])
        lpa_ends.append({'current_boundary_index': end['current_boundary_index'],
                         'center_m': ring.mean(axis=0), 'literal_current_ring_fixed': True,
                         'protected_proximal_role_not_anatomically_established': True,
                         'nearest_native_lower_arterial_endpoint_rings': pair_rows[:6],
                         'closest_fixed_heart_recipient_witnesses': end['closest_fixed_heart_recipient_witnesses'],
                         'actual_shared_arterial_join_established': False})
    endpoint_lpa_witnesses = []
    for curve in curves:
        for station in [0, curve['stations'] - 1]:
            center = arrays[curve['key'] + '_centers'][station]
            closest = lpa_tree.closest_point(center)
            endpoint_lpa_witnesses.append({'curve': curve['index'], 'station': station,
                                           'center_m': center, 'LPA_nearest_surface_face': closest['face_id'],
                                           'LPA_surface_distance_m': closest['distance_m'],
                                           'not_an_attachment_or_lumen_claim': True})

    # Actual current recipient geometry; no source-label fit is performed.
    targets = {}
    union_vertices, union_faces, offset, target_bindings = [], [], 0, []
    target_ids = ['TA98PULMONARY:male:left-inferior:' + territory for territory in TERRITORIES.values()]
    for identity in target_ids:
        part = parts[identity]
        target = reader.read_part(atlas, part, frame=FRAME, catalogue_path=CATALOGUE, catalogue_sha=EXPECTED)
        targets[identity] = target
        keep = np.concatenate([np.arange(group['start'] // 3, (group['start'] + group['count']) // 3)
                               for group in part['tissueDrawGroups'] if group.get('neighbor') not in target_ids])
        union_vertices.append(target.vertices)
        union_faces.append(target.faces[keep] + offset)
        offset += len(target.vertices)
        target_bindings.append({'id': identity, 'geometry_sha256': target.geometry_sha,
                                'bank_sha256': target.source['bank_sha256']})
    vertices, inverse = np.unique(np.concatenate(union_vertices), axis=0, return_inverse=True)
    union = Mesh(vertices, inverse[np.concatenate(union_faces)], FRAME)
    assert topology(union)['closed_oriented_manifold']
    placement, winding_arrays = [], {}
    for label, identity in PARTS.items():
        source_curves = [curve for curve in curves if curve['part_id'] == identity]
        centers = np.concatenate([arrays[curve['key'] + '_centers'] for curve in source_curves])
        outer = np.concatenate([arrays[curve['key'] + '_literal_sections'].reshape(-1, 3) for curve in source_curves])
        target = targets['TA98PULMONARY:male:left-inferior:' + TERRITORIES[label]]
        whole_winding = np.abs(winding_numbers(union.vertices, union.faces, centers))
        named_winding = np.abs(winding_numbers(target.vertices, target.faces, centers))
        control_winding = np.abs(winding_numbers(target.vertices, target.faces, outer))
        placement.append({'label_only_not_accepted_identity': label, 'source_id': identity,
                          'courses': [curve['index'] for curve in source_curves],
                          'whole_current_lower_lobe_centres_inside': int((whole_winding > .5).sum()),
                          'centre_count': len(centers),
                          'nominal_named_territory_id': 'TA98PULMONARY:male:left-inferior:' + TERRITORIES[label],
                          'nominal_named_territory_centres_inside': int((named_winding > .5).sum()),
                          'nominal_named_territory_outer_controls_inside': int((control_winding > .5).sum()),
                          'outer_control_count': len(outer),
                          'geometric_mismatch_not_individual_anatomical_fault': True})
        winding_arrays[label + '_whole_abs_winding'] = whole_winding
        winding_arrays[label + '_named_abs_winding'] = named_winding
        winding_arrays[label + '_outer_named_abs_winding'] = control_winding
    winding_path = OUTPUT / 'source-recipient-winding-samples.npz'
    np.savez_compressed(winding_path, **winding_arrays)
    term_rows = []
    for term in inventory['terms']:
        code = term['authority']['id']
        if code != 'A12.2.01.201' and code not in {f'A12.2.01.{number}' for number in range(214, 221)}:
            continue
        fields = term['source_fields']
        all_ids = fields['resolved_direct_part_ids']
        visible = [identity for identity in all_ids if identity in parts and not parts[identity].get('suppressed') and not parts[identity].get('defaultHidden') and identity not in hidden]
        term_rows.append({'authority': term['authority'], 'recorded_source_reference_status': term['current_representation_status'],
                          'source_binding_ids_in_existing_report': all_ids,
                          'actually_active_direct_binding_ids': visible,
                          'retained_hidden_or_suppressed_binding_ids': [identity for identity in all_ids if identity not in visible],
                          'nominal_native_named_sources_are_not_binding_proof': True})
    lingular = upper['all_nine_lingular_closure']
    for path, binding in lingular['all_banks'].items():
        assert digest(path) == binding['sha256']
    bronchial_file = CURRENT / 'current-left-lower-bronchial-landmarks-v1/evidence.json'
    assert digest(bronchial_file) == '9e75680801f2deafe12228ce8fa179cac61a8f0dcfaf6e623b455ba2df7d17fd'
    bronchial = json.loads(bronchial_file.read_text())
    assert bronchial['source_catalogue_sha256'] == EXPECTED
    report = {'schema': 'current-left-lower-pulmonary-arterial-dependency-interface-inventory/v1',
              'source_catalogue_sha256': EXPECTED, 'whole44_inventory': str(INVENTORY),
              'whole44_inventory_sha256': digest(INVENTORY), 'whole44_inclusive_terms': 44,
              'relevant_current_authority_terms': term_rows, 'source_native_graph': str(graph_file),
              'source_native_graph_sha256': digest(graph_file),
              'old_graph_catalogue_binding_is_not_relabelled': True,
              'old_graph_reuse_qualification': 'Only literal source node membership reused after all21 course geometry hashes and current source banks replay; old target transforms and anatomical assignments not reused.',
              'native_literal_sections_sha256': digest(arrays_file), 'current_sources': source_rows,
              'no_standalone_native_basal_part_mesh': 'ZA:Basal part of left pulmonary artery' not in parts,
              'source_courses': curves, 'source_center_node_incidence_components': components,
              'all_source_station_ring_incidence_metrics': incidences,
              'current_LPA_open_end_interfaces': lpa_ends,
              'all_source_endpoints_to_current_LPA_surface_witnesses': endpoint_lpa_witnesses,
              'protected_main_PT_bifurcation_rightPA_heart': [row for row in source_rows if row['id'] in ['ZA:Pulmonary trunk', 'ZA:Bifurcation of pulmonary trunk', 'ZA:Right pulmonary artery']],
              'protected_lingular_whole_closure': lingular,
              'protected_lingular_rule': 'Preserve all9 coupled pieces plus LPA actual outer33-control seam and all3 reciprocal finite cap interfaces literally. Existing open references are not whole finite arteries. No independent movement of LPA collar, upper branches or lingular ports.',
              'current_recipient_lobe_bindings': target_bindings,
              'current_source_to_recipient_sampled_placement': placement,
              'recipient_winding_samples_sha256': digest(winding_path),
              'current_bronchial_landmark_evidence': str(bronchial_file),
              'current_bronchial_landmark_evidence_sha256': digest(bronchial_file),
              'six_actual_current_bronchial_entries': bronchial['landmarks'],
              'arterial_source_to_named_bronchial_entry_witnesses': bronchial['name_matched_native_arterial_geometric_comparisons'],
              'A7_source_semantic_qualification': bronchial['A7_source_semantic_qualification'],
              'bronchial_landmark_rule': 'Actual current lower-lobar and segmental shared entry cycles are recipient context only. No nearest bronchial point establishes an arterial attachment, source branch identity or target blood lumen.',
              'scope_recommendation': 'The conflicting A7 source has geometric node dependencies on the S6-labeled source, which also meets A9/A10. Review those16 source courses together; A8 has a separate5-course source cohort near the current LPA distal ring. Full LPA/lingular closure and mainPT/bifurcation stay protected. Do not independently tweak the vein around an unregistered named artery, or label an unsupported coarse source basal trunk as all descendants.',
              'new_fit_or_geometry_or_catalogue_change': False, 'anatomical_acceptance': False,
              'script_sha256': digest(__file__), 'elapsed_seconds': time.monotonic() - began}
    output = OUTPUT / 'inventory.json'
    atomic_json(output, report)
    assert digest(CATALOGUE) == EXPECTED
    print(json.dumps({'inventory': str(output), 'sha256': digest(output), 'course_count': len(curves),
                      'source_components': components, 'source_placement': placement,
                      'cross_source_node_pairs': sum(row['cross_source_object'] for row in incidences),
                      'elapsed_seconds': report['elapsed_seconds']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
