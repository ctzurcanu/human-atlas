"""Source-bound ownership/install plan and two read-only wall-section probes."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.simple_compartment_section import polygon_contains

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
OUTPUT = BASE / 'basal-three-owner-install-plan-v1'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FIELDS = BASE / 'whole-unbranched-root-stem-biarc-v1/proper-profile-transport-v1/course-fields.npz'
WALL = BASE / 'biarc-field-finite-wall-v1/candidate-stage-v2/finite-basal-wall.npz'
WALL_PROOF = WALL.parent / 'proof.json'
NEW_ID = 'TA98VENOUS:male:left-common-basal-vein:finite-wall'
PARENT_ID = 'ZA:Left inferior pulmonary vein'
SUPERIOR_ID = 'ZA:Superior basal vein of left lung'
INFERIOR_ID = 'ZA:Inferior basal vein of left lung'
V6_ID = 'ZA:Superior vein of left lung'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def unit(vector):
    return vector / np.linalg.norm(vector)


def section_probe(mesh, center, normal):
    scalar = (mesh.vertices - center) @ normal
    classified = np.abs(scalar) <= 2e-13
    scalar[classified] = 0.
    try:
        clipped = clip_scalar_surface(mesh.vertices, mesh.faces, scalar, 0.)
        loops = boundary_loops(clipped['negative_faces'])
        axis = np.eye(3)[np.argmin(np.abs(normal))]
        horizontal = unit(np.cross(normal, axis))
        basis = np.column_stack([horizontal, np.cross(normal, horizontal)])
        rows = []
        for ring in loops:
            world = clipped['positions'][ring]
            xy = (world - center) @ basis
            area = .5 * float(np.sum(xy[:, 0] * np.roll(xy[:, 1], -1) - xy[:, 1] * np.roll(xy[:, 0], -1)))
            rows.append({'controls': len(ring), 'signed_projected_area_m2': area,
                         'encloses_source_section_center': polygon_contains(xy, np.zeros(2)),
                         'maximum_center_distance_m': float(np.linalg.norm(world - center, axis=1).max())})
        return {'intersection_succeeded': True, 'contour_count': len(rows), 'contours': rows,
                'simple_two_contour_annular_section_only': len(rows) == 2 and all(row['encloses_source_section_center'] for row in rows),
                'source_controls_classified_on_plane': int(classified.sum()),
                'caps_emitted': 0, 'partition_emitted': False}
    except ValueError as error:
        return {'intersection_succeeded': False, 'error': str(error), 'partition_emitted': False}


def main():
    started = time.monotonic()
    assert digest(CATALOGUE) == '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
    assert digest(FIELDS) == '858aff5f40d8f496949d860f43db4e18b8897346079782f887e8ab3336446c96'
    assert digest(WALL) == 'e5e65bc043e92de2a179a5ef65f8db4ea6b5ef5c7b0267cb019862cb2f8e810e'
    OUTPUT.mkdir(exist_ok=False)
    atlas = json.loads(CATALOGUE.read_text())
    parts = {part['id']: part for part in atlas['parts']}
    concepts = {concept['id']: concept for concept in atlas['concepts']}
    inventory_file = BASE / 'venous-inventory-v2.json'
    inventory = json.loads(inventory_file.read_text())
    terms = [{ 'term': row['authority']['id'], 'name': row['authority']['name_en'],
               'parent': row['authority']['parent_id'], 'classification': row['classification'],
               'current_status': row['current']['representation_status'],
               'current_source_ids': row['current']['resolved_direct_part_ids'],
               'source_candidates_only': row['proposed_source_ids'], 'accepted': False}
             for row in inventory['terms'] if row['authority']['id'].startswith('A12.3.02.4')]
    assert len(terms) == 10
    sources = [{key: row[key] for key in ['id', 'geometry_sha256', 'bank_sha256', 'attributes_sha256', 'record_sha256']}
               for row in inventory['sources'] if row['id'] in [PARENT_ID, SUPERIOR_ID, INFERIOR_ID, V6_ID]]
    current = [{ 'part': parts[identity], 'concept': concepts[identity]} for identity in [PARENT_ID, SUPERIOR_ID, INFERIOR_ID, V6_ID]]
    with np.load(FIELDS, allow_pickle=False) as archive:
        fields = {key: archive[key].copy() for key in archive.files}
    slices = {int(curve): (int(first), int(last)) for curve, first, last in fields['curve_slices']}
    mesh = Mesh.load(WALL)
    probes = []
    for curve, station, owner in [(131, 5, SUPERIOR_ID), (60, 8, INFERIOR_ID)]:
        begin, end = slices[curve]
        native = fields['all_source_centres'][begin:end]
        source_tangent = unit(np.gradient(native, axis=0, edge_order=2)[station])
        normal = unit(fields['all_station_rotations'][begin + station] @ source_tangent)
        center = fields['all_centres'][begin + station]
        probes.append({'curve': curve, 'station': station, 'distal_owner_id': owner,
                       'center_m': center.tolist(), 'normal_toward_distal': normal.tolist(),
                       'actual_outer_profile': fields['all_outer'][begin + station].tolist(),
                       'actual_inner_profile': fields['all_inner'][begin + station].tolist(),
                       'read_only_full_wall_plane_probe': section_probe(mesh, center, normal),
                       'ownership_validated': False,
                       'qualification': 'This is a candidate anatomical confluence boundary. Plane contours alone do not establish correct local cut ownership, nonoverlap or branch incidence.'})
    plan = {
        'schema': 'source-bound-left-basal-three-owner-install-plan/v1',
        'state': 'Plan only; no partition, catalogue, app, model or database mutation authorized.',
        'source_catalogue_sha256': digest(CATALOGUE), 'inventory_sha256': digest(inventory_file),
        'authority_cache_sha256': inventory['authority_cache_sha256'],
        'inclusive_401_subtree': terms, 'inclusive_terms': 10,
        'currently_missing_terms': [row['term'] for row in terms if row['current_status'] == 'missing-direct-representation'],
        'current_sources': sources, 'current_parts_and_concepts': current,
        'candidate_wall': str(WALL), 'candidate_wall_sha256': digest(WALL),
        'candidate_wall_proof_sha256': digest(WALL_PROOF), 'source_fields_sha256': digest(FIELDS),
        'candidate_cut_probes': probes,
        'physical_ownership': [
            {'term': 'A12.3.02.405', 'proposed_part_id': NEW_ID,
             'owns': 'Finite proximal common-basal stem and communicating superior/inferior basal confluence between literal current LIPV basal annulus and validated daughter cut cycles.',
             'excludes': 'All distal eight-course-tree tissue beyond the two declared cuts; V6; atrial/heart root; unsupported intra/intersegmental or wall-coat labels.',
             'concept_elements': [NEW_ID]},
            {'term': 'A12.3.02.406', 'reuse_part_id': SUPERIOR_ID,
             'owns': 'Distal curve131 beyond a validated confluence cut near profile5, plus its downstream132/133/134 tissue and distal mouths, each exclusively.'},
            {'term': 'A12.3.02.410', 'reuse_part_id': INFERIOR_ID,
             'owns': 'Distal curve60 from validated profile8 boundary, plus downstream61/62/63 tissue and distal mouths, each exclusively.'},
        ],
        'partition_requirements': [
            'Finish complete canonical wall/context/render review before partitioning. Existing G8/self0 proof is preliminary, not install acceptance.',
            'Validate actual posed local sections with one outer and one inner closed cycle, a usable annular tissue patch, disjoint cut surfaces and correct daughter-side/source ancestry. Multiple full-plane contours require local component ownership; do not cut unrelated distal branches.',
            'Partition one coherent finite wall into exactly three nonoverlapping owned cells. Source boundary real faces belong to exactly one cell; no stacked duplicate surfaces.',
            'Store each internal annular tissue cap once in one common geometry bank and use it reciprocally with opposite oriented faces/normal signs and identical normalword/UV charts. A cap is shared tissue interface, never a lumen disk or separate anatomical solid.',
            'Prove union equals repaired parent wall, pairwise occupied overlap zero, three cells closed/correct genus/self0 and open communicating lumen/blood paths. Preserve all8 distal mouths and7 source branch incidences.',
        ],
        'fixed_parent_and_others': {
            'LIPV_id': PARENT_ID, 'all_current_5960_faces_and_attributes_literal': True,
            'heart_seam_source_faces': [0, 3585], 'heart_seam_kind': 'Three nonplanar source-owned endocardium/epicardium/LA myocardium patches; no fabricated atrial plane.',
            'V6_interface_and_402_source_unchanged': True, 'skeleton_heart_lungs_bronchi_other_vessels_unchanged': True,
            'parent_metadata_only_allowed_after_review': 'If a shared-bank repack is needed, preserve exact current parent storage/pose and all source fields. Split only existing free draw-group range to label the exact24-face basal annulus neighbor; retain all heart groups and V6 fields.',
        },
        'minimal_catalogue_delta_after_all_gates': {
            'parts': 'Add405 finite cell; replace geometry of existing406 and410 IDs with their exclusive finite distal cells; retain current401 ID and all source geometry/heart/V6 fields. Repack current401 only if required for literal reciprocal shared storage, with decoded full-field equality proof.',
            'new405_metadata': {'id': NEW_ID, 'conceptId': NEW_ID, 'ta98Term': 'A12.3.02.405',
                               'name': 'Common basal vein of left lung', 'system': 'venous',
                               'sourceId': 'Constructed left common basal pulmonary vein tissue region',
                               'pairedTissue': True, 'sectionSolid': True,
                               'sectionAssembly': 'left-basal-pulmonary-venous-finite-wall:male',
                               'partOf': PARENT_ID, 'approved_bank_and_frame_marker': 'Must bind actual reviewed final bank SHA/frame; unknown now, never a provisional wildcard.'},
            '406_and410_metadata': 'Keep IDs/source provenance; explicit exact ta98Term406/410, conceptId own ID, partOf405concept, pairedTissue/sectionSolid and normalSigns/UV/draw-group reciprocal neighbors. Do not bind all three terms to the same complete wall.',
            'concept401': {'id': PARENT_ID, 'elements': [PARENT_ID, V6_ID, NEW_ID, SUPERIOR_ID, INFERIOR_ID],
                           'componentElements': [PARENT_ID, V6_ID, NEW_ID, SUPERIOR_ID, INFERIOR_ID],
                           'reason': 'Truthful whole LIPV source assembly; retains root when components are exposed. Existing402 source identity remains qualified and unchanged.'},
            'concept405': {'id': NEW_ID, 'ta98Term': 'A12.3.02.405', 'ta98Kind': 'exact', 'elements': [NEW_ID], 'partOf': PARENT_ID},
            'explicit_ta98_bindings': 'Add exact405 binding to405cell only; retain or make explicit406/410 exact bindings to respective exclusive cells. All IDs, terms and concept elements must match; no hidden or suppressed carriers.',
            'ordinary_source_collections_to_append_405': ["ZA:collection:Left inferior pulmonary vein'", 'ZA:collection:Pulmonary veins', 'ZA:collection:Trunk', 'ZA:collection:5: Cardiovascular system'],
            'derived_hierarchy': 'Rebuild source/code-bound hierarchy and relation caches from final reviewed proposal. Preserve Heart whole concept membership and402 record; regenerate derived TA98 memberships instead of manual aliases.',
            'unchanged_missing_descendants': ['A12.3.02.403', 'A12.3.02.404', 'A12.3.02.407', 'A12.3.02.408', 'A12.3.02.409'],
        },
        'minimal_app_controls_after_final_bank_review': {
            'Systems': 'app/anatomy-hierarchy.ts:hierarchyEntries requires a physical405 part grouped by conceptId; a concept-only alias is insufficient. Keep venous system and current Trunk/Pulmonary veins/source group membership.',
            'term_binding': 'app/ta98-term-bindings.ts:catalogueTa98TermParts needs exact matching concept/part term IDs, unique actual elements, sectionSolid and active visible records.',
            'parent': 'app/ta98-modeled-parent.ts:ta98ModeledParent: whole401 explicit containment supports405→401.406/410→finite405 needs a narrow audited family/source/bank/frame guard recognizing405 nativeZA frame; partOf alone cannot bypass parent containment. Do not widen405.elements to full tree or game an ID prefix.',
            'components': 'app/anatomical-components.ts and ta98ModeledComponentIds need a narrow reviewed-family rule analogous to finite lingual trunk: Show components for405 returns405ownedcell plus406/410, preserving confluence tissue. Apply to ordinary concept and TA98 node without claiming405direct geometry is the entire tree. Do not use componentElements outside concept.elements.',
            'interfaces': 'Reciprocal tissueDrawGroups neighbors plus pairedTissue/sectionSolid/normalSigns/UV; scene.tsx and tissue-interface-visibility.ts suppress shared tissue only when peers visible and poses coincide. No opening-as-solid or tunica-coat controls.',
            'target_tests': ['Systems and TA98 hierarchy each expose one distinct405control and separate406/410 controls.', 'Selecting405 direct geometry selects only confluencecell; Show components retains405 plus both distalowners.', '406/410 parent navigates405;405 parent navigates401; left/right and bank/frame tampering rejected.', 'Whole401 components retain current root,V6 and all basalowners without duplicates; all other term/source memberships unchanged.', 'Nine external port field replay, shared-interface visibility, section patency and camera/explosion invariants.'],
        },
        'installation_blocked_until': 'Full canonical geometry, nineports/patency, current fixed-context, render/cross-section, all source-field and packed-bank replay, final exact family guard/tests, source-bound job review and explicit install authorization.',
        'new_distinct_term_if_later_reviewed': 'At most405 only for this three-owner construction; no completed-subtree or anatomy acceptance claim.',
        'no_partition_catalogue_app_model_DB_edits': True, 'anatomical_acceptance': False,
        'elapsed_seconds': time.monotonic() - started,
    }
    write(OUTPUT / 'plan.json', plan)
    with (ROOT / 'reports/mesh-audit/REVIEW.md').open('a') as stream:
        stream.write('\n\nLeft common-basal ownership/install plan (source346bed): .405 is only finite communicating confluence/stem; .406/.410 exclusively retain their distal courses and branches. All10 .401 descendants inventoried;403/404/407/408/409 remain unresolved. Candidate131profile5 and60profile8 cuts are read-only section probes, not accepted partitions. Literal currentLIPV heart patches/V6 and all other anatomy protected. Ordinary/TA98 parent+component controls require a narrow actualfamily guard, never whole-tree405 aliases. No partition/atlas/app/DB edits. Machine plan: `' + str((OUTPUT / 'plan.json').relative_to(ROOT)) + '`; SHA`' + digest(OUTPUT / 'plan.json') + '`.\n')
    assert digest(CATALOGUE) == plan['source_catalogue_sha256'] and digest(WALL) == plan['candidate_wall_sha256']
    print(json.dumps({'plan': str(OUTPUT / 'plan.json'), 'sha256': digest(OUTPUT / 'plan.json'),
                      'probes': [{'curve': row['curve'], 'station': row['station'],
                                  'probe': row['read_only_full_wall_plane_probe']} for row in probes],
                      'elapsed_seconds': plan['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
