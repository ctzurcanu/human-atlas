"""One source-bound rigid root correction, without a fit or parameter search.

This is an exterior-course candidate. Crossing a finite parent wall is NOT a
communicating lumen; no catalogue is modified by this script.
"""
from pathlib import Path
import json
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.audit import audit
from geometry_tools.common import ROOT, atomic_json, atomic_npz, digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.algorithms.axial_endpoint_field import endpoint_translation

P = ROOT / '.local-models/geometry-workbench/circulatory-completion/pancreatic-public-cohort-v1'
SOURCE = P / 'source-stage-v1/registration-v2'
OUT = P / 'greater-root-rigid-correction-v1'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '0fa165fa0b5a769242b8c406d34df291c2ce7080bc32684e8679611ee61c41b2'
FRAME = 'model:local-ta98:assembled-metres'


def main():
    started = time.monotonic()
    endpoint = '--endpoint-field' in sys.argv[1:]
    global OUT
    if endpoint:
        OUT = P / 'greater-root-endpoint-correction-v2'
    assert digest_file(CATALOGUE) == EXPECTED
    assert not OUT.exists()
    OUT.mkdir()
    lineage_path = SOURCE / 'splenic-root-witnesses-v1/greater-current-wall-lineage.json'
    # Independently measured line: near outer .306577mm, near inner .492669mm,
    # far inner 3.222004mm. A single .8mm translation retains the entire cast.
    direction = np.array([.3616935483, .010940924, .9322328429])
    direction /= np.linalg.norm(direction)
    delta = .0008 * direction
    plan = dict(
        schema='pancreatic-greater-root-rigid-correction/v1',
        source_catalogue_sha256=EXPECTED,
        inventory=str(P / 'current-inventory.json'),
        inventory_sha256=digest_file(P / 'current-inventory.json'),
        source_lineage=str(lineage_path), source_lineage_sha256=digest_file(lineage_path),
        changed_candidate='FJ3433', moving_field='quintic endpoint field' if endpoint else 'one proper rigid translation',
        translation_metres=delta.tolist(), maximum_translation_metres=.001,
        parent_and_all_other_candidate_geometry_fixed=True,
        all_organs_skeleton_heart_basilar_pons_fixed=True,
        calibre_and_shape_must_remain_literal=not endpoint,
        parallel_section_profiles_must_remain_literal=endpoint,
        parameters_searched=0,
        interfaces=[
            'Greater exterior must reach the actual native parent wall.',
            'The existing greater-to-caudal source contact must survive.',
            'No pancreatic, accessory pancreatic or bile-duct crossing.',
            'Exterior overlap does not establish an arterial lumen or ostium.',
        ],
        anatomical_acceptance=False, install_authorized_by_this_script=False,
    )
    z = np.load(SOURCE / 'FJ3433-rigid-candidate.npz')
    old_vertices = z['vertices'].copy()
    arrays = {key: z[key].copy() for key in z.files}
    field_proof = None
    if endpoint:
        _, _, axes = np.linalg.svd(old_vertices - old_vertices.mean(0), full_matrices=False)
        axis = axes[0]
        if axis[1] < 0:
            axis = -axis
        low, high = (old_vertices @ axis).min(), (old_vertices @ axis).max()
        start, stop = low + .45 * (high - low), low + .97 * (high - low)
        plan.update(schema='pancreatic-greater-root-endpoint-correction/v2',
                    axis=axis.tolist(), field_start_metres=float(start), field_stop_metres=float(stop),
                    fixed_tail_region='All controls at or below45% of the actual posed principal-axis extent',
                    supersedes_failed_rigid_plan_sha256=digest_file(P / 'greater-root-rigid-correction-v1/plan.json'))
        atomic_json(OUT / 'plan.json', plan)
        vertices, normals, field_proof = endpoint_translation(
            old_vertices, z['normals_world_float'], axis, start, stop, delta)
        arrays['normals_world_float'] = normals
        arrays['endpoint_field_axis'] = axis
        arrays['endpoint_field_planes'] = np.array([start, stop])
    else:
        atomic_json(OUT / 'plan.json', plan)
        vertices = old_vertices + delta
    arrays['vertices'] = vertices
    arrays['post_registration_translation'] = delta
    atomic_npz(OUT / 'FJ3433-corrected.npz', **arrays)
    mesh = Mesh(vertices, z['faces'], FRAME)
    mesh.save(OUT / 'greater-preview.npz')
    shape_exact = np.max(np.abs((vertices - vertices[0]) - (old_vertices - old_vertices[0])))
    checks = audit(mesh, require_closed=True, expected_components=1, expected_genus=[0],
                   max_pairs=250000, max_tests=2000000)
    atomic_json(OUT / 'greater-audit.json', checks)
    catalogue = json.loads(CATALOGUE.read_text())
    parts = {row['id']: row for row in catalogue['parts']}
    reader = AtlasReader()
    rows = []
    targets = [('native-splenic', 'TA98PART:ZA:Splenic artery:trunk'),
               ('main-pancreatic-duct', 'ZA:Pancreatic duct'),
               ('accessory-pancreatic-duct', 'ZA:Accessory pancreatic duct')]
    # Select only the previously audited exact native bile-duct IDs.
    previous = json.loads((SOURCE / 'duct-context-terminal-v2.json').read_text())
    atomic_json(OUT / 'prior-duct-report-binding.json',
                dict(path=str(SOURCE / 'duct-context-terminal-v2.json'),
                     sha256=digest_file(SOURCE / 'duct-context-terminal-v2.json'),
                     top_level_keys=sorted(previous)))
    for name, pid in targets:
        other = reader.read_part(catalogue, parts[pid], frame=FRAME, catalogue_path=CATALOGUE)
        contacts, summary = pair_contacts(mesh, other, 0., max_pairs=250000, max_tests=2000000)
        atomic_npz(OUT / (name + '-contacts.npz'), **contacts)
        assert summary['complete_within_pair_budget']
        rows.append(dict(name=name, part_id=pid, geometry_sha256=other.geometry_sha, **summary))
    for pid in ['FJ3430', 'FJ3444', 'FJ3419']:
        q = np.load(SOURCE / (pid + '-rigid-candidate.npz'))
        other = Mesh(q['vertices'], q['faces'], FRAME)
        contacts, summary = pair_contacts(mesh, other, 0., max_pairs=250000, max_tests=2000000)
        atomic_npz(OUT / (pid + '-contacts.npz'), **contacts)
        assert summary['complete_within_pair_budget']
        rows.append(dict(name=pid, geometry_sha256=other.geometry_sha, **summary))
    get = lambda name: next(r for r in rows if r['name'] == name)
    passed = (checks['passed'] and (field_proof is not None or shape_exact < 5e-16)
              and get('native-splenic')['crossing_face_pairs'] > 0
              and get('FJ3419')['crossing_face_pairs'] > 0
              and all(not get(name)['crossing_face_pairs'] and not get(name)['coplanar_overlap_face_pairs']
                      for name in ['main-pancreatic-duct', 'accessory-pancreatic-duct']))
    proof = dict(schema=plan['schema'], source_catalogue_sha256=EXPECTED,
                 plan_sha256=digest_file(OUT / 'plan.json'),
                 corrected_npz_sha256=digest_file(OUT / 'FJ3433-corrected.npz'),
                 translation_metres=delta.tolist(), relative_shape_max_error_metres=float(shape_exact),
                 normals_unchanged=not endpoint, face_order_unchanged=True, rows=rows,
                 endpoint_field=field_proof,
                 exterior_course_gates_passed=bool(passed),
                 remaining=['Full bile-duct recheck', 'Explicit communicating finite-wall ostia',
                            'Whole fixed-context and actual packed-normal review'],
                 anatomical_acceptance=False, installed=False,
                 elapsed_seconds=time.monotonic() - started)
    assert digest_file(CATALOGUE) == EXPECTED
    atomic_json(OUT / 'proof.json', proof)
    print(json.dumps(dict(passed=bool(passed), displacement_mm=.8,
                         parent_crossing_pairs=get('native-splenic')['crossing_face_pairs'],
                         caudal_crossing_pairs=get('FJ3419')['crossing_face_pairs'],
                         elapsed_seconds=proof['elapsed_seconds'])))


if __name__ == '__main__':
    main()
