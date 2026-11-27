"""Bounded current fixed-owner contacts for the exact canonical biarc wall."""
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
from geometry_tools.spatial import BVH

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
OLD = BASE / 'collector-left-common-basal-v1/independent-context-v1'
WALL = BASE / 'biarc-field-finite-wall-v1/candidate-stage-v2/finite-basal-wall.npz'
OUTPUT = WALL.parent / 'independent-fixed-context-v1'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
REPLACED = {'ZA:Inferior basal vein of left lung', 'ZA:Superior basal vein of left lung'}
FRAME = 'model:local-ta98:assembled-metres'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def triangle_key(triangle):
    return tuple(sorted(tuple(point) for point in triangle))


def main():
    began = time.monotonic()
    assert digest(CATALOGUE) == EXPECTED
    assert digest(WALL) == 'e5e65bc043e92de2a179a5ef65f8db4ea6b5ef5c7b0267cb019862cb2f8e810e'
    OUTPUT.mkdir(exist_ok=False)
    atlas = json.loads(CATALOGUE.read_text())
    parts = {part['id']: part for part in atlas['parts']}
    hidden = set(atlas.get('defaultHidden', []))
    old = json.loads((OLD / 'priority-stage.json').read_text())
    remaining = json.loads((OLD / 'remaining-context-stage.json').read_text())
    assert old['source_catalogue_sha256'] == EXPECTED == remaining['source_catalogue_sha256']
    previous = {row['id']: row for row in old['completed_owner_results'] + remaining['completed_owner_results']}
    mesh = Mesh.load(WALL)
    assert mesh.frame == FRAME
    low, high = mesh.vertices.min(axis=0), mesh.vertices.max(axis=0)
    old_bounds = np.asarray(old['candidate_bounds_m'])
    contained = bool(np.all(low >= old_bounds[0]) and np.all(high <= old_bounds[1]))
    previous_owners = {row['id']: row for row in old['aabb_source_owners']}
    owners = dict(previous_owners)
    missing_bounds = []
    for identity, part in parts.items():
        bounds = part.get('bounds')
        if bounds is None:
            if not part.get('suppressed') and not part.get('defaultHidden') and identity not in hidden:
                missing_bounds.append(identity)
            continue
        bounds = np.asarray(bounds, float)
        if not (np.all(bounds[1] >= low) and np.all(bounds[0] <= high)):
            continue
        if identity not in owners:
            owners[identity] = {'id': identity, 'system': part['system'],
                                'record_suppressed': bool(part.get('suppressed')),
                                'record_default_hidden': bool(part.get('defaultHidden')),
                                'catalogue_default_hidden': identity in hidden,
                                'dormant': bool(part.get('suppressed') or part.get('defaultHidden') or identity in hidden),
                                'prospective_replaced_native_source': identity in REPLACED,
                                'added_for_new_candidate_AABB': True}
    active = [row for row in owners.values() if not row['dormant'] and row['id'] not in REPLACED]
    report = {
        'schema': 'canonical-biarc-wall-current-fixed-context-audit/v1',
        'source_catalogue_sha256': EXPECTED, 'candidate_archive_sha256': digest(WALL),
        'candidate_geometry_sha256': mesh.geometry_sha,
        'canonical_proof_sha256': digest(WALL.parent / 'proof.json'),
        'priority_owner_manifest_sha256': digest(OLD / 'priority-stage.json'),
        'remaining_owner_manifest_sha256': digest(OLD / 'remaining-context-stage.json'),
        'script_sha256': digest(__file__), 'candidate_bounds_m': [low.tolist(), high.tolist()],
        'new_AABB_contained_in_old': contained,
        'new_AABB_owner_policy': 'Retain all original active117 minus two explicitly replaced basal sources; append all current metadata-AABB overlapping owners. Current record flags plus catalogue defaultHidden determine activation.',
        'original_fixed_active_owners': sum(not row['dormant'] and row['id'] not in REPLACED for row in previous_owners.values()),
        'added_owner_ids': [identity for identity in owners if identity not in previous_owners],
        'owners': list(owners.values()), 'active_fixed_owners': len(active),
        'active_parts_without_bounds_not_certified': missing_bounds,
        'maximum_seconds_before_next_owner': 32., 'pair_test_budget_each': 200000,
        'pair_result_budget_each': 100000, 'completed_owner_results': [],
        'complete': False, 'occupied_volume_not_certified': True,
        'clinical_or_anatomical_acceptance': False, 'no_geometry_or_catalogue_changes': True,
    }
    output = OUTPUT / 'proof.json'
    write(output, report)
    reader, first_tree = AtlasReader(ROOT), BVH(mesh)
    for ordinal, owner in enumerate(active):
        if time.monotonic() - began >= 32.:
            report['stopped_reason'] = 'Declared32s internal stage bound; no automatic retry.'
            break
        identity = owner['id']
        context = reader.read_part(atlas, parts[identity], frame=FRAME,
                                   catalogue_path=CATALOGUE, catalogue_sha=EXPECTED)
        if identity in previous:
            assert context.geometry_sha == previous[identity]['context_geometry_sha256'], identity
            assert context.source['bank_sha256'] == previous[identity]['bank_sha256'], identity
        contacts, summary = pair_contacts(mesh, context, 0., first_bvh=first_tree,
                                           max_tests=200000, max_pairs=100000)
        relations = contacts['relation']
        pairs = contacts['pairs']
        legitimate = np.zeros(len(pairs), dtype=bool)
        if identity == 'ZA:Left inferior pulmonary vein':
            for index, (first, second) in enumerate(pairs):
                one, two = mesh.triangles[first], context.triangles[second]
                if triangle_key(one) != triangle_key(two):
                    continue
                normal_one = np.cross(one[1] - one[0], one[2] - one[0])
                normal_two = np.cross(two[1] - two[0], two[2] - two[0])
                legitimate[index] = normal_one @ normal_two < 0.
        crossing = (relations == 'crossing') & ~legitimate
        overlap = (relations == 'coplanar-overlap') & ~legitimate
        row = {'id': identity, 'system': owner['system'],
               'context_geometry_sha256': context.geometry_sha,
               'bank_sha256': context.source['bank_sha256'],
               'actual_bounds_m': [context.vertices.min(axis=0).tolist(), context.vertices.max(axis=0).tolist()],
               'actual_triangle_contacts': summary,
               'literal_opposed_shared_annulus_pairs_removed_from_defects': int(legitimate.sum()),
               'remaining_crossing_pairs': int(crossing.sum()),
               'remaining_coplanar_overlap_pairs': int(overlap.sum()),
               'witness_face_pairs': pairs[crossing | overlap][:20].tolist(),
               'occupied_volume_not_certified': True}
        if len(pairs):
            path = OUTPUT / f'owner-{ordinal:03d}-contacts.npz'
            np.savez_compressed(path, **contacts, literal_opposed_shared_annulus_pair=legitimate)
            row['contact_archive'] = str(path)
            row['contact_archive_sha256'] = digest(path)
        report['completed_owner_results'].append(row)
        report['elapsed_seconds'] = time.monotonic() - began
        write(output, report)
        if summary['truncated']:
            report['stopped_reason'] = f'Finite pair budget incomplete for {identity}; no broader budget.'
            break
    else:
        report['complete'] = not missing_bounds
        report['stopped_reason'] = 'All declared fixed owners processed.' if not missing_bounds else 'Some active source records have no bounds; full body scope remains unproven.'
    done = {row['id'] for row in report['completed_owner_results']}
    report['remaining_owner_ids'] = [row['id'] for row in active if row['id'] not in done]
    report['remaining_crossing_total'] = sum(row['remaining_crossing_pairs'] for row in report['completed_owner_results'])
    report['remaining_coplanar_overlap_total'] = sum(row['remaining_coplanar_overlap_pairs'] for row in report['completed_owner_results'])
    report['elapsed_seconds'] = time.monotonic() - began
    assert digest(CATALOGUE) == EXPECTED
    write(output, report)
    print(json.dumps({'proof': str(output), 'sha256': digest(output), 'complete': report['complete'],
                      'owners_done': len(done), 'owners_total': len(active),
                      'new_AABB_owners': len(report['added_owner_ids']),
                      'crossing_total': report['remaining_crossing_total'],
                      'overlap_total': report['remaining_coplanar_overlap_total'],
                      'elapsed_seconds': report['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
