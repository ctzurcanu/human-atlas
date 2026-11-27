"""Bounded full-source liver contact audit; no fitting or installation."""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from geometry_tools.model import AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.distance import exact_intersection_segment

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
WORK = BASE / 'male-coeliac-cystic-composition-v1'
OUT = WORK / 'gastric-liver-source-contact-v12'
CURRENT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
CANDIDATE = WORK / 'candidate-v4/candidate-atlas.json'
OLD = BASE / 'male-coeliac-four-port-v1/binary-candidate-v2/candidate-atlas.json'
HELPER = BASE / 'male-left-gastric-far-wall-helper-v3'
MOVING = 'ZA:Left gastric artery'
FIXED = ['TA98LIVER:male:segment:III', 'TA98PART:ZA:Liver:left', 'TA98PART:ZA:Liver:capsule']


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def key(t):
    return tuple(sorted(map(tuple, t.tolist())))


def distance_to_segment(points, segment):
    a, b = segment
    delta = b - a
    s = np.einsum('ij,j->i', points - a, delta) / np.dot(delta, delta)
    return np.linalg.norm(points - (a + np.clip(s, 0, 1)[:, None] * delta), axis=1)


def main():
    started = time.monotonic()
    OUT.mkdir(exist_ok=False)
    reader = AtlasReader(ROOT)
    atlases = [json.loads(p.read_text()) for p in [CURRENT, CANDIDATE, OLD]]
    rows = [{p['id']: p for p in a['parts']} for a in atlases]
    def mesh(i, pid):
        return reader.read_part(atlases[i], rows[i][pid], frame='model:local-ta98:assembled-metres')
    before, after = mesh(0, MOVING), mesh(1, MOVING)
    assert len(before.faces) == 1440 and len(after.faces) == 2644
    native_lookup = {key(t): i for i, t in enumerate(before.triangles)}
    source = np.load(HELPER / 'surfaces-and-fields.npz')
    no, ni = len(source['outer_faces']), len(source['inner_faces'])
    assert (no, ni) == (1284, 1288)
    ancestry = np.r_[source['outer_source_face'], source['inner_source_face'], np.full(72, -1)]
    role = np.asarray(['outer'] * no + ['inner'] * ni + ['annular'] * 72)
    literal = np.asarray([native_lookup.get(key(t), -1) for t in after.triangles])
    assert int((literal >= 0).sum()) == 1143
    before_bvh, after_bvh = BVH(before), BVH(after)
    reports = []
    for index, pid in enumerate(FIXED):
        fixed, original_fixed = mesh(0, pid), mesh(2, pid)
        assert fixed.geometry_sha == original_fixed.geometry_sha
        assert fixed.attributes_sha == original_fixed.attributes_sha
        tree = BVH(fixed)
        results, summaries = {}, {}
        for name, moving, bvh in [('before', before, before_bvh), ('after', after, after_bvh)]:
            data, summary = pair_contacts(moving, fixed, 0., first_bvh=bvh, second_bvh=tree,
                                          max_pairs=250000, max_tests=2000000)
            assert summary['complete_within_pair_budget']
            mask = np.isin(data['relation'], ['crossing', 'coplanar-overlap'])
            results[name] = data['pairs'][mask]
            summaries[name] = summary
            np.savez_compressed(OUT / f'{index}-{name}-full-witnesses.npz',
                                pairs=results[name], moving_triangles=moving.triangles[results[name][:, 0]],
                                fixed_triangles=fixed.triangles[results[name][:, 1]])
        original_pairs = set(map(tuple, results['before'].tolist()))
        evidence = []
        for face_id, fixed_face_id in results['after']:
            source_face = int(ancestry[face_id])
            literal_id = int(literal[face_id])
            paired = (source_face, int(fixed_face_id)) in original_pairs
            record = dict(after_face_id=int(face_id), fixed_face_id=int(fixed_face_id),
                          role=str(role[face_id]), source_native_face_id=source_face,
                          literal_native_outer_face_id=literal_id,
                          same_source_face_and_same_fixed_triangle_already_crossed=paired,
                          exact_literal_outer_contact_reused=(literal_id, int(fixed_face_id)) in original_pairs)
            if paired and role[face_id] == 'inner':
                old_segment = exact_intersection_segment(before.triangles[source_face], fixed.triangles[fixed_face_id])
                new_segment = exact_intersection_segment(after.triangles[face_id], fixed.triangles[fixed_face_id])
                if len(old_segment['endpoints']) == len(new_segment['endpoints']) == 2:
                    d = distance_to_segment(new_segment['endpoints'], old_segment['endpoints'])
                    record['inner_intersection_endpoint_distance_to_native_outer_segment_m'] = d.tolist()
                    record['old_outer_intersection_segment_m'] = old_segment['endpoints'].tolist()
                    record['new_inner_intersection_segment_m'] = new_segment['endpoints'].tolist()
            evidence.append(record)
        unmatched = [x for x in evidence if not x['exact_literal_outer_contact_reused'] and
                     not x['same_source_face_and_same_fixed_triangle_already_crossed']]
        inner = [x for x in evidence if x['role'] == 'inner']
        distances = [x for q in inner for x in q.get('inner_intersection_endpoint_distance_to_native_outer_segment_m', [])]
        reports.append(dict(fixed_part_id=pid,
                            fixed_geometry_sha256=fixed.geometry_sha,
                            fixed_attributes_sha256=fixed.attributes_sha,
                            fixed_bank_path=fixed.source['bank_path'],
                            fixed_bank_sha256=fixed.source['bank_sha256'],
                            old_and_current_fixed_geometry_and_attributes_exact=True,
                            full_before=summaries['before'], full_after=summaries['after'],
                            full_before_pairs=len(results['before']), full_after_pairs=len(results['after']),
                            after_literal_outer_contacts=sum(q['exact_literal_outer_contact_reused'] for q in evidence),
                            after_inner_contacts=len(inner),
                            inner_contacts_with_same_native_parent_and_fixed_triangle=sum(q['same_source_face_and_same_fixed_triangle_already_crossed'] for q in inner),
                            after_annular_contacts=sum(q['role'] == 'annular' for q in evidence),
                            unexplained_new_support_contacts=unmatched,
                            actual_inner_vs_outer_intersection_endpoint_distance_m=dict(min=min(distances), max=max(distances)) if distances else None,
                            face_evidence=evidence))
    report = dict(schema='full-gastric-liver-source-contact-ancestry/v1',
                  current_catalogue=dict(path=str(CURRENT), sha256=sha(CURRENT)),
                  candidate_catalogue=dict(path=str(CANDIDATE), sha256=sha(CANDIDATE)),
                  old_coeliac_catalogue=dict(path=str(OLD), sha256=sha(OLD)),
                  helper_arrays=dict(path=str(HELPER / 'surfaces-and-fields.npz'), sha256=sha(HELPER / 'surfaces-and-fields.npz')),
                  helper_proof=dict(path=str(HELPER / 'helper-proof.json'), sha256=sha(HELPER / 'helper-proof.json')),
                  native_gastric_geometry_sha256=before.geometry_sha,
                  new_gastric_geometry_sha256=after.geometry_sha,
                  literal_native_outer_faces=1143,
                  full_before_faces=1440, full_after_faces=2644,
                  fixed_parts=reports,
                  all_added_inner_contacts_inherit_same_native_outer_face_and_fixed_triangle=all(not q['unexplained_new_support_contacts'] for q in reports),
                  qualification='An added inner surface can cross liver near an already intersecting native exterior without any changed outer course. These are still actual physical wall/liver contacts, not certified arterial entry, not zero overlap, and not a clearance pass. Exact face-parent correspondence does not make the spatially displaced inner intersection identical to the old outer intersection.',
                  geometry_modified=False, anatomical_acceptance=False, install_allowed=False,
                  elapsed_seconds=time.monotonic() - started)
    path = OUT / 'decision.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(artifact=str(path), sha256=sha(path),
                          reports=[{k: v for k, v in q.items() if k in ['fixed_part_id', 'full_before_pairs', 'full_after_pairs', 'after_literal_outer_contacts', 'after_inner_contacts', 'inner_contacts_with_same_native_parent_and_fixed_triangle', 'unexplained_new_support_contacts', 'actual_inner_vs_outer_intersection_endpoint_distance_m']} for q in reports],
                          elapsed_seconds=report['elapsed_seconds']), indent=2))


if __name__ == '__main__':
    main()
