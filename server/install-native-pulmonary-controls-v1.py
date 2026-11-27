"""Verify existing native identities, then optionally publish their controls locally."""
from pathlib import Path
import copy, hashlib, json, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.local_catalogue_revision import ENC, SHA, install_revision

B = ROOT / '.local-models/geometry-workbench/circulatory-completion/native-pulmonary-controls-v1'
MASTER = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = 'b26860231a10b4b4382802751d6a8c15a7491e8baa83d629d73b93f5abb61707'
PROPOSED = '66b4a76c63043a38639d483dd744faa3789b4e468685e2746f70f738b1c1d11b'

def main():
    start = time.monotonic()
    assert SHA(MASTER) == EXPECTED
    source = json.loads(MASTER.read_bytes())
    candidate_path = B / 'proposed-atlas.json'
    assert SHA(candidate_path) == PROPOSED
    candidate = json.loads(candidate_path.read_bytes())
    assert ENC(candidate) == candidate_path.read_bytes()
    evidence = json.loads((B / 'identity-proof.json').read_bytes())
    rows = evidence['rows']; assert len(rows) == len({r['term'] for r in rows}) == 16
    assert source['parts'] == candidate['parts'] and source['chunks'] == candidate['chunks']
    allowed = {'concepts', 'nativePulmonaryControlIntegration'}
    assert {k for k in source.keys() | candidate.keys() if source.get(k) != candidate.get(k)} == allowed
    old = {c['id']: c for c in source['concepts']}; new = {c['id']: c for c in candidate['concepts']}
    assert old.keys() == new.keys()
    assert {k for k in old if old[k] != new[k]} == {r['concept_id'] for r in rows}
    hierarchy_path = B / 'future-relations-workspace-v1/public/assets/ta98.json'
    hierarchy = {n['id']: n for n in json.loads(hierarchy_path.read_bytes())['nodes']}
    parts = {p['id']: p for p in candidate['parts']}
    for row in rows:
        cid, pid, term = row['concept_id'], row['part_id'], row['term']
        assert new[cid]['ta98Term'] == term and new[cid]['ta98Kind'] == 'exact'
        assert new[cid]['elements'] == [pid] and not parts[pid].get('suppressed')
        assert new[cid]['coverageLimitation'] == row['qualification']
        assert 'unverified' in row['qualification']
        unchanged = dict(new[cid])
        for key in ['ta98Term', 'ta98Kind', 'ta98Method', 'coverageLimitation']:
            unchanged.pop(key)
        assert unchanged == old[cid]
        assert cid in hierarchy[term]['matches']
        assert SHA(row['bank_path']) == row['bank_sha256']
    result = dict(schema='native-pulmonary-control-install-check/v1', passed=True,
                  source_sha256=EXPECTED, proposed_sha256=PROPOSED, terms=[r['term'] for r in rows],
                  unchanged_parts=len(parts), unchanged_geometry=True,
                  hierarchy_sha256=SHA(hierarchy_path), new_triangles=0,
                  finite_wall_or_interface_acceptance=False, elapsed_seconds=time.monotonic()-start)
    if '--install' in sys.argv:
        cacheproof = B / 'publication-caches-ready-v1.json'
        assert cacheproof.is_file(), 'Complete future resolver cache proof is required before publication'
        cache_evidence = json.loads(cacheproof.read_bytes())
        assert cache_evidence['ready'] and cache_evidence['sourceMaleMasterSha256'] == EXPECTED
        assert cache_evidence['futureHierarchySha256'] == SHA(hierarchy_path)
        assert cache_evidence['actualHierarchyBuilderSha256'] == SHA(ROOT / 'scripts/build-ta98-hierarchy.py')
        for row in cache_evidence['rows']:
            assert SHA(ROOT / row['cache']) == row['cacheSha256']
            cache = json.loads((ROOT / row['cache']).read_bytes())
            assert cache['sourceSha256'] == row['sourceSha256']
            for file, digest in cache['resolverInputs'].items():
                actual = (B / 'future-relations-workspace-v1' / file) if file == 'app/data/ta98-parents.json' else ROOT / file
                assert SHA(actual) == digest
        expected_female = SHA(ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json')
        def mutate(sex, atlas, before):
            if sex == 'female':
                assert hashlib.sha256(ENC(before)).hexdigest() == expected_female
                return dict(unchanged=True)
            assert hashlib.sha256(ENC(before)).hexdigest() == EXPECTED
            atlas.clear(); atlas.update(copy.deepcopy(candidate))
            return dict(native_pulmonary_identity_controls=16, new_geometry_bytes=0,
                        source_surfaces_retained=True, geometry_acceptance=False)
        result['models'] = install_revision('male-native-pulmonary-controls-local-v1', mutate)
        assert SHA(MASTER) == PROPOSED
        assert SHA(ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json') == expected_female
        result['installed'] = True
    output = B / ('installed-proof.json' if '--install' in sys.argv else 'controls-and-hierarchy-check.json')
    assert not output.exists()
    output.write_bytes(ENC(result))
    print(json.dumps(result))

if __name__ == '__main__':
    main()
