"""Publish verified native blood-vessel identity controls, with geometry literal."""
from pathlib import Path
import copy, json, sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.local_catalogue_revision import ENC, SHA, install_revision

B = ROOT / '.local-models/geometry-workbench/circulatory-completion/native-blood-controls-v2'
MASTER = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
SOURCE_SHA = '66b4a76c63043a38639d483dd744faa3789b4e468685e2746f70f738b1c1d11b'
NEXT_SHA = '043e699fff254b38460d0f008046a11f2291365f57cd84bd2246a8de4c574c54'

def main():
    assert SHA(MASTER) == SOURCE_SHA
    source = json.loads(MASTER.read_bytes()); candidate_path = B / 'proposed-atlas.json'
    assert SHA(candidate_path) == NEXT_SHA
    candidate = json.loads(candidate_path.read_bytes())
    assert ENC(candidate) == candidate_path.read_bytes()
    proof = json.loads((B / 'merge-proof.json').read_bytes())
    assert proof['source_sha256'] == SOURCE_SHA and proof['proposed_sha256'] == NEXT_SHA
    assert source['parts'] == candidate['parts'] and source['chunks'] == candidate['chunks']
    changed_top = {k for k in source.keys() | candidate.keys() if source.get(k) != candidate.get(k)}
    assert changed_top == {'concepts', 'nativeNamedArterialControlIntegration', 'nativeLingularVenousControlIntegration'}
    old = {c['id']: c for c in source['concepts']}; new = {c['id']: c for c in candidate['concepts']}
    assert old.keys() == new.keys()
    delta = {cid for cid in old if old[cid] != new[cid]}
    assert delta == set(proof['changed_concepts']) and len(delta) == 21
    assert len({new[cid]['ta98Term'] for cid in delta}) == 11
    hierarchy_path = B / 'future-relations-workspace-v1/public/assets/ta98.json'
    hierarchy = {n['id']: n for n in json.loads(hierarchy_path.read_bytes())['nodes']}
    parts = {p['id']: p for p in candidate['parts']}; hidden = set(candidate.get('defaultHidden', []))
    for cid in delta:
        item = new[cid]; term = item['ta98Term']
        assert item['ta98Kind'] == 'exact' and item.get('coverageLimitation')
        assert cid in hierarchy[term]['matches']
        assert item['elements'] == old[cid]['elements']
        assert item['elements'] and all(p in parts and not parts[p].get('suppressed') and p not in hidden for p in item['elements'])
        before = dict(item)
        for key in ['ta98Term', 'ta98Kind', 'ta98Method', 'coverageLimitation']:
            before.pop(key, None)
        assert before == old[cid]
    result = dict(schema='native-blood-controls-publication/v2', passed=True,
                  source_sha256=SOURCE_SHA, proposed_sha256=NEXT_SHA,
                  control_concepts=21, term_count=11, terms=proof['terms'],
                  all_geometry_material_visibility_literal=True, new_triangles=0,
                  geometry_or_interface_acceptance=False, hierarchy_sha256=SHA(hierarchy_path))
    if '--install' in sys.argv:
        cacheproof = json.loads((B / 'publication-caches-ready-v1.json').read_bytes())
        assert cacheproof['ready'] and cacheproof['futureHierarchySha256'] == SHA(hierarchy_path)
        for row in cacheproof['rows']:
            path = ROOT / row['cache']; assert SHA(path) == row['cacheSha256']
            cache = json.loads(path.read_bytes()); assert cache['sourceSha256'] == row['sourceSha256']
            for file, digest in cache['resolverInputs'].items():
                actual = B / 'future-relations-workspace-v1' / file if file == 'app/data/ta98-parents.json' else ROOT / file
                assert SHA(actual) == digest
        female = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json'; female_sha = SHA(female)
        def mutate(sex, atlas, before):
            if sex == 'female':
                return dict(unchanged=True)
            assert ENC(before) == MASTER.read_bytes() and SHA(MASTER) == SOURCE_SHA
            atlas.clear(); atlas.update(copy.deepcopy(candidate))
            return dict(native_control_concepts=21, distinct_terms=11, geometry_literal=True,
                        walls_lumens_packing_and_recipient_fit_unverified=True)
        result['models'] = install_revision('male-native-blood-controls-local-v2', mutate)
        assert SHA(MASTER) == NEXT_SHA and SHA(female) == female_sha
        result['installed'] = True
    output = B / ('installed-proof.json' if '--install' in sys.argv else 'identity-and-hierarchy-proof.json')
    assert not output.exists(); output.write_bytes(ENC(result)); print(json.dumps(result))

if __name__ == '__main__':
    main()
