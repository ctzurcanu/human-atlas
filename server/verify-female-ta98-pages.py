import collections
import gzip
import json
from pathlib import Path

root = Path(__file__).resolve().parent.parent / ".local-models"
source = root / 'female-ta98-complete'
target = root / 'female-ta98-pages-by-system-preview'
a = json.loads((source / 'atlas.json').read_text())
b = json.loads((target / 'atlas.json').read_text())
original = a['parts']
old = [(source / c['url'].removeprefix('/local-models/female-ta98-complete/')).read_bytes() for c in a['chunks']]
new = [(target / c['url'].removeprefix('/models/female-ta98/')).read_bytes() for c in b['chunks']]
assert len(original) == len(b['parts'])
widths = {'positions': 12, 'normals': 6, 'indices': 4, 'uvs': 8,
          'normalSigns': 1, 'tissueChart': 6, 'gastricLumen': 1, 'colors': 3}
systems = collections.defaultdict(set)
for p, q in zip(b['parts'], original):
    assert p['id'] == q['id']
    systems[p['chunk']].add(p['system'])
    for key, width in widths.items():
        assert (key in p) == (key in q)
        if key in p:
            n = (p['indexCount'] if key == 'indices' else p['vertexCount']) * width
            assert new[p['chunk']][p[key]:p[key]+n] == old[q['chunk']][q[key]:q[key]+n], (p['id'], key)
assert all(len(s) == 1 for s in systems.values())
assert all(len(data) < 50_000_000 for data in new)
for p in b['parts']:
    field = p.get('hairField')
    if field:
        path = target / field['url'].removeprefix('/models/female-ta98/')
        assert len(path.read_bytes()) == field['bytes']
        assert gzip.decompress(Path(str(path) + '.gz').read_bytes()) == path.read_bytes()
missing = []
def check(value):
    if isinstance(value, dict):
        for item in value.values(): check(item)
    elif isinstance(value, list):
        for item in value: check(item)
    elif isinstance(value, str) and value.startswith('/models/female-ta98/'):
        if not (target / value.removeprefix('/models/female-ta98/').split('?')[0]).is_file(): missing.append(value)
    elif isinstance(value, str) and value.startswith('/local-models/'):
        missing.append(value)
check(b)
assert not missing, missing[:10]
print('Verified:', len(original), 'parts,', len(new), 'single-system chunks, all geometry bytes, and referenced assets.')
