"""Create a private static preview, with geometry grouped by system."""
import hashlib
import gzip
import json
import os
from pathlib import Path

root = Path(__file__).resolve().parent.parent / ".local-models"
source = root / 'female-ta98-complete'
target = root / 'female-ta98-pages-by-system-preview'
if target.exists():
    raise SystemExit('Preview already exists; choose a fresh target before staging.')
atlas = json.loads((source / 'atlas.json').read_text())
old_url = '/local-models/female-ta98-complete/'
new_url = '/models/female-ta98/'
chunks = atlas['chunks']
old_banks = [(source / c['url'].removeprefix(old_url)).read_bytes() for c in chunks]
assert all(len(data) == chunks[i]['bytes'] for i, data in enumerate(old_banks))
widths = {'positions': 12, 'normals': 6, 'indices': 4, 'uvs': 8,
          'normalSigns': 1, 'tissueChart': 6, 'gastricLumen': 1, 'colors': 3}
parts = atlas['parts']
target.mkdir()
(target / 'geometry').mkdir()
limit = 50_000_000
banks = []
bank = bytearray()
offsets = {}
system = None
def flush():
    global bank, offsets
    if bank:
        digest = hashlib.sha256(bank).hexdigest()
        path = 'geometry/' + system + '-' + digest[:20] + '.bin'
        (target / path).write_bytes(bank)
        banks.append({'url': new_url + path, 'bytes': len(bank), 'sha256': digest})
        bank = bytearray()
        offsets = {}
for part in sorted(parts, key=lambda p: (p['system'], p['id'])):
    if part['system'] != system:
        flush()
        system = part['system']
    old_index = part['chunk']
    fields = [(key, part[key], (part['indexCount'] if key == 'indices' else part['vertexCount']) * width)
              for key, width in widths.items() if key in part]
    if sum(length + 3 for _, _, length in fields) > limit:
        raise ValueError('One part exceeds the chunk limit: ' + part['id'])
    if len(bank) + sum(length + 3 for _, _, length in fields) > limit:
        flush()
    for key, offset, length in fields:
        old = old_banks[old_index]
        if offset < 0 or offset + length > len(old):
            raise ValueError('Invalid source geometry range: ' + part['id'])
        identity = (old_index, offset, length)
        if identity not in offsets:
            bank.extend(b'\0' * (-len(bank) % 4))
            offsets[identity] = len(bank)
            bank.extend(old[offset:offset + length])
        part[key] = offsets[identity]
    part['chunk'] = len(banks)
flush()
atlas['chunks'] = banks
def relocate(value):
    if isinstance(value, dict): return {key: relocate(item) for key, item in value.items()}
    if isinstance(value, list): return [relocate(item) for item in value]
    if isinstance(value, str) and value.startswith(old_url): return new_url + value[len(old_url):]
    return value
atlas = relocate(atlas)
for path in (source / 'textures').iterdir():
    dst = target / 'textures' / path.name
    dst.parent.mkdir(parents=True, exist_ok=True)
    os.link(path, dst)
for part in atlas['parts']:
    field = part.get('hairField')
    if field:
        path = target / field['url'].removeprefix(new_url)
        with path.open('rb') as source_file, (Path(str(path) + '.gz')).open('wb') as compressed:
            with gzip.GzipFile(filename='', mode='wb', fileobj=compressed, mtime=0) as out:
                out.write(source_file.read())
encoded = json.dumps(atlas, separators=(',', ':'), allow_nan=False)
(target / 'atlas.json').write_text(encoded)
(target / 'atlas-local.json').write_text(encoded.replace(new_url, '/local-models/female-ta98-pages-by-system-preview/'))
print('Local preview:', len(parts), 'parts in', len(banks), 'system chunks; largest:', max(b['bytes'] for b in banks))
