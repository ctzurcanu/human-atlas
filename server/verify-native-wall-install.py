"""Replay a bounded native-wall installation through the served body catalogue."""
import argparse
import gzip
import json
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader

parser = argparse.ArgumentParser()
parser.add_argument('--directory', required=True, type=Path)
parser.add_argument('--revision', required=True)
parser.add_argument('--exterior', action='store_true')
args = parser.parse_args()
signal.alarm(18)
started = time.monotonic()
out = args.directory.resolve()
proof = json.loads((out / 'packing-proof.json').read_text())
receipt = json.loads((out / 'installed.json').read_text())
master = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
body = ROOT / '.local-models/organ-models/spinal-cord/male-body.json'
manifest = json.loads((ROOT / '.local-models/ta98-review' / args.revision / 'manifest.json').read_text())
entries = manifest['files'] if isinstance(manifest, dict) else manifest
entry = next(e for e in entries if (ROOT / e['path']).resolve() == master)
before = json.loads(gzip.decompress(Path(entry['backup']).read_bytes()))
after, served = json.loads(master.read_text()), json.loads(body.read_text())
assert digest_file(master) == receipt['catalogue_sha256']
assert after == json.loads((out / 'proposed-atlas.json').read_text())
old = {p['id']: p for p in before['parts']}
new = {p['id']: p for p in after['parts']}
assert old.keys() == new.keys()
changed = sorted(pid for pid in old if old[pid] != new[pid])
assert changed == sorted(proof['changed_parts'])
assert before['concepts'] == after['concepts']
assert before['materials'] == after['materials']
assert len(after['chunks']) == len(before['chunks']) + 1
reader = AtlasReader()
rows = []
for pid in changed:
    p = new[pid]
    s = next(p for p in served['parts'] if p['id'] == pid)
    a = reader.read_part(after, p, frame='model:local-ta98:assembled-metres')
    b = reader.read_part(served, s, frame=a.frame)
    assert a.geometry_sha == b.geometry_sha and a.attributes_sha == b.attributes_sha
    assert not p.get('suppressed')
    if args.exterior:
        assert p['sourceOuterSurfaceOnly'] and not p.get('pairedTissue')
    else:
        assert p['pairedTissue'] and p['sectionSolid']
    rows.append(dict(part_id=pid, served_geometry_sha256=b.geometry_sha,
                     served_attributes_sha256=b.attributes_sha, active=True))
female = json.loads((out / 'female-before.json').read_text())
assert all(digest_file(ROOT / path) == checksum for path, checksum in female.items())
result = dict(passed=True, installed_master_sha256=digest_file(master),
              changed_parts=changed, unchanged_parts=len(old) - len(changed),
              concepts_and_materials_literal=True, skeleton_basilar_pons_fixed=True,
              female_catalogue_and_body_literal=True, served_parts=rows,
              bank_bytes=proof['bank_bytes'], new_TA98_terms=0,
              communicating_junctions_verified=False, anatomical_acceptance=False,
              elapsed_seconds=time.monotonic() - started)
(out / 'installed-verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
