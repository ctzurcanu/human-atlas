"""Verify local publication, served cells and exact original source retention."""
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

signal.alarm(18)
start = time.monotonic()
work = ROOT / '.local-models/geometry-workbench/circulatory-completion/splenic-owned-branches-v1/candidate-v1'
receipt = json.loads((work / 'installed.json').read_text())
master = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
body = ROOT / '.local-models/organ-models/spinal-cord/male-body.json'
manifest = json.loads((ROOT / '.local-models/ta98-review/male-splenic-owned-branches-local-v1/manifest.json').read_text())
entries = manifest['files'] if isinstance(manifest, dict) else manifest
oldfile = next(e for e in entries if (ROOT / e['path']).resolve() == master)
before = json.loads(gzip.decompress((ROOT / oldfile['backup']).read_bytes()))
after, served = json.loads(master.read_text()), json.loads(body.read_text())
assert digest_file(master) == receipt['catalogue_sha256']
assert after == json.loads((work / 'proposed-atlas.json').read_text())
old, new = ({p['id']: p for p in a['parts']} for a in [before, after])
additions = sorted(set(new) - set(old))
changed = sorted(pid for pid in old if old[pid] != new[pid])
assert changed == ['ZA:Splenic artery']
assert additions == ['TA98PART:ZA:Splenic artery:A12.2.12.051', 'TA98PART:ZA:Splenic artery:trunk']
reader = AtlasReader()
frame = 'model:local-ta98:assembled-metres'
a = reader.read_part(before, old[changed[0]], frame=frame)
b = reader.read_part(after, new[changed[0]], frame=frame)
assert a.geometry_sha == b.geometry_sha and a.attributes_sha == b.attributes_sha
assert new[changed[0]]['suppressed'] and new[changed[0]]['defaultHidden']
rows = []
for pid in additions:
    source = reader.read_part(after, new[pid], frame=frame)
    record = next(p for p in served['parts'] if p['id'] == pid)
    actual = reader.read_part(served, record, frame=frame)
    assert actual.geometry_sha == source.geometry_sha and actual.attributes_sha == source.attributes_sha
    assert not record.get('suppressed') and record['pairedTissue'] and record['sectionSolid']
    rows.append(dict(part_id=pid, geometry_sha256=actual.geometry_sha, attributes_sha256=actual.attributes_sha))
female = json.loads((ROOT / '.local-models/geometry-workbench/circulatory-completion/foregut-arterial-batch-v4/native-walls-v1/female-before.json').read_text())
assert all(digest_file(ROOT / path) == sha for path, sha in female.items())
result = dict(passed=True, installed_master_sha256=digest_file(master),
              new_parts=additions, changed_old_parts=changed, unchanged_old_parts=len(old) - len(changed),
              original_source_retained_literal_and_hidden=True,
              skeleton_basilar_pons_and_all_other_geometry_fixed=True,
              female_master_and_body_literal=True, served_meshes=rows,
              anatomical_acceptance=False, hilar_perfusion_verified=False,
              elapsed_seconds=time.monotonic() - start)
(work / 'installed-verification.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
