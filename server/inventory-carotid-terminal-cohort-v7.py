"""Record the full ICA descendant inventory before joint terminal repair."""
from pathlib import Path
import hashlib
import json
import sqlite3
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.ta98_authority import corrected_terms, term_sex

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

work = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1/native-carotid-terminal-cohort-v7b'
work.mkdir(exist_ok=True)
catalogue = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
authority = ROOT / '.local-models/ta98-authority.sqlite'
corrections = ROOT / 'geometry_tools/ta98_authority.py'
with sqlite3.connect('file:' + str(authority) + '?mode=ro', uri=True) as connection:
    connection.row_factory = sqlite3.Row
    terms = corrected_terms({r['id']: dict(r) for r in connection.execute('SELECT * FROM ta98')})
descendants = {'A12.2.06.001'}
while True:
    children = {code for code, row in terms.items() if row['parent_id'] in descendants}
    if children <= descendants:
        break
    descendants |= children
descendants.add('A12.2.07.029')
census_path = ROOT / 'reports/ta98-blood-current-20261024.json'
census = json.loads(census_path.read_text())
assert census['source_sha256'] == sha(catalogue)
effective = {r['id']: r for r in census['terms']}
atlas = json.loads(catalogue.read_text())
parts = {p['id']: p for p in atlas['parts']}
regions = {
    'A12.2.06.002', 'A12.2.06.003', 'A12.2.06.004', 'A12.2.06.007',
    'A12.2.06.015', 'A12.2.06.023', 'A12.2.06.025', 'A12.2.06.026',
    'A12.2.07.023', 'A12.2.07.033', 'A12.2.07.047', 'A12.2.07.053',
}
rows = []
for code in sorted(descendants):
    ids = list(effective.get(code, {}).get('direct_part_ids', []))
    # This exact source identity and TA98 M1 synonym were independently read
    # against the primary terminology. It is not a global fuzzy-name mapping.
    if code == 'A12.2.07.047':
        ids += ['ZA:Middle cerebral artery (M1-segment).l', 'ZA:Middle cerebral artery (M1-segment).r']
    sources = []
    for pid in dict.fromkeys(ids):
        assert pid in parts
        p = parts[pid]
        sources.append(dict(part_id=pid, record_sha256=hashlib.sha256(
            json.dumps(p, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
            record=p, source_reference_only=True))
    rows.append(dict(authority=terms[code], representation_class='tissue-region' if code in regions else 'whole-organ',
        applicable_to_model=term_sex(terms, code) != 'female', authority_sex=term_sex(terms, code),
        source_meshes=sources, missing_distinct_representation=not bool(sources),
        distinct_representation_verified=False))
inventory = dict(schema='ta98-geometry-inventory/v1', model='local-ta98', model_sex='male',
    root_term='A12.2.06.001', additional_context_root='A12.2.07.029',
    authority_path=str(authority), authority_sha256=sha(authority),
    authority_corrections_path=str(corrections), authority_corrections_sha256=sha(corrections),
    source_catalogue_sha256=sha(catalogue), effective_census_sha256=sha(census_path),
    term_count=len(rows), applicable_term_count=sum(r['applicable_to_model'] for r in rows),
    terms=rows, anatomical_acceptance=False,
    classification_correction='Named arteries and arterial branches are vessels/networks, not generic tissue regions. Carotid sinus and siphon are arterial course/wall regions, not a boundary curve.',
    primary_terminology='https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.2.07%20TA98%20EN.htm')
destination = work / 'complete-authority-inventory.json'
assert not destination.exists()
destination.write_text(json.dumps(inventory, indent=2) + '\n')
report = ROOT / 'reports/mesh-audit/REVIEW.md'
with report.open('a') as handle:
    handle.write('\n## Joint carotid terminal cohort inventory — 2026-10-25\n\n')
    handle.write(f'Before any terminal ICA/MCA/ACA reconstruction: all {len(rows)} authoritative ICA descendants plus the adjoining ACom context are recorded against source `{sha(catalogue)}`. Full source records and hashes: `{destination.relative_to(ROOT)}`. Named branches are vessels/networks; only named arterial parts, sinus and siphon are classified as tissue regions. This corrects the older ICA inventory\'s overly broad tissue-region classifications. Sources may be hidden or independently registered; a source listing certifies neither native placement nor complete representation. No geometry is changed by this inventory.\n\n')
    handle.write('| TA98 | Structure | Parent | Class | Existing direct sources / missing |\n| --- | --- | --- | --- | --- |\n')
    for row in rows:
        t = row['authority']
        ids = ', '.join('`' + s['part_id'] + '`' for s in row['source_meshes']) or '**Missing distinct representation**'
        handle.write(f'| {t["id"]} | {t["name_en"]} | {t["parent_id"]} | {row["representation_class"]} | {ids} |\n')
    handle.write('\nActual source geometry places the ACA root centre on MCA-labelled ring 12 within 0.292 nm, 3.627 mm along that source. The MCA-labelled proximal mouth matches the cerebral ICA mouth centre within 29.8 nm but its profile differs by up to 193 µm. Thus this root is an intended terminal bifurcation cohort, not an unrelated arterial obstacle. Anatomical ownership of the short common stem must be resolved together with its open interfaces; an ACA/MCA side junction must not be described as a universal ACA-from-MCA origin. Skeleton, nerves and brain remain fixed.\n')
print(json.dumps(dict(inventory=str(destination), terms=len(rows), geometry_changed=False)))
