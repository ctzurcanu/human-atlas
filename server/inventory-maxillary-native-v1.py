"""Source-bound complete TA98 maxillary inventory before native reconstruction."""
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


def run():
    work = ROOT / '.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1'
    work.mkdir(parents=True, exist_ok=True)
    destination = work / 'complete-authority-inventory.json'
    assert not destination.exists(), 'Retain previous source-bound evidence'
    catalogue = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    authority = ROOT / '.local-models/ta98-authority.sqlite'
    census_path = ROOT / 'reports/ta98-blood-current-20261024.json'
    census = json.loads(census_path.read_text())
    source_sha = sha(catalogue)
    assert census['source_sha256'] == source_sha
    with sqlite3.connect('file:' + str(authority) + '?mode=ro', uri=True) as conn:
        conn.row_factory = sqlite3.Row
        terms = corrected_terms({r['id']: dict(r) for r in conn.execute('SELECT * FROM ta98')})
    descendants = {'A12.2.05.053'}
    while True:
        children = {code for code, term in terms.items() if term['parent_id'] in descendants}
        if children <= descendants:
            break
        descendants |= children
    assert len(descendants) == 38
    effective = {t['id']: t for t in census['terms']}
    atlas = json.loads(catalogue.read_text())
    parts = {p['id']: p for p in atlas['parts']}
    rows = []
    for code in sorted(descendants):
        ids = effective.get(code, {}).get('direct_part_ids', [])
        sources = [dict(part_id=pid, record=parts[pid],
                        record_sha256=hashlib.sha256(json.dumps(parts[pid], sort_keys=True,
                                                               separators=(',', ':')).encode()).hexdigest(),
                        source_reference_only=True) for pid in ids]
        rows.append(dict(authority=terms[code], representation_class='whole-organ',
                         applicable_to_model=term_sex(terms, code) != 'female',
                         authority_sex=term_sex(terms, code), source_meshes=sources,
                         missing_distinct_representation=not bool(sources),
                         distinct_representation_verified=False))
    inventory = dict(schema='ta98-geometry-inventory/v1', model='local-ta98', model_sex='male',
                     root_term='A12.2.05.053', authority_path=str(authority),
                     authority_sha256=sha(authority),
                     authority_corrections_sha256=sha(ROOT / 'geometry_tools/ta98_authority.py'),
                     source_catalogue_sha256=source_sha, effective_census_sha256=sha(census_path),
                     term_count=len(rows), applicable_term_count=sum(t['applicable_to_model'] for t in rows),
                     terms=rows, anatomical_acceptance=False,
                     classification_note='All listed terms are named vessels or vessel networks. No arterial cast establishes a microscopic feature or finite wall. Clinical mandibular/pterygoid/pterygopalatine segments are not extra TA98 descendants in this authority.',
                     primary_terminology='https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.2.05%20TA98%20EN.htm')
    destination.write_text(json.dumps(inventory, indent=2) + '\n')
    with (ROOT / 'reports/mesh-audit/REVIEW.md').open('a') as handle:
        handle.write('\n## Native male maxillary artery — complete inventory, 2026-10-25\n\n')
        handle.write(f'All 38 authoritative descendants, parents and direct sources are recorded before construction against male master `{source_sha}`. Full source records/hashes: `{destination.relative_to(ROOT)}`. A source reference can be hidden or independently registered and is not a verified recipient interface. The native maxillary source has multiple authored courses; their identity and communication require actual source/recipient review. No geometry changes or accepted new terms follow from this inventory.\n\n')
        handle.write('| TA98 | Structure | Parent | Class | Existing direct sources / missing |\n| --- | --- | --- | --- | --- |\n')
        for row in rows:
            term = row['authority']
            ids = ', '.join('`' + s['part_id'] + '`' for s in row['source_meshes']) or '**Missing distinct representation**'
            handle.write(f'| {term["id"]} | {term["name_en"]} | {term["parent_id"]} | whole-organ | {ids} |\n')
        handle.write('\nNamed arteries and dental/peridental networks require actual 3D courses. An opening is an interface, not an added solid. Preserve the literal source controls, normals and UVs when recovering components; do not classify a terminal branch by nearest tooth alone. Skeleton, teeth, existing recipients, basilar artery and female geometry remain fixed.\n')
    print(json.dumps(dict(inventory=str(destination), terms=len(rows), source_sha256=source_sha,
                          missing=sum(t['missing_distinct_representation'] for t in rows), geometry_changed=False)))


if __name__ == '__main__':
    run()
