"""Read-only anatomical baseline; writes an inventory, never geometry or catalogues."""
from pathlib import Path
import hashlib
import json
import sqlite3
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.ta98_authority import corrected_terms

AUTHORITY = ROOT / '.local-models/ta98-authority.sqlite'
METADATA = ROOT / 'app/data/ta98-metadata.json'
SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas'
REPORT = ROOT / 'reports/ta98-tongue-inventory.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def classification(code):
    if code in ('A04.1.05.001', 'A05.1.04.001'):
        return 'whole-organ'
    suffix = code.rsplit('.', 1)[-1]
    if suffix in ('004', '005', '006', '007'):
        return 'surface-region'
    if suffix in ('009', '018', '019'):
        return 'boundary-curve'
    if suffix == '010':
        return 'point-landmark'
    if suffix == '011':
        return 'wall-layer'
    if suffix == '020':
        return 'opening'
    if suffix in ('013', '014', '015', '016', '017', '023'):
        return 'microscopic-structure'
    return 'tissue-region'


def representation(code):
    if code == 'A04.1.05.001':
        return 'Assembly of the actual muscle members; no additional overlapping solid.'
    if code == 'A05.1.04.001':
        return 'Complete organ assembly with correctly owned external surface and internal tissue.'
    if code == 'A05.1.04.020':
        return 'Blind surface pit/orifice with a selectable site; not a separate solid or an invented through-hole.'
    if code == 'A05.1.04.021':
        return 'Variant-only developmental remnant; keep absent from the normal adult assembly unless specifically represented as a variant.'
    kind = classification(code)
    return {
        'surface-region': 'Selectable subset of the correct existing surface; not an extra whole-organ shell.',
        'boundary-curve': 'Curve anchored to the corresponding margin or sulcus on the actual organ.',
        'point-landmark': 'Point anchored to the actual apex.',
        'wall-layer': 'Finite mucosal coating with a correct inner boundary; outer tongue geometry alone is insufficient.',
        'microscopic-structure': 'Scale-appropriate relief or grouped small structures; distinguish type and distribution without one mesh per repeated element.',
        'tissue-region': 'Owned finite tissue or a qualified regional structure; preserve adjoining interfaces and attachments.',
    }[kind]


def main():
    started = time.monotonic()
    with sqlite3.connect('file:' + str(AUTHORITY) + '?mode=ro', uri=True) as conn:
        conn.row_factory = sqlite3.Row
        authority = corrected_terms({row['id']: dict(row) for row in conn.execute('SELECT * FROM ta98')})
    # The muscle, tonsil and palatoglossus parents are canonical cross-references,
    # not descendants of the tongue root. Preserve their actual parents.
    codes = sorted(code for code in authority if code.startswith('A05.1.04.') or code == 'A04.1.05.001')
    assert len(codes) == 36
    metadata = json.loads(METADATA.read_text())
    by_term = {sex: {} for sex in ('male', 'female')}
    catalogues, bank_hashes = {}, {}
    for sex in ('male', 'female'):
        path = SOURCE / sex / 'atlas.json'
        atlas = json.loads(path.read_text())
        parts = {part['id']: part for part in atlas['parts']}
        hidden = set(atlas.get('defaultHidden', []))
        matched = {}
        for concept in atlas['concepts']:
            mapped = metadata['byConceptMatch'].get(concept['id'], {})
            term = concept.get('ta98Term') or mapped.get('term') or metadata['byConcept'].get(concept['id'], {}).get('ta98')
            if term not in codes:
                continue
            for pid in concept['elements']:
                if pid not in parts:
                    continue
                matched.setdefault((term, pid), []).append(dict(concept_id=concept['id'], kind=concept.get('ta98Kind') or mapped.get('kind'), method=concept.get('ta98Method') or mapped.get('method')))
        for part in atlas['parts']:
            if part.get('ta98Term') in codes:
                matched.setdefault((part['ta98Term'], part['id']), []).append(dict(part_record_term=True, kind=part.get('ta98Kind')))
        for (term, pid), mappings in sorted(matched.items()):
            part = parts[pid]
            bank = atlas['chunks'][part['chunk']]
            url = bank['url']
            bank_path = ROOT / '.local-models' / url.removeprefix('/local-models/')
            if url not in bank_hashes:
                assert bank_path.stat().st_size == bank['bytes']
                bank_hashes[url] = dict(path=str(bank_path), bytes=bank_path.stat().st_size, sha256=digest(bank_path))
            by_term[sex].setdefault(term, []).append(dict(part_id=pid, name=part['name'], concept_id=part['conceptId'], mappings=mappings,
                hidden=pid in hidden, suppressed=bool(part.get('suppressed')), triangles=part['indexCount'] // 3,
                bounds=part['bounds'], bank_url=url, bank_sha256=bank_hashes[url]['sha256'],
                review_status=part.get('reviewStatus'), coverage_limitation=part.get('coverageLimitation'), record=part))
        body = ROOT / f'.local-models/organ-models/spinal-cord/{sex}-body.json'
        catalogues[sex] = dict(path=str(path), sha256=digest(path), ordinary_body_catalogue=str(body), ordinary_body_sha256=digest(body) if body.exists() else None)
    rows = []
    for code in codes:
        sources = {sex: by_term[sex].get(code, []) for sex in ('male', 'female')}
        rows.append(dict(authority=authority[code], representation_class=classification(code), proposed_representation=representation(code),
            required_in_typical_adult=code != 'A05.1.04.021', canonical_parent_preserved=True,
            source_meshes=sources, no_source_geometry={sex: not sources[sex] for sex in sources},
            distinct_representation_verified=False, missing_or_unverified={sex: True for sex in sources},
            acceptance='No anatomical or attachment acceptance established by this inventory; a source identity is not proof of representation.'))
    report = dict(schema='ta98-tongue-inventory/v1', authority=dict(path=str(AUTHORITY), sha256=digest(AUTHORITY),
        source=metadata['source'], corrections_path=str(ROOT / 'geometry_tools/ta98_authority.py'), corrections_sha256=digest(ROOT / 'geometry_tools/ta98_authority.py')),
        metadata=dict(path=str(METADATA), sha256=digest(METADATA)), source_catalogues=catalogues, source_banks=bank_hashes,
        term_count=len(rows), normal_adult_applicable_terms=35, terms=rows,
        assembly=dict(whole='A05.1.04.001', muscle_group='A04.1.05.001',
            cross_references={'A05.1.04.022': 'A13.2.02.001', 'A05.1.04.023': 'A05.1.04.022', 'A05.1.04.110': 'A04.1.06.001'},
            fixed_skeleton=True, interfaces=['Genioglossus to mandibular genial spines and tongue interior', 'Hyoglossus/chondroglossus/ceratoglossus to the appropriate hyoid sites and tongue', 'Styloglossus to the styloid region and lateral tongue', 'Palatoglossus to soft palate and lateral tongue', 'Intrinsic muscle territories, lingual septum and aponeurosis', 'Mucosa, surface regions, papillae and posterior lymphoid tissue', 'Frenulum and floor of mouth']),
        source_only_term_counts={sex: sum(bool(row['source_meshes'][sex]) for row in rows) for sex in by_term},
        no_source_term_counts={sex: sum(row['no_source_geometry'][sex] for row in rows) for sex in by_term},
        anatomically_accepted_terms=0, geometry_changed=False, catalogue_changed=False,
        body_presentation='Keep a small whole-body assembly. Detailed representations can use a separate linked model; its mesh total must not enter the body hierarchy totals.',
        elapsed_seconds=round(time.monotonic() - started, 4))
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    lines = ['# TA98 tongue inventory', '', 'Baseline inventory only. No geometry or catalogue is changed, and no anatomical acceptance is claimed.', '',
        'The authority contains **36 relevant terms**, including the separate muscle group, lingual tonsil and palatoglossus cross-references. Their canonical parents are preserved. The thyroglossal duct is a variant and is not required in the ordinary adult assembly.', '',
        'Both models have source evidence for only **5 terms**; **31 terms per model have no matching source geometry**. Source evidence includes container membership and hidden, misregistered sources, so it is not a completion count.', '',
        '| TA98 | Term | Parent | Representation | Male source meshes | Female source meshes |', '|---|---|---|---|---|---|']
    for row in rows:
        a = row['authority']
        cells = []
        for sex in ('male', 'female'):
            ss = row['source_meshes'][sex]
            cells.append('; '.join('`' + s['part_id'] + '`' + (' (hidden)' if s['hidden'] else '') + (' (suppressed)' if s['suppressed'] else '') for s in ss) or '**missing**')
        lines.append(f"| {a['id']} | {a['name_en']} | {a['parent_id']} | {row['representation_class']} | {' | '.join(cells)} |")
    lines += ['', '## Current source limitations', '',
        '- Male: the active whole tongue has 343 vertices / 682 triangles. Paired genioglossus and hyoglossus are present. The two native styloglossus sources are hidden and explicitly await alignment/overlap review; do not simply unhide them.',
        '- Female: the native whole tongue has 205 vertices / 304 triangles. Genioglossus and paired hyoglossus/styloglossus are present. The hidden HRA dorsal source is large and misregistered; it does not establish an installed dorsum.',
        '- The existing oral-muscle plan explicitly leaves intrinsic muscles, palatoglossus, chondroglossus/ceratoglossus boundaries, attachments and fibre architecture unverified or missing.',
        '- Surface regions use owned surface subsets, sulci/margins use curves, the apex uses a point, and the foramen caecum is a blind surface pit. None should become an unrelated solid.',
        '- Mucosa needs finite tissue ownership. Septum/aponeurosis are connective tissue. Papilla types need correct scale and distribution, preferably grouped relief rather than a mesh for every repeated element.',
        '- The male assembly is completed first with the skeleton fixed, then transferred and fitted to female mandibular, hyoid, styloid and palatal attachments. Existing source muscle envelopes alone do not establish those attachments.',
        '- Keep the body assembly economical. A separate detailed organ model can own additional meshes and hierarchies, with its count excluded from the body model.', '',
        'Full source records, catalogue/binary hashes and individual missing/unverified flags are in [reports/ta98-tongue-inventory.json](reports/ta98-tongue-inventory.json). Reproduce with the bounded preparation script:', '',
        '```sh', 'python3 -B -m geometry_tools.runtime --seconds 20 -- python3 -B server/prepare-tongue-inventory.py', '```', '']
    (ROOT / 'TA98-TONGUE.md').write_text('\n'.join(lines))
    print(json.dumps(dict(report=str(REPORT), terms=36, source_only_term_counts=report['source_only_term_counts'], no_source_term_counts=report['no_source_term_counts'], geometry_changed=False, elapsed_seconds=report['elapsed_seconds'])))


if __name__ == '__main__':
    main()
