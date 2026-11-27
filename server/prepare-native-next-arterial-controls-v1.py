"""Bind unchanged, already displayed native arterial sources to their identities."""
from pathlib import Path
import copy, json, sqlite3, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.audit import topology
from geometry_tools.ta98_authority import corrected_terms

B = ROOT / '.local-models/geometry-workbench/circulatory-completion/native-next-arterial-controls-v1'
SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED = '66b4a76c63043a38639d483dd744faa3789b4e468685e2746f70f738b1c1d11b'
ROWS = [
    ('A12.2.08.054', 'superficial cervical artery', 'ZA:(Superficial branch of transverse cervical artery).',
     'The named superficial transverse-cervical branch is the native superficial cervical arterial source.'),
    ('A12.2.08.058', 'dorsal scapular artery', 'O3M:Dorsal scapular artery (Deep br of transverse cervical a.).',
     'The publisher explicitly identifies the dorsal scapular artery and its deep transverse-cervical origin variant.'),
    ('A12.2.16.069', 'plantar digital arteries proper', 'ZA:Proper plantar digital arteries.',
     'The named native proper plantar digital arterial network corresponds to the collective TA98 identity.'),
    ('A12.2.08.021', 'pontine arteries', 'ZA:Medial pontine branches of basilar artery.',
     'The named medial pontine source is a partial component of the collective pontine arterial identity.'),
    ('A12.2.08.021', 'pontine arteries', 'ZA:Lateral pontine branches of basilar artery.',
     'The named lateral pontine source is a partial component of the collective pontine arterial identity.'),
    ('A12.2.09.036', 'princeps pollicis artery', 'O3M:Princeps pollicis arteries.',
     'The publisher explicitly identifies the princeps pollicis arterial source.'),
    ('A12.2.09.037', 'radialis indicis artery', 'O3M:Radialis indicis.',
     'The publisher explicitly identifies the radialis indicis arterial source.'),
    ('A12.2.09.006', 'acromial branch', 'O3M:Thoracoacromial artery Acromial br.',
     'The publisher names both the acromial branch and its thoraco-acromial parent.'),
    ('A12.2.09.009', 'deltoid branch', 'O3M:Thoracoacromial artery Deltoid br.',
     'The publisher names both the deltoid branch and its thoraco-acromial parent.'),
]
ROOTS = ['A12.2.08.053', 'A12.2.08.058', 'A12.2.16.068', 'A12.2.07.081',
         'A12.2.09.027', 'A12.2.09.005']

def main():
    began = time.monotonic()
    assert digest_file(SOURCE) == EXPECTED
    source = json.loads(SOURCE.read_bytes()); candidate = copy.deepcopy(source)
    cp = {c['id']: c for c in candidate['concepts']}
    parts = {p['id']: p for p in source['parts']}; hidden = set(source.get('defaultHidden', []))
    authority_path = ROOT / '.local-models/ta98-authority.sqlite'
    db = sqlite3.connect(authority_path); db.row_factory = sqlite3.Row
    authority = corrected_terms({r['id']: dict(r) for r in db.execute('SELECT * FROM ta98')})
    metadata = json.loads((ROOT / 'app/data/ta98-metadata.json').read_bytes())['byConceptMatch']
    reader = AtlasReader(); rows = []
    for term, name, prefix, basis in ROWS:
        assert authority[term]['name_en'] == name
        for side in ['l', 'r']:
            pid = prefix + side; p = parts[pid]; c = cp[p['conceptId']]
            assert p['system'] == 'arterial' and not p.get('suppressed') and not p.get('defaultHidden') and pid not in hidden
            assert c['elements'] == [pid] and not c.get('ta98Term')
            mesh = reader.read_part(source, p, frame='model:local-ta98:assembled-metres')
            shape = topology(mesh)
            qualification = (f"{basis} Existing native outside surface; {len(shape['components'])} source component(s), "
                             f"{shape['open_edges']} open edges. Identity/control only: coordinates, source visibility, "
                             "branch junctions and recipient interfaces are unchanged. Finite walls, lumen continuity, "
                             "complete territory coverage and anatomical fit remain unverified.")
            c.update(ta98Term=term, ta98Kind='exact', ta98Method='reviewed-native-named-arterial-identity',
                     coverageLimitation=qualification)
            rows.append(dict(term=term, concept_id=c['id'], part_id=pid, identity_basis=basis,
                             qualification=qualification, geometry_sha256=mesh.geometry_sha,
                             attributes_sha256=mesh.attributes_sha, source_id=p['sourceId'],
                             provenance=p['provenance'], existing_topology=shape,
                             bank_path=mesh.source['bank_path'], bank_sha256=mesh.source['bank_sha256']))
    candidate['nativeNamedArterialControlIntegration'] = dict(schema='native-named-arterial-controls/v1',
                terms=sorted({r['term'] for r in rows}), partIds=[r['part_id'] for r in rows],
                sourceCatalogueSha256=EXPECTED, identityOnly=True, geometryUnchanged=True,
                anatomicalAcceptance=False)
    assert candidate['parts'] == source['parts'] and candidate['chunks'] == source['chunks']
    B.mkdir(parents=True, exist_ok=False)
    source_copy = B / 'source-atlas.json'; source_copy.write_bytes(SOURCE.read_bytes())
    proposed = B / 'proposed-atlas.json'
    proposed.write_text(json.dumps(candidate, separators=(',', ':')) + '\n')
    part_by_term = {}
    for c in source['concepts']:
        code = c.get('ta98Term') or metadata.get(c['id'], {}).get('term')
        if code:
            part_by_term.setdefault(code, set()).update(c['elements'])
    inventory_paths = []; all_terms = {}
    for root in ROOTS:
        codes = {root}
        while True:
            added = {code for code, row in authority.items() if row.get('parent_id') in codes}
            if added <= codes: break
            codes |= added
        inventory_rows = []
        for code in sorted(codes):
            source_ids = sorted(pid for pid in part_by_term.get(code, set()) if pid in parts)
            planned = [r['part_id'] for r in rows if r['term'] == code]
            classification = 'whole-organ'
            if authority[code].get('immaterial'):
                classification = 'opening' if any(t in authority[code]['name_en'] for t in ['opening', 'orifice', 'ostium']) else 'boundary-curve'
            row = dict(authority=authority[code], representation_class=classification,
                       source_meshes=[dict(part_id=pid, source_catalogue=str(source_copy)) for pid in sorted(set(source_ids + planned))],
                       planned_representation=dict(part_ids=planned, scope='Existing native outside-surface identity controls only; all geometry and interfaces remain literal.'),
                       missing_distinct_representation=not bool(source_ids or planned),
                       distinct_representation_verified=False, representation_provisional=True)
            inventory_rows.append(row); all_terms[code] = row
        inventory = dict(schema='ta98-geometry-inventory/v1', root_term=root, model='local-ta98', model_sex='male',
                         source_catalogue=str(source_copy), source_catalogue_sha256=EXPECTED,
                         authority_path=str(authority_path), authority_sha256=digest_file(authority_path),
                         authority_corrections_path=str(ROOT / 'geometry_tools/ta98_authority.py'),
                         authority_corrections_sha256=digest_file(ROOT / 'geometry_tools/ta98_authority.py'),
                         checked_in_report=str(ROOT / 'reports/mesh-audit/REVIEW.md'), terms=inventory_rows,
                         term_count=len(inventory_rows), applicable_term_count=len(inventory_rows))
        path = B / 'inventories' / (root + '.json'); atomic_json(path, inventory); inventory_paths.append(str(path))
    proof = dict(schema='native-named-arterial-control-proof/v1', source_catalogue_sha256=EXPECTED,
                 proposed_catalogue_sha256=digest_file(proposed), rows=rows, inventory_paths=inventory_paths,
                 distinct_new_direct_reference_terms=len({r['term'] for r in rows}),
                 existing_part_records_literal=True, existing_banks_literal=True, visibility_literal=True,
                 new_geometry_bytes=0, new_triangles=0, new_mesh_records=0,
                 anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)
    atomic_json(B / 'identity-proof.json', proof)
    report = ROOT / 'reports/mesh-audit/REVIEW.md'
    with report.open('a') as stream:
        stream.write('\n\n### Native named arterial controls — 2026-10-26\n\n')
        stream.write('Eight collective/branch identities use eighteen existing displayed native sources. No geometry, visibility, finite walls or recipient interfaces are changed or accepted. Exact source-qualified identity does not establish a complete communicating vessel. Full relevant authority closures follow.\n\n')
        stream.write('| TA98 | Parent | Term | Existing source references | Representation |\n|---|---|---|---|---|\n')
        for code, row in sorted(all_terms.items()):
            a = row['authority']; ids = ', '.join(p['part_id'] for p in row['source_meshes']) or 'none'
            stream.write(f"| {code} | {a.get('parent_id','')} | {a['name_en']} | {ids} | {'missing' if row['missing_distinct_representation'] else 'unverified source reference'} |\n")
    print(json.dumps(dict(proposed_catalogue=str(proposed), sha256=digest_file(proposed), terms=proof['distinct_new_direct_reference_terms'], concepts=len(rows), inventory_terms=len(all_terms), elapsed_seconds=time.monotonic()-began)))

if __name__ == '__main__': main()
