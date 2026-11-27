"""Expose already-visible named native sources; never change/activate geometry."""
from pathlib import Path
import copy, json, sqlite3, sys, time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_bytes, atomic_json, digest_file
from geometry_tools.ta98_authority import corrected_terms

SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
SHA = 'b26860231a10b4b4382802751d6a8c15a7491e8baa83d629d73b93f5abb61707'
OUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/native-pulmonary-controls-v1'
PAIRS = [
    ('A12.2.01.102', 'Superior lobar artery of right lung'),
    ('A12.2.01.103', 'Apical segmental artery of right lung'),
    ('A12.2.01.110', 'Middle lobar artery of right lung'),
    ('A12.2.01.113', 'Inferior lobar artery of right lung'),
    ('A12.2.01.114', 'Superior segmental artery of right lung'),
    ('A12.2.01.116', 'Anterior basal segmental artery of right lung'),
    ('A12.2.01.117', 'Lateral basal segmental artery of right lung'),
    ('A12.2.01.119', 'Posterior basal segmental artery of right lung'),
    ('A12.2.01.204', 'Apical segmental artery of left lung'),
    ('A12.2.01.205', 'Anterior segmental artery of left lung'),
    ('A12.2.01.208', 'Posterior segmental artery of left lung'),
    ('A12.2.01.215', 'Superior segmental artery of left lung'),
    ('A12.2.01.217', 'Anterior basal segmental artery of left lung'),
    ('A12.2.01.218', 'Lateral basal segmental artery of left lung'),
    ('A12.2.01.219', 'Medial basal segmental artery of left lung'),
    ('A12.2.01.220', 'Posterior basal segmental artery of left lung'),
]

def main():
    start = time.monotonic(); assert digest_file(SOURCE) == SHA
    before = json.loads(SOURCE.read_bytes()); after = copy.deepcopy(before)
    parts = {p['id']: p for p in after['parts']}; concepts = {c['id']: c for c in after['concepts']}
    prior = ROOT / '.local-models/geometry-workbench/circulatory-completion/lateral-circumflex-femoral-native-v1/source-ownership-proposal-v4'
    census = json.loads((prior / 'census-installed.json').read_bytes()); rows = {r['id']: r for r in census['terms']}
    resolved = json.loads((prior / 'resolved-installed.json').read_bytes())
    visible = {r['id']: r['available'] for r in resolved['effectivePartVisibility']}
    gap = json.loads((ROOT / '.local-models/geometry-workbench/circulatory-completion/next-ready-batch-v1.json').read_bytes())
    gap_parts = {r['source_part_id']: r for r in gap['rows']}
    evidence = []
    for term, name in PAIRS:
        pid = 'ZA:' + name; p = parts[pid]; c = concepts[pid]; g = gap_parts[pid]
        assert visible[pid] and c['elements'] == [pid] and p['conceptId'] == pid
        assert not rows[term]['active_exact_source_part_ids']
        assert p['sourceId'] == name and p['system'] == 'arterial' and p['license'] == 'CC BY-SA 4.0'
        assert 'Pulmonary arteries' in p['groups'] and not g['literal_faces_shared_with_other_batch_parts']
        limitation = ('Existing named native pulmonary arterial outside surface, independently selectable under TA98. '
                      f"Source retains {len(g['topology']['components'])} component(s) and {g['topology']['open_edges']} open edges. "
                      'This identity/control installation changes no coordinates, vessel course, source visibility or recipient interface. '
                      'Finite wall thickness, branch-junction tissue ownership, lumen continuity and lobar perfusion remain unverified.')
        c.update(ta98Term=term, ta98Kind='exact', ta98Method='reviewed-native-source-pulmonary-identity', coverageLimitation=limitation)
        evidence.append(dict(term=term, concept_id=pid, part_id=pid, source_id=p['sourceId'],
                             geometry_sha256=g['geometry_sha256'], attributes_sha256=g['attributes_sha256'],
                             bank_path=g['bank_path'], bank_sha256=g['bank_sha256'],
                             source_groups=p['groups'], existing_topology=g['topology'], qualification=limitation))
    after['nativePulmonaryControlIntegration'] = dict(schema='native-pulmonary-source-controls/v1',
        sourceCatalogueSha256=SHA, conceptIds=[r['concept_id'] for r in evidence], partIds=[r['part_id'] for r in evidence],
        terms=[r['term'] for r in evidence], geometryUnchanged=True, sourceVisibilityUnchanged=True,
        anatomyAccepted=False, finiteWallOrPatencyClaim=False)
    assert before['parts'] == after['parts'] and before['chunks'] == after['chunks']
    assert before['materials'] == after['materials'] and before['defaultHidden'] == after['defaultHidden']
    atomic_bytes(OUT/'source-atlas.json', SOURCE.read_bytes())
    atomic_bytes(OUT/'proposed-atlas.json', (json.dumps(after,separators=(',',':'))+'\n').encode())
    authority = ROOT/'.local-models/ta98-authority.sqlite'; connection = sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True)
    connection.row_factory = sqlite3.Row
    terms = corrected_terms({r['id']:dict(r) for r in connection.execute('select * from ta98')})
    closure = {'A12.2.01.001'}
    while True:
        enlarged = closure | {k for k,v in terms.items() if v.get('parent_id') in closure}
        if enlarged == closure: break
        closure = enlarged
    assert len(closure) == 44
    selected = {r['term']:r['part_id'] for r in evidence}; inventory_rows=[]
    for code in sorted(closure):
        old = rows[code]; ids = sorted(set(old['active_exact_source_part_ids']+old['suppressed_exact_source_part_ids']+([selected[code]] if code in selected else [])))
        inventory_rows.append(dict(authority=terms[code], representation_class='luminal-content' if code.endswith('.002') else 'tissue-region' if code in ['A12.2.01.003','A12.2.01.202'] else 'whole-organ',
            distinct_representation_verified=False, representation_provisional=True,
            missing_distinct_representation=not bool(old['active_exact_source_part_ids']),
            source_meshes=[dict(part_id=pid,source_catalogue=str(OUT/'source-atlas.json')) for pid in ids],
            planned_representation=dict(part_ids=[selected[code]] if code in selected else [], scope='Native outside-surface identity control only; no geometry, wall, patency or anatomical completion claim.')))
    inventory = dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',root_term='A12.2.01.001',
        term_count=44,applicable_term_count=44,terms=inventory_rows,source_catalogue=str(OUT/'source-atlas.json'),source_catalogue_sha256=SHA,
        authority_path=str(authority),authority_sha256=digest_file(authority),authority_corrections_path=str(ROOT/'geometry_tools/ta98_authority.py'),
        authority_corrections_sha256=digest_file(ROOT/'geometry_tools/ta98_authority.py'),checked_in_report=str(ROOT/'reports/mesh-audit/REVIEW.md'))
    atomic_json(OUT/'source-bound-inventory.json',inventory)
    proof = dict(schema='native-pulmonary-identity-control-proof/v1',source_catalogue_sha256=SHA,proposed_catalogue_sha256=digest_file(OUT/'proposed-atlas.json'),
        existing_part_records_literal=True,existing_banks_literal=True,materials_literal=True,visibility_literal=True,
        distinct_new_direct_reference_terms=16,new_mesh_records=0,new_triangles=0,new_geometry_bytes=0,anatomical_acceptance=False,installed=False,
        rows=evidence,elapsed_seconds=time.monotonic()-start)
    atomic_json(OUT/'identity-proof.json',proof)
    spec = dict(model='local-ta98',proposed_catalogue=str(OUT/'proposed-atlas.json'),expected_catalogue_sha256=SHA,changed_parts=[],closed_parts=[],require_closed=False,
        inventory=str(OUT/'source-bound-inventory.json'),max_pairs=250000,max_tests=2000000,interface_checks=[],
        assembly=dict(whole_organ=dict(name='Unchanged native pulmonary arterial source controls',authority_roots=['A12.2.01.001'],part_ids=[r['part_id'] for r in evidence]),
            adjoining_interfaces=[dict(first=r['part_id'],kind='retained-native-pulmonary-source-outside-interface-patency-unverified',second='ZA:'+('Right' if 'right lung' in r['part_id'] else 'Left')+' pulmonary artery') for r in evidence],
            fixed_skeleton=True,fixed_all_geometry=True,fixed_all_source_visibility=True,anatomical_acceptance=False),
        evidence=dict(identity_proof=str(OUT/'identity-proof.json')),
        remaining_limits=['Sixteen existing named native outside-surface controls; no geometry changes or new triangles.',
                          'Native disconnected/open source casts retain their original defects; no finite wall or communicating lumen claim.',
                          'Lobar recipient interfaces, perfusion and global packing remain unaccepted.'])
    atomic_json(OUT/'install-spec-draft.json',spec)
    assert digest_file(SOURCE) == SHA
    print(json.dumps(dict(output=str(OUT),proposed_sha256=proof['proposed_catalogue_sha256'],references=16,elapsed_seconds=proof['elapsed_seconds'])))

if __name__ == '__main__': main()
