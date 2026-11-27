"""Deterministic local donor identity inventory; never geometry acceptance.

Run under geometry_tools.runtime. Exact labels are retrieval evidence, and
source hierarchy may disambiguate them. No approximate or nearest-name match.
"""
from __future__ import annotations
import argparse, collections, hashlib, json, re, sqlite3, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.ta98_authority import corrected_terms, term_sex

DEFAULT_CENSUS = ROOT / '.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3/proposal-v2/census-installed.json'
DEFAULT_RESOLVED = DEFAULT_CENSUS.with_name('resolved-installed.json')
DEFAULT_SOURCE = Path('/Users/christiantzurcanu/Documents/dev/fetch/viewer/data/male/catalog.json')
DEFAULT_OUTPUT = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1/ta98-inventory'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def normal(label):
    label = label.lower().strip().replace('’', "'").replace('-', '')
    label = re.sub(r'\.(l|r)$', lambda m: ' left' if m[1] == 'l' else ' right', label)
    return re.sub(r'\s+', ' ', re.sub(r'[^a-z0-9]+', ' ', label)).strip()

def side(label):
    sides = sorted(set(re.findall(r'\b(?:left|right)\b', normal(label))))
    return sides[0] if len(sides) == 1 else ('mixed' if sides else None)

def unsided(label):
    return re.sub(r'\s+', ' ', re.sub(r'\b(?:left|right)\b', '', normal(label))).strip()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--source-directory', type=Path, help='Hash-bound corrected local import directory, e.g. source-v2; original catalog remains --source')
    parser.add_argument('--identity-rules', type=Path, help='Explicit reviewed source-name/hierarchy and authority-parent rules; never generic parent mapping')
    parser.add_argument('--census', type=Path, default=DEFAULT_CENSUS)
    parser.add_argument('--resolved', type=Path, default=DEFAULT_RESOLVED)
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    authority = ROOT / '.local-models/ta98-authority.sqlite'
    with sqlite3.connect('file:' + str(authority) + '?mode=ro', uri=True) as con:
        con.row_factory = sqlite3.Row
        terms = corrected_terms({r['id']: dict(r) for r in con.execute('select * from ta98')})
        synonyms = [dict(r) for r in con.execute('select id,synonym,lang from synonyms')]
    census = json.loads(args.census.read_bytes())
    resolved = json.loads(args.resolved.read_bytes())
    installed = ROOT / census['source_catalogue']
    assert sha(installed) == census['source_catalogue_sha256'] == resolved['sourceSha256'], 'Installed revision changed'
    assert sha(authority) == census['authority_cache_sha256'], 'Authority revision changed'
    children = collections.defaultdict(list)
    for code, row in terms.items():
        children[row.get('parent_id')].append(code)
    def closure(code):
        found, stack = {code}, [code]
        while stack:
            for child in children[stack.pop()]:
                if child not in found:
                    found.add(child); stack.append(child)
        return found
    excluded_heart = closure('A12.1.00.001')
    excluded_lymphatic = set().union(*(closure(c) for c in census['authority_roots']['excluded_lymphatic']))
    scope = closure('A12.0.00.000') - excluded_heart - excluded_lymphatic - closure('A13.0.00.000')
    prior = {r['id']: r for r in census['terms']}
    assert set(prior) == scope and len(scope) == 1178
    male = {c for c in scope if term_sex(terms, c) != 'female'}
    assert len(male) == 1149
    ancestors = {}
    for code in scope:
        seen, current = set(), terms[code].get('parent_id')
        while current in terms and current not in seen:
            seen.add(current); current = terms[current].get('parent_id')
        ancestors[code] = seen
    lexicon = collections.defaultdict(list)
    def add(code, label, evidence):
        if code in male and label:
            record = (code, evidence, label)
            if record not in lexicon[normal(label)]: lexicon[normal(label)].append(record)
    for code in sorted(male): add(code, terms[code]['name_en'], 'authority-English-label')
    for row in synonyms:
        if row['lang'].lower() in {'en', 'eng', 'english'}: add(row['id'], row['synonym'], 'authority-English-synonym')
    # TA98 often gives a short branch name below its parent vessel. An exact
    # grammatical expansion makes that explicit without a proximity guess.
    for code in sorted(male):
        label = terms[code]['name_en']
        for ancestor in sorted(ancestors[code] & male):
            parent_label = terms[ancestor]['name_en']
            if re.search(r'\b(?:artery|arteries|vein|veins|trunk|plexus|sinus)\b', parent_label):
                add(code, label + ' of ' + parent_label, 'authority-label-plus-vascular-ancestor')
    # Exact current identities supply explicit, traceable linguistic aliases.
    # They do not transfer the installed donor's shape acceptance to this donor.
    for binding in resolved['conceptBindings']:
        if binding.get('kind') == 'exact': add(binding['term'], binding.get('name'), 'installed-exact-concept-label')
    catalog = json.loads(args.source.read_bytes())
    root = next(n for n in catalog['sceneGraph'] if n['objectId'] == 'human_20_male_cardiovascular_system-root')
    records, all_nodes, heart_records = [], {}, []
    def walk(node, parents, group):
        if node.get('displayName') in {'Arterial system', 'Venous system', 'Heart'}: group = node['displayName']
        key = node['objectId']; all_nodes[key] = dict(node=node, parents=parents, group=group)
        if node.get('geometryId'):
            row = dict(source_id=key, source_geometry=node['geometryId'], source_name=node.get('displayName'), source_group=group,
                source_parent_id=parents[-1]['objectId'] if parents else None,
                source_parent_name=parents[-1].get('displayName') if parents else None,
                source_hierarchy=[dict(id=p['objectId'], name=p.get('displayName')) for p in parents], source_side=side(node.get('displayName', '')))
            (heart_records if group == 'Heart' else records).append(row)
        for child in node.get('objects', []): walk(child, parents + [node], group)
    walk(root, [], None)
    assert len(records) == 1699 and len(heart_records) == 46
    assert collections.Counter(r['source_group'] for r in records) == {'Arterial system': 919, 'Venous system': 780}
    audited_heart_records = list(heart_records)
    source_import_inputs = []
    source_scope_policy = None
    if args.source_directory:
        directory = args.source_directory.resolve()
        source_import_inputs = [directory / name for name in ['manifest.json', 'instances.json', 'source-scope-policy.json', 'original-cardiovascular-hierarchy.json']]
        source_manifest = json.loads((directory / 'manifest.json').read_bytes())
        imported_instances = json.loads((directory / 'instances.json').read_bytes())['instances']
        source_scope_policy = json.loads((directory / 'source-scope-policy.json').read_bytes())
        assert source_manifest['counts']['instances'] == len(imported_instances) == 1724
        assert source_manifest['counts']['uniqueGeometryAssets'] == 1165
        assert source_manifest['sourceFiles'][str(args.source.resolve())]['sha256'] == sha(args.source)
        selected = {r['id']: r for r in imported_instances}
        pool = {r['source_id']: r for r in records + heart_records}
        assert set(selected).issubset(pool)
        policy_included = {r['sourceId'] for r in source_scope_policy['leaves'] if r['include']}
        assert len(policy_included) == 25
        assert policy_included == set(selected) & {r['source_id'] for r in heart_records}
        records = [pool[key] for key in selected]
        heart_records = [r for r in heart_records if r['source_id'] not in selected]
        assert len(heart_records) == 21
        for row in records:
            imported = selected[row['source_id']]
            row.update(source_original_group=row['source_group'], source_group=imported['bloodSystem'], source_vascular_role=imported['vascularRole'], source_geometry_canonical_id=imported['geometryId'], exact_geometry_ta98_binding_asserted=False)
        for key, item in all_nodes.items():
            names = {p.get('displayName') for p in item['parents']} | {item['node'].get('displayName')}
            if key in selected: item['group'] = selected[key]['bloodSystem']
            elif 'Coronary arteries' in names: item['group'] = 'Arterial system'
            elif 'Coronary veins' in names: item['group'] = 'Venous system'
    def candidates(name, group, context=False):
        if context:
            name = re.sub(r'^(?:branches|tributaries) of ', '', name, flags=re.I)
            name = re.sub(r' and branches$', '', name, flags=re.I)
        exact = list(lexicon.get(normal(name), []))
        if not exact and side(name) in {'left', 'right'}:
            exact = [v for v in lexicon.get(unsided(name), []) if side(terms[v[0]]['name_en']) is None]
        prefix = {'Arterial system': 'A12.2.', 'Venous system': 'A12.3.', 'Arterial remnant companions': 'A12.2.'}.get(group)
        source_side = side(name)
        def compatible(code):
            if source_side not in {'left', 'right'}: return True
            # Only a nearest explicitly sided authority ancestor can qualify a
            # short child label; this preserves right/left pulmonary closures.
            current, seen = code, set()
            while current in terms and current not in seen:
                seen.add(current)
                qualifier = side(terms[current]['name_en'])
                if qualifier in {'left', 'right'}: return qualifier == source_side
                current = terms[current].get('parent_id')
            return True
        return sorted(set(v for v in exact if (not prefix or v[0].startswith(prefix)) and compatible(v[0])))
    # Resolve only explicit linguistic matches; source ancestor identities select
    # among duplicate labels only when one TA98 candidate has stronger ancestry.
    node_matches = {}
    for key, item in all_nodes.items():
        hits = candidates(item['node'].get('displayName', ''), item['group'], context=True)
        codes = {v[0] for v in hits}
        if len(codes) == 1: node_matches[key] = next(iter(codes))
    mapping = collections.defaultdict(list)
    reviewed_rules = json.loads(args.identity_rules.read_bytes()) if args.identity_rules else None
    explicit_rules = {r['source_id']:r for r in reviewed_rules['in_denominator_rules']} if reviewed_rules else {}
    for row in records:
        hits = candidates(row['source_name'], row['source_group'])
        explicit = explicit_rules.get(row['source_id'])
        if explicit:
            code = explicit['authority']['id']
            assert code in male and row['source_name'] == explicit['source_name'] and row['source_hierarchy'] == explicit['source_hierarchy']
            assert terms[code]['name_en'] == explicit['authority']['name_en'] and terms[code]['parent_id'] == explicit['authority']['parent_id']
            assert terms[terms[code]['parent_id']]['name_en'] == explicit['authority']['parent_name']
            assert not hits, 'Reviewed rules may only recover previously unmatched labels'
            hits = [(code, 'reviewed-context-rule:'+explicit['rule_id'], row['source_name'])]
        codes = sorted({v[0] for v in hits})
        path_codes = {node_matches[p['id']] for p in row['source_hierarchy'] if p['id'] in node_matches}
        scores = {c: len(ancestors[c] & path_codes) for c in codes}
        chosen = codes[0] if len(codes) == 1 else None
        hierarchy_used = False
        if len(codes) > 1:
            best = max(scores.values())
            winners = [c for c in codes if scores[c] == best]
            if best > 0 and len(winners) == 1: chosen = winners[0]; hierarchy_used = True
        row.update(status='candidate-identity' if chosen else ('ambiguous-identity' if codes else 'unmapped'),
            ta98_id=chosen, ta98_name=terms[chosen]['name_en'] if chosen else None,
            confidence='exact-linguistic-and-source-hierarchy-candidate' if hierarchy_used else ('exact-linguistic-candidate' if chosen else None),
            evidence=[dict(ta98_id=c, matched_label=label, method=method, exact_parent_matches=scores[c]) for c, method, label in hits],
            candidate_ta98_ids=codes, source_hierarchy_disambiguated=hierarchy_used,
            geometry_acceptance=False, finite_wall_lumen_patency_verified=False)
        if chosen: mapping[chosen].append(row['source_id'])
        if explicit:
            assert chosen == explicit['authority']['id']
            row.update(confidence='reviewed-source-name-and-authority-context-candidate', reviewed_context_rule=explicit)
    instances_per_geometry = collections.Counter(r['source_geometry'] for r in records)
    for row in records: row['geometry_shared_by_instance_count'] = instances_per_geometry[row['source_geometry']]
    rows = []
    for code in sorted(scope):
        r, old = terms[code], prior[code]
        source_ids = sorted(mapping[code])
        rows.append(dict(id=code, name_en=r['name_en'], name_la=r['name_la'], parent_id=r['parent_id'], parent_name=terms.get(r['parent_id'], {}).get('name_en'),
            source_parent_id=r.get('source_parent_id', r['parent_id']), male_applicable=code in male, sex_qualification=term_sex(terms, code),
            planned_representation_class='non-patent-arterial-remnant-companion' if args.source_directory and code == 'A12.2.01.202' else old['planned_representation_class'], authority_immaterial=bool(r['immaterial']), variant=bool(r['variant']),
            existing_installed_status=old['representation_status'], existing_installed_direct_part_ids=old['resolved_direct_part_ids'],
            existing_installed_suppressed_exact_part_ids=old['suppressed_exact_source_part_ids'],
            donor_candidate_source_ids=source_ids, donor_status='candidate-identity' if source_ids else ('not-male-applicable' if code not in male else 'missing-donor-identity'),
            anatomical_acceptance=False, wall_lumen_patency_and_interfaces_verified=False))
    covered = {code for code, source_ids in mapping.items() if source_ids}
    existing = {c for c in male if prior[c]['resolved_direct_part_ids']}
    summary = dict(male_applicable_terms=len(male), all_sex_terms=len(scope), donor_blood_mesh_instances=len(records),
        donor_blood_unique_geometry_ids=len(instances_per_geometry), excluded_heart_mesh_instances=len(heart_records),
        source_instance_groups=dict(sorted(collections.Counter(r['source_group'] for r in records).items())),
        source_instance_status=dict(sorted(collections.Counter(r['status'] for r in records).items())),
        donor_candidate_distinct_ta98_terms=len(covered), donor_missing_ta98_terms=len(male-covered),
        donor_candidate_term_percent=round(100*len(covered)/len(male), 3),
        existing_installed_direct_reference_terms=len(existing), donor_candidate_terms_additional_to_installed=len(covered-existing),
        potential_union_identity_reference_terms=len(covered|existing), potential_union_percent=round(100*len(covered|existing)/len(male),3),
        accepted_anatomy_terms=0)
    if args.source_directory:
        summary.update(audited_source_heart_mesh_instances=46, included_source_heart_vascular_instances=24,
            included_source_heart_non_patent_remnant_companions=1, included_venous_valve_companions=1,
            excluded_cardiac_or_ostial_source_instances=21)
    output = dict(schema='biodigital-male-blood-ta98-candidate-inventory/v1',
        inputs={str(p): sha(p) for p in [args.source, args.census, args.resolved, authority, installed] + source_import_inputs},
        scope_policy=census['scope_policy'], authority_roots=census['authority_roots'], summary=summary,
        identity_policy='Exact normalized English authority/synonym, authority label of vascular ancestor, or installed exact concept labels; source arterial/venous hierarchy and explicit sided authority ancestry constrain scope; ancestor identity resolves duplicate labels only with a unique positive score. Source hierarchy grammatical wrappers branches/tributaries of and and branches may be removed for ancestor matching. No fuzzy match. Display-name laterality retained independently of shared geometry IDs.',
        acceptance_policy='Candidate labels establish no geometry, wall, layer, lumen, junction, organ-interface, fit, completeness or redistribution acceptance. Openings are not solids; microscopic and wall terms stay explicit.',
        terms=rows, source_mesh_instances=sorted(records, key=lambda r: r['source_id']),
        excluded_heart_mesh_instances=sorted(heart_records, key=lambda r:r['source_id']))
    if args.source_directory:
        output.update(schema='biodigital-male-blood-ta98-candidate-inventory/v2',
            donor_scope_policy=source_scope_policy,
            audited_original_heart_mesh_instances=sorted(audited_heart_records, key=lambda r:r['source_id']),
            scope_correction='Source Heart is an organizational grouping: include 24 great/coronary vascular leaves and one explicitly non-patent arterial remnant; exclude 21 cardiac/ostial source leaves. Original names and hierarchy remain provenance. Authoritative term scope remains unchanged.')
    if reviewed_rules:
        stage=reviewed_rules.get('inventoryStage','v3')
        assert stage in ['v3','v4']
        output.update(schema='biodigital-male-blood-ta98-candidate-inventory/'+stage, reviewed_identity_rules=reviewed_rules,
            cross_root_identity_candidates_outside_current_denominator=reviewed_rules['cross_root_out_of_denominator'],
            contextual_rule_policy='Exact reviewed source IDs/names and complete source hierarchy, bound to literal authority term and parents. No generic parent substitution, fuzzy match or runtime binding. Explicit cross-root vessels remain separate from the historical denominator.')
        output['inputs'][str(args.identity_rules)] = sha(args.identity_rules)
        summary.update(reviewed_context_rule_instances=len(explicit_rules), cross_root_candidate_instances_outside_denominator=len(reviewed_rules['cross_root_out_of_denominator']), cross_root_candidate_distinct_terms_outside_denominator=len({r['authority']['id'] for r in reviewed_rules['cross_root_out_of_denominator']}))
    args.output.mkdir(parents=True, exist_ok=True)
    destination = args.output / 'inventory.json'
    assert not destination.exists(), 'Never overwrite retained evidence'
    destination.write_text(json.dumps(output, ensure_ascii=False, separators=(',', ':'), sort_keys=True) + '\n')
    (args.output / 'summary.json').write_text(json.dumps(summary, sort_keys=True, indent=2) + '\n')
    print(json.dumps(dict(output=str(destination), **summary), sort_keys=True))

if __name__ == '__main__': main()
