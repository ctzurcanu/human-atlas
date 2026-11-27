"""Validate corrected donor identity evidence; optionally refresh one review section.

No geometry, installed catalogue, authority, source or historical inventory writes.
"""
import argparse, collections, hashlib, json, re, time
from pathlib import Path

started=time.monotonic();root=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source-directory',type=Path,required=True)
parser.add_argument('--inventory-directory',type=Path,required=True)
parser.add_argument('--update-review',action='store_true')
args=parser.parse_args();sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
path=args.inventory_directory/'inventory.json';inventory=json.loads(path.read_bytes());summary=inventory['summary']
source_path=args.source_directory/'instances.json';source=json.loads(source_path.read_bytes())['instances']
policy=json.loads((args.source_directory/'source-scope-policy.json').read_bytes())
historical=root/'.local-models/geometry-workbench/biodigital-male-circulation-v1/ta98-inventory/inventory.json'
historical_sha=sha(historical);assert historical_sha=='066e1560f8720832adaa286e29556b306c47ce5509ec9081bc9650bd4bcba976'
records=inventory['source_mesh_instances'];terms=inventory['terms']
assert len(records)==len(source)==1724 and len(terms)==1178
assert {r['source_id'] for r in records}=={r['id'] for r in source}
assert len({r['source_geometry'] for r in records})==1165
assert collections.Counter(r['source_group'] for r in records)=={'Arterial system':933,'Venous system':790,'Arterial remnant companions':1}
assert sum(r['male_applicable'] for r in terms)==1149
assert all(not r['anatomical_acceptance'] and not r['wall_lumen_patency_and_interfaces_verified'] for r in terms)
assert all(not r['geometry_acceptance'] and not r['finite_wall_lumen_patency_verified'] and not r['exact_geometry_ta98_binding_asserted'] for r in records)
covered={r['ta98_id'] for r in records if r['status']=='candidate-identity'}
stage=inventory['schema'].rsplit('/',1)[-1]
assert len(covered)==summary['donor_candidate_distinct_ta98_terms']=={'v2':482,'v3':492,'v4':550}[stage]
assert {r['id'] for r in terms if r['donor_candidate_source_ids']}==covered
included={r['sourceId'] for r in policy['leaves'] if r['include']};excluded={r['sourceId'] for r in policy['leaves'] if not r['include']}
assert len(included)==25 and len(excluded)==21 and included.issubset({r['source_id'] for r in records}) and excluded.isdisjoint({r['source_id'] for r in records})
assert len(inventory['excluded_heart_mesh_instances'])==21 and len(inventory['audited_original_heart_mesh_instances'])==46
assert next(r for r in terms if r['id']=='A12.2.01.202')['planned_representation_class']=='non-patent-arterial-remnant-companion'
for filename,digest in inventory['inputs'].items():assert sha(filename)==digest,filename
if stage in ['v3','v4']:
 previous_stage='v2' if stage=='v3' else 'v3';previous_sha='ed03708e161a9cdc20a914b660243e2ed322933955f6c95a4857968b0dcf0490' if stage=='v3' else '6c9160f030c698da9c842340e698e1f1ccf8c0d3f4f22ec5cb862fc37cf94a31'
 prior_path=root/('.local-models/geometry-workbench/biodigital-male-circulation-v1/ta98-inventory-'+previous_stage+'/inventory.json');assert sha(prior_path)==previous_sha
 prior=json.loads(prior_path.read_bytes());before={r['source_id']:r for r in prior['source_mesh_instances']};rules=inventory['reviewed_identity_rules'];assert len(rules['in_denominator_rules'])==(17 if stage=='v3' else 298) and len(rules['cross_root_out_of_denominator'])==(4 if stage=='v3' else 8)
 changed={r['source_id'] for r in records if r['ta98_id']!=before[r['source_id']]['ta98_id']};assert changed=={r['source_id'] for r in rules['in_denominator_rules'] if before[r['source_id']]['status']=='unmapped'}
 for rule in rules['in_denominator_rules']:
  old=before[rule['source_id']];assert rule['source_name']==old['source_name'] and rule['source_hierarchy']==old['source_hierarchy'] and rule['candidate_identity_only'] and not rule['runtime_binding'] and not rule['anatomical_acceptance']
  if old['status']=='candidate-identity':assert old['ta98_id']==rule['authority']['id']
  else:assert old['status']=='unmapped'
 for rule in rules['cross_root_out_of_denominator']:assert rule['authority']['id'] not in {r['id'] for r in terms} and rule['out_of_current_denominator'] and rule['candidate_identity_only'] and not rule['runtime_binding']
outside_section_unchanged=None
if args.update_review:
 assert stage=='v2', 'Focused v3/v4 evidence does not automatically replace the published v2 review section'
 review=root/'reports/mesh-audit/REVIEW.md';text=review.read_text();heading='### BioDigital male whole-circulation donor inventory — 2026-10-26';start=text.index(heading);following=re.search(r'^#{1,3} ',text[start+len(heading):],re.M);end=start+len(heading)+following.start() if following else len(text)
 intro=f'''{heading}

Corrected source-v2 scope retains 1,724 named instances (933 arterial, 790 venous and one explicitly non-patent arterial-remnant companion), using 1,165 shared source assets. All 46 leaves under the source's organizational Heart group are audited: 24 great/coronary vascular leaves and ligamentum arteriosum are included; 21 cardiac or ostial reference leaves remain excluded. The coronary venous scope includes the explicitly typed Vieussens valve component. No chamber, myocardium, cardiac valve apparatus or ostial reference surface is imported. Original source names, hierarchy, raw attributes and scope rationales remain preserved. Lymphatic geometry remains excluded.

The complete authority closure remains 1,178 terms: 1,149 male-applicable terms and 29 sex exclusions. Exact linguistic and source-hierarchy retrieval finds 482 distinct candidate identities (41.950%); 823 source instances have candidate identities, 10 are ambiguous and 891 remain unmapped. The potential identity union with current direct references is 616/1,149 (53.612%); 667 male terms lack a donor candidate identity. These are identity candidates, with **zero anatomically accepted terms**. Names and candidate matches establish no wall, lumen, patency, junction, fit, packing, completeness or redistribution acceptance. Openings are not solid plugs, remnant tissue is not a patent artery, and source casts do not establish wall layers. Local review assets and all original target/source geometry remain unchanged by this inventory refresh.

Corrected evidence: `.local-models/geometry-workbench/biodigital-male-circulation-v1/ta98-inventory-v2/inventory.json`, SHA256 `{sha(path)}`. Explicit source scope is `.local-models/geometry-workbench/biodigital-male-circulation-v1/source-v2/source-scope-policy.json`. Historical v1 remains immutable at `ta98-inventory/inventory.json`, SHA256 `{historical_sha}`; its earlier 1,699-instance / 466-candidate count omitted vascular leaves under Heart. The 90% identity target is not achieved.

| TA98 | Parent | Term | Planned representation | Source instance IDs / missing | Status |
|---|---|---|---|---|---|
'''
 lines=[]
 for row in terms:
  cells=[row['id'],row['parent_id'],row['name_en'],row['planned_representation_class'],'; '.join(row['donor_candidate_source_ids']) if row['donor_candidate_source_ids'] else 'missing distinct donor identity',row['donor_status']]
  lines.append('| '+' | '.join(str(c or '').replace('|','\\|').replace('\n',' ') for c in cells)+' |')
 section=intro+'\n'.join(lines)+'\n\n';updated=text[:start]+section+text[end:];assert updated[:start]==text[:start] and updated[start+len(section):]==text[end:];review.write_text(updated);outside_section_unchanged=True
assert sha(historical)==historical_sha
report={'status':'pass-corrected-source-identity-inventory','inventorySha256':sha(path),'historicalV1Sha256Unchanged':historical_sha,'sourceInstances':1724,'sharedSourceAssets':1165,'authorityTerms':1178,'maleApplicableTerms':1149,'candidateIdentityTerms':len(covered),'missingDonorIdentityTerms':summary['donor_missing_ta98_terms'],'candidateIdentityInstances':summary['source_instance_status']['candidate-identity'],'ambiguousIdentityInstances':10,'unmappedInstances':summary['source_instance_status']['unmapped'],'potentialUnionIdentityTerms':summary['potential_union_identity_reference_terms'],'crossRootCandidateInstancesOutsideDenominator':summary.get('cross_root_candidate_instances_outside_denominator',0),'crossRootCandidateDistinctTermsOutsideDenominator':summary.get('cross_root_candidate_distinct_terms_outside_denominator',0),'acceptedAnatomyTerms':0,'exactGeometryBindingsInstalled':0,'auditedSourceHeartLeaves':46,'includedSourceHeartLeaves':25,'excludedSourceHeartLeaves':21,'onlyRequestedReviewSectionUpdated':outside_section_unchanged,'geometryOrInstalledCatalogueMutated':False,'elapsedSeconds':time.monotonic()-started}
(args.inventory_directory/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
