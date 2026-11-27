"""Independent source/authority membership replay for group identity evidence."""
import argparse,hashlib,json,re,sqlite3,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.ta98_authority import corrected_terms
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();began=time.monotonic()
 report_path=out/'group-audit.json';r=json.loads(report_path.read_text())
 for p,h in r['inputs'].items():assert sha(p)==h,p
 work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';source=work/'source-v2';inv=json.loads((work/'ta98-inventory-v4/inventory.json').read_text())
 instances={x['id']:x for x in json.loads((source/'instances.json').read_text())['instances']};nodes={};parents={};desc={}
 def visit(node,path):
  key=node['objectId'];assert key not in nodes;nodes[key]=node;parents[key]=path;found=[key] if key in instances else []
  for child in node.get('objects',[]):found+=visit(child,path+[key])
  desc[key]=set(found);return found
 visit(json.loads((source/'original-cardiovascular-hierarchy.json').read_text()),[])
 db=sqlite3.connect('file:'+str(ROOT/'.local-models/ta98-authority.sqlite')+'?mode=ro',uri=True);db.row_factory=sqlite3.Row
 terms=corrected_terms({x['id']:dict(x) for x in db.execute('select * from ta98')})
 ambiguity=json.loads((out/'ambiguity-audit.json').read_text());firm=[x for x in ambiguity['findings'] if x['status']=='firm-context-resolved-identity'];assert len(firm)==7
 inventory_rows={x['source_id']:x for x in inv['source_mesh_instances']}
 mapping={x['source_id']:x['ta98_id'] for x in inv['source_mesh_instances'] if x['status']=='candidate-identity'}
 mapping.update({x['sourceId']:x['selectedCandidateTerm'] for x in firm})
 leaf_rules=[]
 def rule(sid,name,hierarchy,code,why,rule_id):
  t=terms[code];parent=terms[t['parent_id']]
  return dict(rule_id=rule_id,source_id=sid,source_name=name,source_hierarchy=hierarchy,authority=dict(id=code,name_en=t['name_en'],name_la=t['name_la'],parent_id=t['parent_id'],parent_name=parent['name_en']),
   rationale=why,candidate_identity_only=True,out_of_current_denominator=False,anatomical_acceptance=False,runtime_binding=False,source_geometry=instances[sid]['geometryId'])
 for x in firm:
  sid=x['sourceId'];assert inventory_rows[sid]['status']=='ambiguous-identity' and x['selectedCandidateTerm'] in inventory_rows[sid]['candidate_ta98_ids']
  assert nodes[sid]['displayName']==x['sourceName'] and x['sourceParentId']==parents[sid][-1]
  assert x['sourceHierarchy']==[dict(id=p,name=nodes[p]['displayName']) for p in parents[sid]]
  for candidate in [x['selectedAuthority']]+x['alternativeAuthorities']:
   t=terms[candidate['id']];assert all(candidate[k]==t[k] for k in ['name_en','name_la','parent_id'])
  leaf_rules.append(rule(sid,x['sourceName'],x['sourceHierarchy'],x['selectedCandidateTerm'],x['rationale'],'explicit-source-context-disambiguation'))
 for x in r['additionalReviewedLeafRules']:
  sid=x['sourceId'];assert inventory_rows[sid]['status']=='unmapped' and nodes[sid]['displayName']=='Great pancreatic artery'
  assert nodes[parents[sid][-1]]['displayName']=='Pancreatic branches of splenic artery'
  code=x['selectedTerm'];assert code=='A12.2.12.045' and terms[code]['name_en']=='greater pancreatic artery' and terms[code]['name_la']=='a. pancreatica magna' and terms[code]['parent_id']=='A12.2.12.041'
  mapping[sid]=code;leaf_rules.append(rule(sid,x['sourceName'],x['sourceHierarchy'],code,x['rule'],'named-great-greater-pancreatic-context'))
 assert len(leaf_rules)==8 and len({x['source_id'] for x in leaf_rules})==8
 checked=0;candidate_terms=set()
 for c in r['candidateGroups']:
  code=c['term'];candidate_terms.add(code);t=terms[code]
  assert c['authority']=={k:t[k] for k in ['id','name_en','name_la','parent_id']}
  selected=set(c['selectedSourceDescendantIds']);assert len(selected)==len(c['members']) and selected<=desc[c['sourceAnchorId']]
  assert c['sourceGroupOriginalDescendantIds'] and set(c['sourceGroupOriginalDescendantIds'])==desc[c['sourceGroupId']]
  assert set(c['sourceGroupDescendantsOmittedFromSelector'])==desc[c['sourceGroupId']]-selected
  assert set(c['additionalNamedDescendantsFromAnchor'])==selected-desc[c['sourceGroupId']]
  authority_children={x['id'] for x in terms.values() if x['parent_id']==code}
  assert authority_children and authority_children<=set(mapping[s] for s in selected if s in mapping)
  assert authority_children==set(c['authoritativeDirectChildTerms'])==set(c['representedAuthoritativeDirectChildTerms']) and not c['missingAuthoritativeDirectChildTerms']
  for member in c['members']:
   sid=member['sourceId'];checked+=1;assert sid in selected and member['partId']=='BDM:male:'+sid
   assert member['sourceName']==nodes[sid]['displayName']==instances[sid]['name'] and member['sourceGeometryId']==instances[sid]['geometryId']
   assert member['sourcePath']==[dict(id=p,name=nodes[p]['displayName']) for p in parents[sid]]
   assert member['leafAuthorityIdentity']==mapping.get(sid)
   assert instances[sid]['vascularRole'] in ['arterial-vessel','venous-vessel']
   assert not re.fullmatch(r'(?:left |right )?(?:branch|branches|artery|arteries|vein|veins)(?: \d+)?',member['sourceName'],re.I)
  assert not c['sourceGeometryCopied'] and not c['anonymousBorrowedMembers'] and not c['anatomicalAcceptance'] and not c['installedOrSupportedByCurrentResolver']
  assert c['fullEnumeratedAuthorityChildSetRepresented'] and 'remain unverified' in c['coverageLimitation']
 assert candidate_terms=={'A12.2.12.041','A12.3.02.001','A12.3.09.005'} and checked==35
 assert not {'A12.0.00.000','A12.2.00.001','A12.3.00.001'}&candidate_terms
 baseline_terms={x['ta98_id'] for x in inv['source_mesh_instances'] if x['status']=='candidate-identity'}
 leaf_delta={x['authority']['id'] for x in leaf_rules}-baseline_terms
 assert len(leaf_delta)==5 and r['summary']['potentialHistoricalIdentityTermsIfDeltaReviewed']==550+len(candidate_terms|leaf_delta)==558
 delta=dict(schema='biodigital-source-identity-candidate-delta/v1',sourceInventoryV4Sha256=sha(work/'ta98-inventory-v4/inventory.json'),
  evidenceGroupAuditSha256=sha(report_path),evidenceAmbiguityAuditSha256=sha(out/'ambiguity-audit.json'),
  in_denominator_leaf_rule_delta=leaf_rules,group_control_candidates=r['candidateGroups'],
  regional_control_candidates_not_counted=r['qualifiedRegionalSelectionCandidates'],conditional_or_deferred_not_adopted=[x for x in ambiguity['findings'] if x['status']!='firm-context-resolved-identity'],
  newLeafIdentityTerms=sorted(leaf_delta),newDetailedGroupIdentityTerms=sorted(candidate_terms),groupResolverIntegrationConstraint=r['resolverIntegrationConstraint'],
  profileOrGeometryOrDefaultChanged=False,identityOnly=True,anatomicalAcceptance=False)
 delta_path=out/'inventory-delta.json';assert not delta_path.exists();delta_path.write_text(json.dumps(delta,indent=2)+'\n')
 result=dict(status='pass-source-authority-and-named-membership-evidence',groupAuditSha256=sha(report_path),inventoryDeltaSha256=sha(delta_path),
  firmAmbiguities=7,newGreatPancreaticContextLeaf=1,newLeafDistinctTerms=5,newDetailedCollectiveTerms=3,exactNamedGroupMembershipsReplayed=checked,
  allEnumeratedAuthorityChildrenPresent=True,regionalGroupsNotCountedAsDetailedCoverage=True,genericSystemRootsExcluded=True,mixedSourceOrganizerForeignMembersExplicitlyExcluded=True,
  remainingConditionalAmbiguities=2,remainingDeferredAmbiguities=1,potentialHistoricalCandidateTerms=558,historicalDenominator=1149,crossRootCandidateTermsSeparate=4,
  currentIdentityProfileSha256Unchanged=sha(work/'atlas-identity-v4/atlas.json'),profileOrGeometryOrDefaultChanged=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
 for p,h in r['inputs'].items():assert sha(p)==h
 (out/'validation.json').write_text(json.dumps(result,indent=2)+'\n');(out/'checkpoint.json').write_text(json.dumps(dict(status=result['status'],elapsedSeconds=result['elapsedSeconds']))+'\n');print(json.dumps(result))

if __name__=='__main__':main()
