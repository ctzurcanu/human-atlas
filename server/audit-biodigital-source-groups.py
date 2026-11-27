"""Evidence-only named source-group identity delta; no profile or geometry edits."""
import argparse,collections,hashlib,json,re,sqlite3,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.ta98_authority import corrected_terms

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def normal(s):return re.sub(r'\s+',' ',re.sub(r'[^a-z0-9]+',' ',re.sub(r'\b(?:left|right)\b','',s.lower().replace('-','')))).strip()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();out=args.output.resolve();began=time.monotonic()
 if (out/'group-audit.json').exists():raise SystemExit('Refusing to overwrite group evidence')
 work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';source=work/'source-v2'
 inventory_path=work/'ta98-inventory-v4/inventory.json';hierarchy_path=source/'original-cardiovascular-hierarchy.json';authority=ROOT/'.local-models/ta98-authority.sqlite';profile_path=work/'atlas-identity-v4/atlas.json'
 assert sha(inventory_path)=='1c809cc053afcb0d318eff11ec6269642a601885315c8ee35aaf147ff90ae7ba'
 paths=[inventory_path,hierarchy_path,source/'instances.json',source/'manifest.json',authority,profile_path,out/'ambiguity-audit.json'];inputs={str(p):sha(p) for p in paths}
 inv=json.loads(inventory_path.read_text());instances={r['id']:r for r in json.loads((source/'instances.json').read_text())['instances']};leaves={r['source_id']:r for r in inv['source_mesh_instances']}
 ambiguity=json.loads((out/'ambiguity-audit.json').read_text());resolved=[r for r in ambiguity['findings'] if r['status']=='firm-context-resolved-identity'];assert len(resolved)==7
 baseline={r['source_id']:r['ta98_id'] for r in inv['source_mesh_instances'] if r['status']=='candidate-identity'}
 baseline.update({r['source_id']:r['authority']['id'] for r in inv['cross_root_identity_candidates_outside_current_denominator']})
 term_by_source={**baseline,**{r['sourceId']:r['selectedCandidateTerm'] for r in resolved}}
 with sqlite3.connect('file:'+str(authority)+'?mode=ro',uri=True) as db:
  db.row_factory=sqlite3.Row;terms=corrected_terms({r['id']:dict(r) for r in db.execute('select * from ta98')});synonyms=[dict(r) for r in db.execute('select * from synonyms')]
 scope={t['id'] for t in inv['terms'] if t['male_applicable']};children=collections.defaultdict(set)
 for c,t in terms.items():children[t['parent_id']].add(c)
 def authority_row(c):return {k:terms[c][k] for k in ['id','name_en','name_la','parent_id']}
 def ancestors(c):
  result=[];seen=set()
  while c in terms and c not in seen:
   seen.add(c);result.append(authority_row(c));c=terms[c]['parent_id']
  return result
 nodes={};parents={};descendants={};raw_leaves={}
 def walk(n,path):
  key=n['objectId'];assert key not in nodes;nodes[key]=n;parents[key]=path;selected=[];raw=[]
  for child in n.get('objects',[]):
   d,r=walk(child,path+[key]);selected+=d;raw+=r
  if n.get('geometryId'):raw.append(key)
  if key in instances:selected.append(key)
  descendants[key]=selected;raw_leaves[key]=raw;return selected,raw
 walk(json.loads(hierarchy_path.read_text()),[])
 assert set(descendants[next(iter(nodes))])==set(instances)
 lex=collections.defaultdict(list)
 for code in scope:lex[normal(terms[code]['name_en'])].append(dict(term=code,label=terms[code]['name_en'],method='authority-English-label'))
 for r in synonyms:
  if r['id'] in scope and r['lang'].lower() in {'en','eng','english'}:lex[normal(r['synonym'])].append(dict(term=r['id'],label=r['synonym'],method='authority-English-synonym'))
 matches=[]
 for key,n in nodes.items():
  if n.get('geometryId') or len(descendants[key])<2:continue
  evidence=lex.get(normal(n['displayName']),[])
  if evidence:matches.append(dict(sourceGroupId=key,sourceGroupName=n['displayName'],sourcePath=[dict(id=p,name=nodes[p]['displayName']) for p in parents[key]],
   exactCandidateAuthorities=[authority_row(c) for c in sorted({e['term'] for e in evidence})],matchingEvidence=evidence,
   sourceDescendantIds=descendants[key],excludedNonvascularRawDescendantIds=sorted(set(raw_leaves[key])-set(instances)),
   newTermsRelativeToLeafInventory=sorted({e['term'] for e in evidence}-set(baseline.values()))))
 # One reviewed comparative synonym in its explicit named splenic branch context.
 great='human_20_male_cardiovascular_system-great_pancreatic_artery_ID';great_row=leaves[great]
 assert great_row['source_name']=='Great pancreatic artery' and great_row['source_parent_name']=='Pancreatic branches of splenic artery'
 assert terms['A12.2.12.045']['name_en']=='greater pancreatic artery' and terms['A12.2.12.045']['name_la']=='a. pancreatica magna' and terms['A12.2.12.045']['parent_id']=='A12.2.12.041'
 great_rule=dict(sourceId=great,sourceName=great_row['source_name'],sourceHierarchy=great_row['source_hierarchy'],selectedTerm='A12.2.12.045',authority=authority_row('A12.2.12.045'),authorityAncestorPath=ancestors('A12.2.12.045'),
  rule='Reviewed explicit great/greater pancreatic comparative synonym, under the named pancreatic branches of splenic artery; unique matching authority vessel and parent. No generic pancreas-vessel mapping.',
  identityOnly=True,anatomicalAcceptance=False)
 term_by_source[great]='A12.2.12.045'
 def membership(ids):
  return [dict(sourceId=s,partId='BDM:male:'+s,sourceName=instances[s]['name'],sourceGeometryId=instances[s]['geometryId'],sourcePath=[dict(id=p,name=nodes[p]['displayName']) for p in parents[s]],
   leafAuthorityIdentity=term_by_source.get(s),vascularRole=instances[s]['vascularRole']) for s in sorted(ids)]
 def candidate(name,term,group,anchor,ids,selection):
  assert len(ids)==len(set(ids)) and set(ids)<=set(descendants[anchor]) and set(ids)<=set(instances)
  core=children[term];represented={term_by_source.get(s) for s in ids};missing=core-represented
  assert core and not missing
  return dict(name=name,term=term,authority=authority_row(term),authorityAncestorPath=ancestors(term),sourceGroupId=group,sourceGroupName=nodes[group]['displayName'],sourceAnchorId=anchor,
   sourceAnchorName=nodes[anchor]['displayName'],selectionRule=selection,sourceGroupOriginalDescendantIds=descendants[group],selectedSourceDescendantIds=sorted(ids),members=membership(ids),
   sourceGroupDescendantsOmittedFromSelector=sorted(set(descendants[group])-set(ids)),additionalNamedDescendantsFromAnchor=sorted(set(ids)-set(descendants[group])),
   authoritativeDirectChildTerms=sorted(core),representedAuthoritativeDirectChildTerms=sorted(core&represented),missingAuthoritativeDirectChildTerms=[],
   fullEnumeratedAuthorityChildSetRepresented=True,anonymousBorrowedMembers=False,sourceGeometryCopied=False,sourceHierarchyUnchanged=True,
   candidateConceptId='BDM:male:identity-group:'+term,proposedRepresentation='source-identity-group',installedOrSupportedByCurrentResolver=False,
   coverageLimitation='Named source members establish this collective identity and the enumerated authority child set. Anatomical completeness, body fit, vessel walls and lumen continuity remain unverified.',anatomicalAcceptance=False)
 prefix='human_20_male_cardiovascular_system-';pulmonary=prefix+'pulmonary_veins_ID';hepatic=prefix+'hepatic_veins_ID';pancreatic=prefix+'pancreatic_branches_of_splenic_artery_ID';celiac=prefix+'celiac_trunk_ID'
 pulmonary_ids=descendants[pulmonary]
 hepatic_ids=[s for s in descendants[hepatic] if term_by_source.get(s) in children['A12.3.09.005']]
 pancreatic_ids=[s for s in descendants[celiac] if term_by_source.get(s) in children['A12.2.12.041']]
 assert len(pulmonary_ids)==27 and len(hepatic_ids)==3 and len(pancreatic_ids)==5
 candidates=[candidate('Pulmonary veins','A12.3.02.001',pulmonary,pulmonary,pulmonary_ids,'All27named vascular descendants of original pure pulmonary-veins group; all four named pulmonary trunks are present.'),
  candidate('Hepatic veins','A12.3.09.005',hepatic,hepatic,hepatic_ids,'Exactly right/intermediate/left hepatic trunks. Original organizer includes24portal/gastric/splenic descendants which are explicitly excluded from this derived selector.'),
  candidate('Pancreatic branches of splenic artery','A12.2.12.041',pancreatic,celiac,pancreatic_ids,'All five explicitly named authority children under common source celiac anchor. Three are in pancreatic organizer; inferior pancreatic and tail arteries are named source siblings, not anonymous borrowed geometry.')]
 regional_codes={'A12.3.08.001','A12.3.11.001','A12.2.07.001','A12.2.09.001','A12.2.16.001'}
 regional=[dict(sourceGroupId=r['sourceGroupId'],sourceGroupName=r['sourceGroupName'],sourceDescendantIds=r['sourceDescendantIds'],candidateAuthorities=[t for t in r['exactCandidateAuthorities'] if t['id'] in regional_codes],
  status='qualified-regional-selection-control-only-not-detailed-coverage',anatomicalCompletenessAsserted=False) for r in matches if any(t['id'] in regional_codes for t in r['exactCandidateAuthorities'])]
 def deferred_group(suffix,code,reason):
  key=prefix+suffix+'_ID';represented={term_by_source.get(s) for s in descendants[key]}
  return dict(sourceGroupId=key,sourceGroupName=nodes[key]['displayName'],authority=authority_row(code),sourceDescendantIds=descendants[key],missingEnumeratedAuthorityChildren=sorted(children[code]-represented),reason=reason,status='not-counted-as-complete-group-identity')
 deferred=[deferred_group('emissary_veins','A12.3.05.301','Source group lacks explicitly established foramen-ovale and internal-carotid plexus identities. Sphenoidal emissary veins cannot be silently relabeled as a plexus.'),
  deferred_group('veins_of_spinal_cord','A12.3.07.023','Mixed source organizer includes lateral sacral, intervertebral and vertebral-plexus vessels. Do not equate every descendant with spinal-cord veins.'),
  deferred_group('lumbar_veins','A12.3.09.003','Exact authority homonym occurs under ascending lumbar and inferior vena cava branches. Source systemic-abdominal organizer does not establish the required drainage distinction.'),
  dict(sourceGroupId=prefix+'cardiovascular_system_ID',status='excluded-generic-root-and-cardiac-mixture',reason='Generic root is not an additional detailed blood-vessel identity; includes excluded cardiac anatomy.')]
 baseline_terms=set(baseline.values());new_leaf_terms=({r['selectedCandidateTerm'] for r in resolved}|{'A12.2.12.045'})-baseline_terms;group_terms={r['term'] for r in candidates}-baseline_terms
 summary=dict(rawOriginalHierarchyNodes=len(nodes),selectedSourceInducedHierarchyNodes=len(set().union(*(set(i['lineage']) for i in instances.values()))),selectedVascularInstances=len(instances),
  eligibleMultiMemberSourceGroups=sum(not n.get('geometryId') and len(descendants[k])>1 for k,n in nodes.items()),strictAuthorityGroupMatches=len(matches),
  proposedDetailedCollectiveGroups=len(candidates),newDetailedGroupTerms=sorted(group_terms),groupMembershipReferences=sum(len(c['members']) for c in candidates),
  firmAmbiguityResolvedInstances=7,newContextLeafInstances=1,newLeafIdentityTerms=sorted(new_leaf_terms),conditionalAmbiguitiesNotAdopted=2,unresolvedAmbiguities=1,
  regionalControlsSeparateNotDetailedCoverage=len(regional),distinctRegionalTermsNotCounted=len(regional_codes),baselineHistoricalLeafTerms=550,baselineCrossRootTerms=4,
  potentialHistoricalIdentityTermsIfDeltaReviewed=550+len(new_leaf_terms|group_terms),potentialDistinctIdentityTermsIncludingCrossRoot=554+len(new_leaf_terms|group_terms),
  newGeometry=0,installedAdditionalBindings=0,defaultProfileChanged=False,anatomicalAcceptance=False)
 result=dict(schema='biodigital-named-source-group-identity-delta/v1',status='candidate-evidence-only',inputs=inputs,summary=summary,
  candidateGroups=candidates,additionalReviewedLeafRules=[great_rule],resolvedAmbiguityEvidence=dict(path=str(out/'ambiguity-audit.json'),sha256=sha(out/'ambiguity-audit.json')),
  qualifiedRegionalSelectionCandidates=regional,deferredOrExcludedGroups=deferred,strictGroupNameProspects=matches,
  resolverIntegrationConstraint='Current source-identity controls require each Part conceptId and ta98Term to equal the group term/concept. These proposed group controls require an explicit guarded aggregate representation that preserves exact leaf identities; do not relabel/copy parts to bypass it.',
  sourceHierarchyAndGeometryAndCurrentProfileUnchanged=True,elapsedSeconds=time.monotonic()-began)
 for p,digest in inputs.items():assert sha(p)==digest
 (out/'group-audit.json').write_text(json.dumps(result,indent=2)+'\n');(out/'checkpoint.json').write_text(json.dumps(dict(status='source-group-evidence-prepared',elapsedSeconds=result['elapsedSeconds']))+'\n')
 print(json.dumps(dict(output=str(out/'group-audit.json'),summary=summary,elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
