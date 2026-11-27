"""Bounded explicit context audit of the ten v4 ambiguous source vessel names.

Writes evidence only. No source, profile, runtime binding or geometry mutation.
"""
import argparse, collections, hashlib, json, sqlite3, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.ta98_authority import corrected_terms
WORK=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
INVENTORY_SHA='1c809cc053afcb0d318eff11ec6269642a601885315c8ee35aaf147ff90ae7ba'

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args();started=time.monotonic()
 output=args.output.resolve();assert output.parent.is_dir() and not output.exists()
 inv_path=WORK/'ta98-inventory-v4/inventory.json';assert sha(inv_path)==INVENTORY_SHA
 inventory=json.loads(inv_path.read_bytes());amb=[r for r in inventory['source_mesh_instances'] if r['status']=='ambiguous-identity'];assert len(amb)==10
 hierarchy_path=WORK/'source-v2/original-cardiovascular-hierarchy.json';graph=json.loads(hierarchy_path.read_bytes());nodes={};ancestry={}
 def walk(n,parents):
  nodes[n['objectId']]=n;ancestry[n['objectId']]=parents
  for child in n.get('objects',[]):walk(child,parents+[dict(id=n['objectId'],name=n['displayName'])])
 walk(graph,[])
 authority_path=ROOT/'.local-models/ta98-authority.sqlite'
 with sqlite3.connect('file:'+str(authority_path)+'?mode=ro',uri=True) as db:
  db.row_factory=sqlite3.Row;terms=corrected_terms({r['id']:dict(r) for r in db.execute('select * from ta98')})
  synonyms=collections.defaultdict(list)
  for r in db.execute('select id,synonym,lang from synonyms'):synonyms[r['id']].append(dict(r))
 def authority(code):
  term=terms[code];parents=[];current=term['parent_id'];seen=set()
  while current in terms and current not in seen:
   seen.add(current);r=terms[current];parents.append({k:r[k] for k in ['id','name_en','name_la','parent_id']});current=r['parent_id']
  return {k:term[k] for k in ['id','name_en','name_la','parent_id','parent_name']}|dict(ancestorPath=parents,synonyms=synonyms[code])
 findings=[]
 for row in amb:
  sid=row['source_id'];name=row['source_name'];assert nodes[sid]['displayName']==name and ancestry[sid]==row['source_hierarchy']
  hierarchy=[x['name'] for x in row['source_hierarchy']];chosen=None;status=None;rationale=None;qualification=None
  if 'calcaneal_branches_of_' in sid:
   assert 'posterior tibial artery' in name.lower() and row['source_parent_name'].startswith('Branches of ') and 'posterior tibial artery' in row['source_parent_name'].lower()
   chosen='A12.2.16.059';assert terms[chosen]['parent_id']=='A12.2.16.055' and terms['A12.2.16.075']['parent_id']=='A12.2.16.074'
   status='firm-context-resolved-identity';rationale='The source explicitly names calcaneal branches OF posterior tibial artery and places them immediately under that artery branch group. The homonym A12.2.16.075 is under lateral malleolar branches OF fibular artery; an ultimate posterior-tibial ancestor is insufficient to select that more specific alternative.'
  elif 'dorsal_digital_arteries' in sid:
   assert 'Arteries of upper limb' in hierarchy and 'hand and wrist' in row['source_parent_name']
   chosen='A12.2.09.035';assert terms[chosen]['parent_id']=='A12.2.09.032' and terms['A12.2.16.053']['parent_id']=='A12.2.16.048'
   status='firm-context-resolved-identity';rationale='The exact dorsal-digital-arteries name occurs in the explicit upper-limb, hand-and-wrist source hierarchy. The alternative A12.2.16.053 is under dorsalis pedis artery in the lower-limb authority subtree. The source group does not establish a separate radial-origin interface, which remains anatomically unverified.'
  elif 'intercapitular_veins' in sid:
   assert 'Veins of upper limb' in hierarchy and 'hand and wrist' in row['source_parent_name']
   chosen='A12.3.08.024';assert terms[chosen]['parent_id']=='A12.3.08.023' and terms['A12.3.11.019']['parent_id']=='A12.3.11.010'
   status='firm-context-resolved-identity';rationale='The source explicitly places these named intercapitular veins in the hand-and-wrist/upper-limb hierarchy. The competing term A12.3.11.019 is under the small saphenous vein in the lower limb. Source palmar organizational placement is retained and does not assert an accepted dorsal-network connection.'
  elif sid.endswith('-left_hepatic_vein_ID'):
   chosen='A12.3.09.008';assert terms[chosen]['name_en']=='left hepatic vein' and terms[chosen]['name_la']=='v. hepatica sinistra'
   assert terms['A12.3.09.007']['name_en']=='intermediate hepatic vein'
   status='firm-context-resolved-identity';rationale='The source name and parent Hepatic veins match the authoritative left hepatic vein (v. hepatica sinistra). The competing A12.3.09.007 is the intermediate hepatic vein; only an inherited installed-concept alias called it left, so that borrowed alias is rejected as authority evidence.'
  elif 'dorsal_scapular_artery' in sid:
   chosen='A12.2.08.058';assert terms[chosen]['name_en']=='dorsal scapular artery' and terms['A12.2.08.057']['name_en']=='deep branch'
   assert row['source_parent_name'].endswith('subclavian artery')
   status='conditional-canonical-identity-origin-distinction-unresolved';rationale='The source exactly names dorsal scapular artery, matching canonical A12.2.08.058. The competing A12.2.08.057 is the deep branch OF transverse cervical artery and has dorsal scapular artery as a synonym. The source only groups this mesh under branches of subclavian artery, so source ancestry does not independently resolve the transverse-cervical-versus-separate-artery origin distinction. Do not count as a firm new resolved term solely by preferring canonical labels.'
   qualification='Canonical source-name candidate only; branch-origin distinction remains unresolved. No accepted arterial origin or body-fit claim.'
  elif sid.endswith('-pancreatic_veins_ID'):
   assert row['source_parent_name']=='Hepatic portal system' and set(row['candidate_ta98_ids'])=={'A12.3.12.022','A12.3.12.029'}
   status='deferred-drainage-context-insufficient';rationale='Both authoritative pancreatic-vein terms are exact English/Latin homonyms, distinguished by drainage parent: superior mesenteric vein versus splenic vein. The source places the mesh only under Hepatic portal system, with Splenic vein as a sibling, so no parent-specific drainage identity is established. Sibling presence is not ancestry or evidence that this entire source mesh represents either or both collections.'
  else:raise AssertionError('Unexpected ambiguous source '+sid)
  if chosen:assert chosen in row['candidate_ta98_ids']
  collective=bool(chosen and terms[chosen].get('composite_property'))
  findings.append(dict(sourceId=sid,sourceName=name,sourceHierarchy=row['source_hierarchy'],sourceParentId=row['source_parent_id'],sourceParentName=row['source_parent_name'],sourceGeometryId=row['source_geometry_canonical_id'],
   priorStatus=row['status'],priorCandidateIds=row['candidate_ta98_ids'],status=status,selectedCandidateTerm=chosen,selectedAuthority=authority(chosen) if chosen else None,
   alternativeAuthorities=[authority(c) for c in row['candidate_ta98_ids'] if c!=chosen],rationale=rationale,sourceContextFactsVerifiedAgainstLiteralGraph=True,
   wholeCollectiveRepresentationAsserted=False,collectiveTerm=collective,identityQualification=qualification or ('Source member(s) identity only; complete collective, body fit and lumen continuity remain unverified.' if collective else 'Source identity only; body fit, wall completeness and lumen continuity remain unverified.'),
   geometryOrFitAccepted=False,runtimeBindingChanged=False,inventoryV4Changed=False))
 firm=[r for r in findings if r['status']=='firm-context-resolved-identity'];conditional=[r for r in findings if r['status'].startswith('conditional-')]
 current={r['ta98_id'] for r in inventory['source_mesh_instances'] if r['status']=='candidate-identity'}
 result=dict(schema='biodigital-v4-ambiguous-source-context-audit/v1',status='complete-bounded-evidence-only',inputs={str(inv_path):sha(inv_path),str(hierarchy_path):sha(hierarchy_path),str(authority_path):sha(authority_path)},
  reviewedAmbiguousInstances=len(findings),firmResolvedInstances=len(firm),firmResolvedDistinctTerms=sorted({r['selectedCandidateTerm'] for r in firm}),newFirmTermsRelativeToV4=sorted({r['selectedCandidateTerm'] for r in firm}-current),
  conditionalCanonicalInstances=len(conditional),conditionalCanonicalDistinctTerms=sorted({r['selectedCandidateTerm'] for r in conditional}),deferredInstances=sum(r['status'].startswith('deferred-') for r in findings),
  findings=findings,sourceOrGeometryOrProfileMutated=False,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-started)
 output.write_text(json.dumps(result,indent=2)+'\n');assert sha(inv_path)==INVENTORY_SHA
 print(json.dumps({k:v for k,v in result.items() if k not in ['findings','inputs']}))

if __name__=='__main__':main()
