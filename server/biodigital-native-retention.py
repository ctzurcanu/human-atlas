"""Plan selective native vessel references absent from qualified donor identity.

No restoration, catalogue mutation, geometry alteration, fitting or acceptance.
"""
import collections,hashlib,json,re,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
started=time.monotonic();work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';out=work/'native-retention-plan-v1'
if out.exists():raise SystemExit('Refusing to overwrite native retention evidence')
master_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';inventory_path=work/'ta98-inventory-v4/inventory.json';resolved_path=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3/proposal-v2/resolved-installed.json';profile_path=work/'atlas-review-v3/atlas.json'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
master_sha=sha(master_path);assert master_sha=='46262d8ab07682fcae3ff0ee634f652abad766c1e010f6e2ec311755aea9abf5'
master=json.loads(master_path.read_bytes());inventory=json.loads(inventory_path.read_bytes());resolved=json.loads(resolved_path.read_bytes());assert resolved['sourceSha256']==master_sha
profile=json.loads(profile_path.read_bytes());profile_parts={p['id']:p for p in profile['parts']};parts={p['id']:p for p in master['parts']};concepts={c['id']:c for c in master['concepts']};hidden=set(master.get('defaultHidden',[]));profile_hidden=set(profile.get('defaultHidden',[]))
resolved_exact=collections.defaultdict(set)
for b in resolved['conceptBindings']:
 if b.get('kind')=='exact':
  for pid in b['partIds']:resolved_exact[(b['conceptId'],pid)].add(b['term'])
out.mkdir();(out/'checkpoint.json').write_text(json.dumps({'status':'inputs-bound','masterSha256':master_sha,'donorInventorySha256':sha(inventory_path)})+'\n')
reader=AtlasReader(ROOT);candidates=[];rejected=[];raw_context=[];geometry_groups=collections.defaultdict(list)
donor_terms={r['id'] for r in inventory['terms'] if r['donor_candidate_source_ids']};missing=[r for r in inventory['terms'] if r['male_applicable'] and not r['donor_candidate_source_ids'] and r['existing_installed_direct_part_ids']]
for term in missing:
 code=term['id']
 for pid in term['existing_installed_direct_part_ids']:
  p=parts.get(pid);reason=None
  if p is None:reason='No actual current master part'
  elif term['planned_representation_class']!='vessel-or-vascular-network' or code.startswith('A12.0.'):reason='Generic container, fluid, layer, valve or other non-vessel scope; do not count a vessel reference as this whole identity'
  elif p['system'] not in ['arterial','venous']:reason='Reference part is not explicitly arterial or venous; reject generic root/digestive/cardiac reference inflation'
  elif p.get('suppressed') or p.get('defaultHidden') or pid in hidden:reason='Not the preferred active native master representation'
  elif pid not in profile_parts:reason='Part absent from current reduced-body review; needs separately justified delivery rather than blind restoration'
  elif not (pid.startswith('ZA:') or p.get('provenance',{}).get('url')=='https://github.com/Z-Anatomy/Models-of-human-anatomy'):reason='Independent donor/reference frame rather than preferred native source; not automatically interchangeable'
  else:
   c=concepts.get(p['conceptId'],{})
   if p.get('ta98Term') and p['ta98Term']!=code:reason='Literal part identity belongs to a different term; ancestor/parent alias is not new distinct geometry'
   elif c.get('ta98Term') and (c.get('ta98Kind')!='exact' or c['ta98Term']!=code):reason='Current own concept has a different identity or only parent placement'
   elif c.get('ta98Term')!=code and code not in resolved_exact.get((p['conceptId'],pid),set()):reason='No exact own-concept/part correspondence; census ancestry alone is insufficient'
   elif not re.search(r'arter|aort|vein|venous|branch|trunk|arteriae|vena',p['name'],re.I):reason='Source name does not identify vascular tissue; keep uncertain generic reference deferred'
  if reason:
   rejected.append({'term':code,'termName':term['name_en'],'partId':pid,'sourceName':p.get('name') if p else None,'reason':reason});continue
  mesh=reader.read_part(master,p,frame='model:local-ta98:assembled-metres');c=concepts[p['conceptId']]
  geometry_sha=mesh.geometry_sha;attribute_sha=mesh.attributes_sha;material_sha=hashlib.sha256(json.dumps(master.get('materials',{}).get(p.get('material')),sort_keys=True,separators=(',',':')).encode()).hexdigest();combined=hashlib.sha256((geometry_sha+attribute_sha+material_sha).encode()).hexdigest()
  native_identity='authored-exact-native-concept' if c.get('ta98Term')==code and c.get('ta98Kind')=='exact' else 'current-resolved-exact-own-concept-source'
  region=bool(p.get('geometryInstance') or p.get('sourceId')!=pid.removeprefix('ZA:') or pid.startswith('TA98PART:'))
  row={'term':code,'termEnglish':term['name_en'],'termLatin':term['name_la'],'authorityParent':{'id':term['parent_id'],'name':term['parent_name']},'partId':pid,'conceptId':p['conceptId'],'name':p['name'],'sourceId':p.get('sourceId'),'system':p['system'],'identityEvidence':native_identity,'sourceConcept':c,'exactResolvedOwnConceptTerms':sorted(resolved_exact.get((p['conceptId'],pid),set())),'meshIdentity':{'worldGeometrySha256':geometry_sha,'attributesSha256':attribute_sha,'materialSha256':material_sha,'combinedRenderIdentitySha256':combined,'bankPath':mesh.source['bank_path'],'bankSha256':mesh.source['bank_sha256'],'vertices':len(mesh.vertices),'faces':len(mesh.faces),'boundsMeters':mesh.bounds.tolist()},'nativeMasterPreferred':True,'reviewV3Suppressed':bool(profile_parts[pid].get('suppressed')),'reviewV3DefaultHidden':bool(profile_parts[pid].get('defaultHidden') or pid in profile_hidden),'potentialNativeRegionOverlay':region,'catalogueDeclaredFiniteWall':bool(p.get('sectionSolid')),'existingCoverageQualification':p.get('coverageLimitation') or c.get('coverageLimitation'),'donorAbsenceQualification':'No qualified donor identity for this TA98 term in inventory-v4; absence of the actual donor geometry is NOT established. Unmapped or combined donor meshes may still contain the vessel/region.','retentionStatus':'candidate-native-reference-for-selective-review-only','restoreAuthorized':False,'anatomicalAcceptance':False,'interfaceCaveat':'Native source is in the fixed current body frame; the BioDigital assembly has a different registered/deformed course. Before visibility restoration, verify source versus native ownership, exact duplicate/near overlap, reciprocal inlet/outlet correspondence and tissue contacts. Native region pieces can overlap the donor whole vessel and must not be blindly added.','geometryOrSourceModified':False}
  candidates.append(row);geometry_groups[combined].append({'term':code,'partId':pid})
  if len(candidates)%40==0:(out/'checkpoint.json').write_text(json.dumps({'status':'native-mesh-identities-saved-in-progress','candidateParts':len(candidates),'elapsedSeconds':time.monotonic()-started})+'\n')
duplicate_groups=[v for v in geometry_groups.values() if len(v)>1]
duplicate_ids={row['partId'] for group in duplicate_groups for row in group}
for row in candidates:row['samePlacedRenderIdentityHasOtherCandidateOwners']=row['partId'] in duplicate_ids
summary={'missingDonorTermsWithHistoricalNativeReferences':len(missing),'candidateNativeParts':len(candidates),'candidateDistinctMissingTerms':len({r['term'] for r in candidates}),'uniquePlacedRenderIdentities':len(geometry_groups),'sharedIdentityGroupsRequiringOwnershipReview':len(duplicate_groups),'rejectedOrDeferredReferences':len(rejected),'candidatePotentialRegionOverlays':sum(r['potentialNativeRegionOverlay'] for r in candidates),'candidateDeclaredFiniteParts':sum(r['catalogueDeclaredFiniteWall'] for r in candidates),'anatomicallyAcceptedTerms':0,'restoredParts':0,'elapsedSeconds':time.monotonic()-started}
plan={'schema':'selective-native-blood-retention-plan/v1','status':'review-plan-only-no-restoration','inputs':{str(p):sha(p) for p in [master_path,inventory_path,resolved_path,profile_path]},'summary':summary,'candidates':candidates,'identicalPlacedRenderIdentityGroups':duplicate_groups,'rejectedOrDeferred':rejected,'selectionPolicy':'Only current preferred native arterial/venous parts, with literal own Part/Concept or current exact own-concept identity matching the missing term. Reject generic root or borrowed ancestor references. All restoration remains conditional on donor/native geometry and interface review.','profileGeometryAndVisibilityChanged':False,'localOnly':True}
(out/'plan.json').write_text(json.dumps(plan,indent=2)+'\n');(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(out/'checkpoint.json').write_text(json.dumps({'status':plan['status'],'elapsedSeconds':summary['elapsedSeconds']})+'\n')
assert sha(master_path)==master_sha;print(json.dumps(summary))
