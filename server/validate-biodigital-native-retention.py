"""Read-only native retention identity/bank replay; no visibility restoration."""
import collections,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
started=time.monotonic();work=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';out=work/'native-retention-plan-v1';path=out/'plan.json';plan=json.loads(path.read_bytes());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for filename,digest in plan['inputs'].items():assert sha(filename)==digest,filename
master_path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';master=json.loads(master_path.read_bytes());parts={p['id']:p for p in master['parts']};concepts={c['id']:c for c in master['concepts']};hidden=set(master.get('defaultHidden',[]))
inventory=json.loads((work/'ta98-inventory-v4/inventory.json').read_bytes());terms={r['id']:r for r in inventory['terms']};donor={r['id'] for r in inventory['terms'] if r['donor_candidate_source_ids']};reader=AtlasReader(ROOT)
geometry_groups=collections.defaultdict(list);surface_groups=collections.defaultdict(list);bank_files={};categories=collections.Counter()
assert len(plan['candidates'])==126 and len({r['term'] for r in plan['candidates']})==74
for row in plan['candidates']:
 p=parts[row['partId']];c=concepts[p['conceptId']];code=row['term'];term=terms[code]
 assert code not in donor and not code.startswith('A12.0.') and p['system'] in ['arterial','venous'] and not(p.get('suppressed') or p.get('defaultHidden') or p['id'] in hidden)
 assert row['termEnglish']==term['name_en'] and row['termLatin']==term['name_la'] and row['authorityParent']=={'id':term['parent_id'],'name':term['parent_name']}
 assert not row['restoreAuthorized'] and not row['anatomicalAcceptance'] and not row['geometryOrSourceModified']
 if p.get('ta98Term'):assert p['ta98Term']==code
 if c.get('ta98Term'):assert c['ta98Term']==code and c['ta98Kind']=='exact'
 else:assert code in row['exactResolvedOwnConceptTerms']
 mesh=reader.read_part(master,p,frame='model:local-ta98:assembled-metres');identity=row['meshIdentity'];assert mesh.geometry_sha==identity['worldGeometrySha256'] and mesh.attributes_sha==identity['attributesSha256'] and mesh.source['bank_sha256']==identity['bankSha256']
 assert len(mesh.vertices)==identity['vertices'] and len(mesh.faces)==identity['faces']
 geometry_groups[mesh.geometry_sha].append(row['partId']);surface_groups[mesh.surface_sha].append(row['partId']);bank_files[identity['bankPath']]=identity['bankSha256']
 if row['catalogueDeclaredFiniteWall']:categories['catalogue-declared-finite-vessel-or-junction']+=1
 elif row['potentialNativeRegionOverlay']:categories['native-region-ownership-review-required']+=1
 else:categories['named-native-source-reference']+=1
duplicate_geometry=[ids for ids in geometry_groups.values() if len(ids)>1];duplicate_surface=[ids for ids in surface_groups.values() if len(ids)>1]
assert all(not r['term'].startswith('A12.0.') for r in plan['candidates'])
assert any(r['term']=='A12.0.00.000' for r in plan['rejectedOrDeferred'])
assert not plan['profileGeometryAndVisibilityChanged'] and plan['summary']['restoredParts']==0 and plan['summary']['anatomicallyAcceptedTerms']==0
result={'status':'pass-native-reference-identity-and-actual-mesh-replay-only','planSha256':sha(path),'masterSha256Unchanged':sha(master_path),'candidateParts':126,'candidateDistinctStillMissingTerms':74,'uniqueWorldGeometryHashes':len(geometry_groups),'uniqueOrientedSurfaceHashes':len(surface_groups),'exactGeometryDuplicateGroups':duplicate_geometry,'exactOrientedSurfaceDuplicateGroups':duplicate_surface,'boundSourceBanks':bank_files,'sourceReferenceCategories':dict(categories),'rejectedOrDeferredGenericBorrowedNonpreferredReferences':len(plan['rejectedOrDeferred']),'potentialNativeRegionOverlaysRequireDonorOwnershipExclusion':sum(r['potentialNativeRegionOverlay'] for r in plan['candidates']),'donorMissingIdentityDoesNotProveDonorGeometryAbsence':True,'interfacesAndNewGlobalRegistrationNotAccepted':True,'actualRetainedOrRestoredParts':0,'newRuntimeBindings':0,'geometryOrProfileModified':False,'elapsedSeconds':time.monotonic()-started}
(out/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='boundSourceBanks'}))
