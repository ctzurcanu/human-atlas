"""Aggregate saved sampled strict crossings; no new spatial queries or edits."""
from pathlib import Path
from collections import Counter,defaultdict
import hashlib,json,re,time
import numpy as np
start=time.monotonic();R=Path(__file__).resolve().parents[1]
B=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1';D=B/'regional-fit-v1/skeleton-samples-v1';O=D/'aggregation-v1'
if O.exists():raise SystemExit('Refusing to overwrite aggregation-v1')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
context=B/'regional-interface-v1/fixed-context';candidate=B/'regional-fit-v1/candidate-v1'
owners=json.loads((context/'context.json').read_bytes())['partOwnership'];fixed=np.load(context/'full.npz',allow_pickle=False)
bone_owner=fixed['face_owner_part_index'];bone_source_face=fixed['face_owner_source_face_id']
insts=json.loads((B/'source-v2/instances.json').read_bytes());instances={x['id']:x for x in insts['instances']}
parts=json.loads((candidate/'instance-ownership.json').read_bytes())['instances'];report=json.loads((D/'report.json').read_bytes());meta=[]
inputs=[D/'report.json',context/'context.json',context/'full.npz',candidate/'instance-ownership.json',B/'source-v2/instances.json']
for system in ['arterial','venous','arterial-remnant-companion']:
 path=candidate/(system+'-ancestry.npz');inputs.append(path);a=np.load(path,allow_pickle=False);ancestry_faces=a['candidate_face_source_face']
 for row in parts:
  if row['candidateMesh']!=system+'.npz':continue
  count=row['candidateFaceCount'];ids=row['candidateFaceStart']+np.unique(np.linspace(0,count-1,min(32,count)).astype(int))
  for i,s in zip(ids,ancestry_faces[ids]):meta.append(dict(instance=row['sourceInstanceId'],name=row['name'],candidateMesh=system,candidateFace=int(i),sourceFace=int(s),sampleIndex=len(meta)))
def family(name):
 s=name.lower()
 for label,terms in [
  ('intracranial-cerebral',('cerebr','pericallosal','callosomarginal','cingular','precuneal','choroid','lenticulostriate','frontobasal','frontopolar','sagittal sinus','cavernous sinus','straight sinus','transverse sinus','sigmoid sinus','basilar','pericallosal')),
  ('vertebral-spinal',('vertebral','spinal','radicular','medullary')),
  ('carotid-jugular',('carotid','jugular','facial','thyroid','lingual','pharyngeal','occipital','temporal','maxillary','alveolar','meningeal','ophthalm','retinal','angular','diploic','emissary')),
  ('aortic-central',('aorta','aortic arch','brachiocephalic','vena cava','pulmonary trunk')),
  ('upper-limb',('subclavian','axillary','brachial','radial','ulnar','palmar','digital arter','digital vein','metacarpal','circumflex humeral','cephalic vein','basilic','interosseous','suprascapular','circumflex scapular','thoracoacromial','dorsal scapular')),
  ('iliac-pelvic',('iliac','obturator','gluteal','pudendal','vesical','sacral','rectal','uterine','prostatic')),
  ('femoral-popliteal',('femoral','popliteal','genicular','saphenous','perforating arter')),
  ('leg-foot',('tibial','fibular','peroneal','plantar','dorsalis pedis','metatarsal','dorsal venous arch of foot')),
  ('thoracoabdominal-branches',('intercostal','lumbar','thoracic','phrenic','epigastric','musculophrenic','azygos')),
  ('visceral-coronary',('coronary','cardiac','renal','mesenteric','hepatic','gastric','splenic','cystic','gonadal','testicular','pancreatic','duodenal','jejunal','ileal','colic','esophageal'))]:
  if any(t in s for t in terms):return label
 return 'other-named-vessel'
def role(name,bone):
 s=name.lower();b=bone.lower()
 if 'nutrient'in s:return 'named-osseous-supply-review','Explicit nutrient vessel name; bone penetration may be intended, but geometry/entry remains unverified.'
 if ('alveolar'in s or 'dental'in s)and any(x in b for x in ['maxilla','mandib','tooth','teeth']):return 'possible-osseous-supply-review','Dental/alveolar vessel and jaw/teeth pairing; exact canal/entry not established.'
 if any(x in s for x in ['diploic','emissary','meningeal'])and any(x in b for x in ['bone','cranium']):return 'possible-foraminal-or-osseous-review','Name suggests cranial osseous/foraminal supply; crossing still needs opening-specific review.'
 if ('vertebral'in s or 'spinal'in s or 'radicular'in s)and any(x in b for x in ['vertebra','atlas','axis','sacrum']):return 'possible-foraminal-route-review','Spinal/vertebral route against vertebral bone; witness does not prove passage through an opening.'
 if 'internal carotid artery'in s and'temporal bone'in b:return 'possible-foraminal-route-review','Internal carotid / temporal-bone pairing may involve the carotid canal; opening-specific support is unverified.'
 canonical=re.sub(r'^(left|right) ','',s)
 trunks=['abdominal aorta','thoracic aorta','aortic arch','inferior vena cava','superior vena cava','pulmonary trunk','basilar artery','costocervical trunk','thyrocervical trunk','brachiocephalic trunk']
 trunks += [stem+' '+tail for stem in ['common carotid','internal carotid','external carotid','internal jugular','external jugular','subclavian','axillary','brachial','femoral','deep femoral','common femoral','popliteal','common iliac','internal iliac','external iliac','anterior tibial','posterior tibial','fibular','brachiocephalic']for tail in ['artery','vein','veins']]
 trunk=canonical in trunks
 if trunk:return 'major-trunk-crossing-priority','Named major vessel trunk crossing actual bone surface; unintended fit conflict is the working hypothesis, not a validated canal.'
 return 'other-vessel-crossing-review','No opening or osseous target inferred from names. Treat as unresolved strict surface crossing.'
def table(counter,label):return [{label:k,'sampledStrictWitnesses':int(v)}for k,v in counter.most_common()]
allrows={};summaries={};flat=[]
for stage in ['initial_similarity','regional_candidate']:
 path=D/(stage+'-witnesses.npz');inputs.append(path);z=np.load(path,allow_pickle=False)
 assert len(meta)==len(z['relation'])==report['results'][stage]['sampledBloodFaces']
 crossing=np.flatnonzero(z['relation']=='crossing');assert np.array_equal(crossing,z['crossing_sample_indices']);rows=[]
 bone_faces=z['bone_face_ids'];triangles=z['sampled_triangles']
 for si in crossing:
  bf=int(bone_faces[si]);pi=int(bone_owner[bf]);bone=owners[pi];m=meta[si];ins=instances[m['instance']];tier,why=role(m['name'],bone['id'])
  row={**m,'stage':stage,'bonePart':bone['id'],'boneName':bone['name'],'anatomicalRegion':bone['region'],'compoundBoneFace':bf,'originalBonePartFace':int(bone_source_face[bf]),'vesselFamily':family(m['name']),'sourceLineage':ins.get('lineage',[]),'reviewTier':tier,'reviewReason':why,'vesselTriangleCentroid_m':triangles[si].mean(0).tolist()};rows.append(row)
 allrows[stage]=rows;flat.extend(rows)
 summaries[stage]=dict(sampledStrictWitnesses=len(rows),distinctVesselInstances=len({r['instance']for r in rows}),byRegion=table(Counter(r['anatomicalRegion']for r in rows),'anatomicalRegion'),byOriginalBonePart=table(Counter(r['bonePart']for r in rows),'bonePart'),byVesselFamily=table(Counter(r['vesselFamily']for r in rows),'vesselFamily'),byReviewTier=table(Counter(r['reviewTier']for r in rows),'reviewTier'),byVesselInstance=[dict(instance=k,name=next(r['name']for r in rows if r['instance']==k),sampledStrictWitnesses=n)for k,n in Counter(r['instance']for r in rows).most_common()])
before={r['sampleIndex']:r for r in allrows['initial_similarity']};after={r['sampleIndex']:r for r in allrows['regional_candidate']}
groups=defaultdict(lambda:{'initial_similarity':0,'regional_candidate':0,'gainedSampleIndices':[],'lostSampleIndices':[],'persistentSampleIndices':[]})
for stage,rows in allrows.items():
 for r in rows:groups[(r['anatomicalRegion'],r['bonePart'],r['vesselFamily'],r['instance'],r['name'],r['reviewTier'])][stage]+=1
for si in sorted(set(before)|set(after)):
 typ='persistentSampleIndices'if si in before and si in after else'gainedSampleIndices'if si in after else'lostSampleIndices';r=after.get(si,before.get(si));groups[(r['anatomicalRegion'],r['bonePart'],r['vesselFamily'],r['instance'],r['name'],r['reviewTier'])][typ].append(si)
group_rows=[]
for key,x in groups.items():
 region,bone,fam,instance,name,tier=key;group_rows.append(dict(anatomicalRegion=region,bonePart=bone,vesselFamily=fam,instance=instance,name=name,reviewTier=tier,**x,delta=x['regional_candidate']-x['initial_similarity']))
group_rows.sort(key=lambda r:(-r['regional_candidate'],-r['initial_similarity'],r['instance']))
pair_sets=defaultdict(lambda:defaultdict(set))
for stage,rows in allrows.items():
 for r in rows:pair_sets[(r['anatomicalRegion'],r['bonePart'],r['vesselFamily'],r['instance'])][stage].add(r['sampleIndex'])
for r in group_rows:
 key=(r['anatomicalRegion'],r['bonePart'],r['vesselFamily'],r['instance']);a=pair_sets[key]['initial_similarity'];b=pair_sets[key]['regional_candidate'];r['gainedSampleIndices']=sorted(b-a);r['lostSampleIndices']=sorted(a-b);r['persistentSampleIndices']=sorted(a&b)
O.mkdir();(O/'witness-ownership.json').write_text(json.dumps(dict(rows=flat),indent=2)+'\n')
(O/'groups.json').write_text(json.dumps(dict(groups=group_rows),indent=2)+'\n')
out=dict(schema='saved-skeleton-sampled-witness-aggregation/v1',inputHashes={str(p.relative_to(R)):sha(p)for p in inputs},summaries=summaries,pairedSamples=dict(gained=len(set(after)-set(before)),lost=len(set(before)-set(after)),persistent=len(set(after)&set(before))),scope=report['scope'],sampleMetaCount=len(meta),exactFullCrossingCounts=False,noNewSpatialQueries=True,anatomicalAcceptance=False,installed=False,skeletonMoved=False,rawEvidenceChanged=False,limits=['Counts are strict triangle crossing witnesses from a deterministic subset and one proposed bone face each; no full crossing count or prevalence estimate.','Named nutrient/dental/spinal/meningeal/emissary pairings are review categories, never confirmed intended passages.','Zero sampled crossings is not clearance, containment, patency or interface proof.','Suggested fit corrections require additional local review and must preserve all fixed skeleton geometry.'],elapsed_seconds=time.monotonic()-start)
(O/'summary.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(dict(summary=str(O/'summary.json'),stages={k:{q:v[q]for q in ['sampledStrictWitnesses','distinctVesselInstances','byRegion','byReviewTier']}for k,v in summaries.items()},paired=out['pairedSamples'],elapsed=out['elapsed_seconds'])))
