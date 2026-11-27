"""Independent current-source validation of one right IOV middle-profile move.

Run through geometry_tools.runtime. This checker never fits, packs or installs.
Only its JSON report is written. All IDs address the current indexed mesh.
"""
import argparse
from pathlib import Path
import json,time,traceback
import numpy as np

from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.source_port_dependencies import check_source_mouth
from geometry_tools.algorithms.shared_source_normals import transport_shared_source_normals
from geometry_tools.algorithms.triangle_cross_candidates import TriangleCrossIndex

ROOT=Path(__file__).resolve().parents[1]
B=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
SOURCE=B/'ocular-IOV-pterygoid-diagnostic-installed-v2938/circulation-atlas.json'
SOURCE_SHA='d6b77ed0dca78fcb0f097aa25324356a0fe3824d0f6327591414cfce24d77768'
LEDGER=SOURCE.parent/'source-port-dependencies.json'
LEDGER_SHA='60554af73051b44dfbdefed87ef467cb8365364b67499135c1eee61a8b27f429'
PLAN=ROOT/'server/right-iov-profile-plan.json'
PLAN_SHA='5672095eea3f566f573d5c57b812d27b4f5b2ea604ce123bf125b8157d14477b'
PID='BDM:male:human_20_male_cardiovascular_system-right_inferior_ophthalmic_vein_ID'
DP='BDM:male:human_20_male_cardiovascular_system-right_descending_palatine_artery_ID'
FRAME='model:local-ta98:assembled-metres'
RING=np.array([107,109,111,113,115,117],int)
BOUNDARIES=np.array([106,108,110,112,114,116,220,221,222,223,224,225],int)
ANNULI=np.arange(174,198)
BONES=['ZA:Sphenoid bone:Bone-4','ZA:Maxilla.r:Bone-5','ZA:Palatine bone.r:Bone-2','ZA:Zygomatic bone.r:Bone-8']
NATIVE_EYE=['ZA:Sclera.r','ZA:Cornea.r','ZA:Iris.r','ZA:Retina.r','ZA:Lens.r','ZA:Vitreous body.r']

def alias_partition(vertices):
 _,first,inverse=np.unique(vertices,axis=0,return_index=True,return_inverse=True)
 return first[inverse]

def polygon_area(points):
 x=points-points.mean(0)
 return float(np.linalg.norm(np.cross(x,np.roll(x,-1,axis=0)).sum(0))/2)

def diameter(points):
 return float(np.linalg.norm(points[:,None]-points[None,:],axis=2).max())

def codec(points,record):
 assert record.get('geometryInstance') is None,'Unexpected instance frame'
 origin=np.asarray(record.get('geometryOrigin',[0,0,0]),float)
 offset=np.asarray(record.get('sourceOffset',[0,0,0]),float)
 scale=float(record.get('geometryScale',1.))
 return ((np.asarray(points)-origin-offset)/scale).astype('<f4').astype(float)*scale+origin+offset

def common_translation(source,target,record):
 """Certify one common displacement before its declared Float32 rounding."""
 delta=target-source;mean=delta.mean(0)
 origin=np.asarray(record.get('geometryOrigin',[0,0,0]),float)+np.asarray(record.get('sourceOffset',[0,0,0]),float)
 scale=float(record.get('geometryScale',1.))
 words=((target-origin)/scale).astype('<f4')
 previous=np.nextafter(words,np.float32(-np.inf)).astype(float)
 following=np.nextafter(words,np.float32(np.inf)).astype(float)
 low=((previous+words.astype(float))/2)*scale+origin-source
 high=((following+words.astype(float))/2)*scale+origin-source
 lower=low.max(0);upper=high.min(0)
 interval=bool(np.all(lower<=upper));translation=(lower+upper)/2 if interval else mean
 replay=bool(np.array_equal(codec(source+translation,record),target))if interval else False
 roundoff=np.max(np.spacing(np.maximum(np.abs(source),np.abs(target))))*16
 residual=float(np.linalg.norm(delta-translation,axis=1).max())
 exact_float64=bool(np.linalg.norm(delta-mean,axis=1).max()<=roundoff)
 return dict(passed=bool(replay or exact_float64),translation_m=translation.tolist(),
  commonFloat32PreimageExists=interval,commonFloat32ReplayExact=replay,
  commonFloat64Translation=exact_float64,preimageLower_m=lower.tolist(),preimageUpper_m=upper.tolist(),
  maximumStoredPointResidualFromProperTranslation_m=residual,
  numericalFloat64Bound_m=float(roundoff))

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--candidate',required=True,type=Path)
 parser.add_argument('--output',required=True,type=Path)
 parser.add_argument('--wall-seconds',type=float,default=120)
 args=parser.parse_args()
 if not 0<args.wall_seconds<=120:parser.error('At most 120 seconds allowed')
 if args.output.exists():raise FileExistsError('Use a new immutable report path')
 tick=time.monotonic()
 report=dict(schema='right-iov-middle-profile-check/v1',complete=False,passed=False,
  readOnly=True,geometryModified=False,anatomicalAcceptance=False,wholeLumenPatencyEstablished=False,
  sourceCataloguePath=str(SOURCE),sourceCatalogueSha256=SOURCE_SHA,
  sourceMouthLedgerPath=str(LEDGER),sourceMouthLedgerSha256=LEDGER_SHA,
  sourcePlanPath=str(PLAN),sourcePlanSha256=PLAN_SHA,
  candidatePath=str(args.candidate.resolve()),sourceCodePath=str(Path(__file__).resolve()),
  sourceCodeSha256=digest_file(__file__),partId=PID,checks={},stages=[])
 checks=report['checks']
 def deadline():
  if time.monotonic()-tick>max(.1,args.wall_seconds-3):raise TimeoutError('Independent checker deadline; no complete acceptance')
 def require(name,value):
  checks[name]=bool(value)
  if not value:raise ValueError(name)
 try:
  require('currentSourceCatalogueHashBound',digest_file(SOURCE)==SOURCE_SHA)
  require('current542LedgerHashBound',digest_file(LEDGER)==LEDGER_SHA)
  plan=json.loads(PLAN.read_text())
  require('explicitCurrentSourcePlanHashBound',digest_file(PLAN)==PLAN_SHA and plan['sourceCatalogueSha256']==SOURCE_SHA)
  maximum_translation=float(plan['assembly']['maximumTranslationNorm_m'])
  require('sourcePlanMaximumTranslationExactly2mm',maximum_translation==.002)
  a=json.loads(SOURCE.read_text());parts={p['id']:p for p in a['parts']}
  ledger=json.loads(LEDGER.read_text());require('exactSourceInterfaceCount542',len(ledger['interfaces'])==542)
  evidence={}
  for name,expected in [('decoded-payload-proof.json','0aad94e92ae4ed3a7bb2c185b0986e8b0b0251e9be0fcaa8eee502a591d7be99'),('mouth-ledger-extension-review.json',None)]:
   path=SOURCE.parent/name;value=json.loads(path.read_text());sha=digest_file(path)
   require('installedReplay_'+name,bool(value.get('complete')and value.get('passed')and(expected is None or sha==expected)))
   if name=='mouth-ledger-extension-review.json':require('installedReplayCatalogueAndLedgerBound',value['catalogueSha256']==SOURCE_SHA and value['sourcePortLedgerSha256']==LEDGER_SHA)
   evidence[str(path)]=sha
  report['inheritedInstalledInterfaceEvidence']=evidence
  reader=AtlasReader(ROOT,max_open_banks=4)
  def load(pid):return reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=SOURCE,catalogue_sha=SOURCE_SHA)
  s=load(PID);q=Mesh.load(args.candidate);record=parts[PID]
  report.update(candidateSha256=digest_file(args.candidate),sourceGeometrySha256=s.geometry_sha,
   sourceAttributesSha256=s.attributes_sha,candidateGeometrySha256=q.geometry_sha,candidateAttributesSha256=q.attributes_sha)
  require('sourceCoordinateFrameLiteral',q.frame==s.frame)
  require('sourceControlLayoutLiteral',q.vertices.shape==s.vertices.shape)
  require('sourceFacesAndOrderLiteral',np.array_equal(q.faces,s.faces))
  require('sourceIndexedRawIdsLiteral',s.source.get('source_vertex_ids')is not None and np.array_equal(q.source.get('source_vertex_ids'),s.source['source_vertex_ids']))
  require('sourceAttributeFieldsLiteral',set(s.attributes)==set(q.attributes))
  require('allNonNormalAttributeWordsLiteral',all(v.dtype==q.attributes[k].dtype and np.array_equal(v,q.attributes[k])for k,v in s.attributes.items()if k!='normals'))
  aliases=alias_partition(s.vertices);moving=np.flatnonzero(np.isin(aliases,aliases[RING]));fixed=np.setdiff1d(np.arange(len(s.vertices)),moving)
  changed_faces=np.flatnonzero(np.isin(s.faces,moving).any(1));normal_support=np.flatnonzero(np.isin(aliases,aliases[np.unique(s.faces[changed_faces])]))
  report['scope']=dict(middleRingRawIds=RING.tolist(),movingAndExactAliasRawIds=moving.tolist(),fixedBoundaryRawIds=BOUNDARIES.tolist(),
   completeSourceAnnulusFaceIds=ANNULI.tolist(),actualAffectedSourceFaceIds=changed_faces.tolist(),normalTransportSupportRawIds=normal_support.tolist())
  require('onlyAuthorizedProfileAliasClosureMoves',np.array_equal(q.vertices[fixed],s.vertices[fixed]))
  require('bothCompleteBoundaryProfilesLiteral',np.array_equal(q.vertices[BOUNDARIES],s.vertices[BOUNDARIES]))
  require('exactTwoSourceAnnuliAffected',np.array_equal(changed_faces,ANNULI))
  require('noCoordinateAliasesCreatedOrBroken',np.array_equal(alias_partition(q.vertices),aliases))
  require('storedFloat32GridIdempotent',np.array_equal(codec(q.vertices,record),q.vertices))
  field=common_translation(s.vertices[moving],q.vertices[moving],record);report['properProfileTranslation']=field
  require('oneCommonProperProfileTranslation',field['passed'])
  motion=float(np.linalg.norm(q.vertices-s.vertices,axis=1).max());report['maximumSourceControlMotion_m']=motion
  require('authorizedTranslationAtMost2mm',motion<=maximum_translation+field['maximumStoredPointResidualFromProperTranslation_m']+field['numericalFloat64Bound_m'])
  residual=field['maximumStoredPointResidualFromProperTranslation_m']+field['numericalFloat64Bound_m']
  pa,pb=polygon_area(s.vertices[RING]),polygon_area(q.vertices[RING]);da,db=diameter(s.vertices[RING]),diameter(q.vertices[RING])
  centered=s.vertices[RING]-s.vertices[RING].mean(0)
  area_bound=float(sum((np.linalg.norm(x)+np.linalg.norm(y))*residual+2*residual*residual for x,y in zip(centered,np.roll(centered,-1,axis=0))))
  report['profileCalibre']=dict(sourceArea_m2=pa,candidateArea_m2=pb,sourceDiameter_m=da,candidateDiameter_m=db,
   relativeAreaChange=pb/pa-1,relativeDiameterChange=db/da-1,areaErrorBoundFromSingleCodec_m2=area_bound,diameterErrorBoundFromSingleCodec_m=2*residual)
  require('profileAreaPreservedThroughDeclaredCodec',abs(pb-pa)<=area_bound+32*np.finfo(float).eps*pa)
  require('profileDiameterPreservedThroughDeclaredCodec',abs(db-da)<=2*residual+32*np.finfo(float).eps*da)
  old_cross=np.cross(s.triangles[:,1]-s.triangles[:,0],s.triangles[:,2]-s.triangles[:,0]);new_cross=np.cross(q.triangles[:,1]-q.triangles[:,0],q.triangles[:,2]-q.triangles[:,0])
  old_area=np.linalg.norm(old_cross,axis=1);new_area=np.linalg.norm(new_cross,axis=1);valid=old_area>0
  reversed_faces=np.flatnonzero(valid&(np.einsum('ij,ij->i',old_cross,new_cross)<=0));collapsed=np.flatnonzero(valid&(new_area<=0))
  report['sourceFaces']=dict(reversedFaceIds=reversed_faces.tolist(),newCollapsedFaceIds=collapsed.tolist(),minimumAffectedAreaRatio=float((new_area[ANNULI]/old_area[ANNULI]).min()))
  require('noReversedSourceFaces',not len(reversed_faces));require('noNewCollapsedSourceFaces',not len(collapsed))
  normals,normal_proof=transport_shared_source_normals({PID:s},{PID:q.vertices},[[PID]],maximum_facet_rotation_degrees=65,wall_seconds=5)
  report['normalReplay']=normal_proof
  require('freshStoredGridNormalTransportExact',q.attributes['normals'].dtype==s.attributes['normals'].dtype and np.array_equal(q.attributes['normals'],normals[PID]))
  unaffected=np.setdiff1d(np.arange(len(s.vertices)),normal_support)
  require('allUnaffectedNormalWordsLiteral',np.array_equal(q.attributes['normals'][unaffected],s.attributes['normals'][unaffected]))
  before,after=topology(s),topology(q);report['topology']=dict(source=before,candidate=after)
  topology_keys=['triangles','used_vertices','components','open_edges','nonmanifold_edges','inconsistent_winding_edges','invalid_vertex_links','duplicate_faces','zero_area_faces']
  require('sourcePhysicalTopologyLiteral',all(before[k]==after[k]for k in topology_keys))
  report['stages'].append('layout-profile-faces-normal-topology')
  deadline();mouths=[]
  affected=[i for i in ledger['interfaces']if PID in(i['first'],i['second'])]
  require('fiveActualRightIOVInterfacesInventoried',len(affected)==5)
  for item in affected:
   deadline();peer_id=item['second']if item['first']==PID else item['first'];peer=load(peer_id)
   first=s if item['first']==PID else peer;second=peer if item['first']==PID else s
   first_q=q if item['first']==PID else peer;second_q=peer if item['first']==PID else q
   baseline=check_source_mouth(first,second,item['firstRing'],item['secondRing']);new=check_source_mouth(first_q,second_q,item['firstRing'],item['secondRing'])
   ring=item['firstRing']if item['first']==PID else item['secondRing']
   require('interfaceCoordinatesLiteral_'+peer_id,np.array_equal(s.vertices[ring],q.vertices[ring]))
   require('interfaceNormalWordsLiteral_'+peer_id,np.array_equal(s.attributes['normals'][ring],q.attributes['normals'][ring]))
   require('interfaceCompleteAndPreserved_'+peer_id,baseline['passed']and new==baseline)
   mouths.append(dict(first=item['first'],second=item['second'],firstRing=item['firstRing'],secondRing=item['secondRing'],peerGeometrySha256=peer.geometry_sha,source=baseline,candidate=new))
  report['mouthReplay']=dict(sourceInterfaceCount=542,freshLocalInterfaceCount=5,inheritedUnchangedInterfaceCount=537,
   invarianceReason='Every other actor is the literal current catalogue; only one interior profile moves. All five incident owner rims, normal words and source topology are literal and freshly checked.',interfaces=mouths)
  require('all542SourceInterfacesPreserved',len(mouths)==5)
  report['stages'].append('five-local-mouths-and-542-invariance')
  def pairs(mesh,peer=None):
   deadline();data,summary=pair_contacts(mesh,peer or mesh,0.,selfcheck=peer is None,max_pairs=100000,max_tests=2000000)
   if not summary['complete_within_pair_budget']:raise RuntimeError('Actual pair budget exhausted')
   hit=data['pairs'][np.isin(data['relation'],['crossing','coplanar-overlap'])]
   return dict(geometrySha256=mesh.geometry_sha,crossingPairs=sorted(map(list,set(map(tuple,hit.tolist())))),summary=summary)
  old_self,new_self=pairs(s),pairs(q);added_self=sorted(set(map(tuple,new_self['crossingPairs']))-set(map(tuple,old_self['crossingPairs'])))
  report['wholeSelf']=dict(source=old_self,candidate=new_self,addedCrossingPairs=added_self)
  require('noAddedWholeOwnerSelfCrossings',not added_self)
  low=s.bounds[0]-.002;high=s.bounds[1]+.002
  require('candidateWithinDeclared2mmSourceBroadphase',np.all(q.bounds[0]>=low)and np.all(q.bounds[1]<=high))
  selected=[];excluded=[];outside=[]
  for pid,p in parts.items():
   if pid==PID:continue
   fixed_native=pid in BONES or pid in NATIVE_EYE
   blood=pid.startswith('BDM:')and p.get('system')in('arterial','venous')
   if not fixed_native and not blood:continue
   if not fixed_native and(p.get('suppressed')or p.get('referenceOnly')):excluded.append(pid);continue
   bounds=np.asarray(p['bounds'],float)
   if not(np.all(bounds[1]>=low)and np.all(bounds[0]<=high)):
    if pid not in BONES:outside.append(pid);continue
   selected.append(pid)
  require('activeDPAndFourNativeBonesIncluded',DP in selected and set(BONES)<=set(selected))
  contexts=[]
  def cross_against(actor,peer,peer_triangles,active):
   deadline();hits=[];tested=0
   if len(active):
    index=TriangleCrossIndex(peer_triangles[active])
    for batch in index.batches(actor.triangles,batch_size=256,maximum_proposals=2000000):
     deadline();tested+=len(batch);proof=__import__('geometry_tools.distance',fromlist=['triangle_pairs']).triangle_pairs(actor.triangles[batch[:,0]],peer_triangles[active[batch[:,1]]])
     good=np.isin(proof['relation'],['crossing','coplanar-overlap']);hits.extend(np.c_[batch[good,0],active[batch[good,1]]].tolist())
   return dict(ownerGeometrySha256=actor.geometry_sha,peerGeometrySha256=peer.geometry_sha,peerAttributesSha256=peer.attributes_sha,crossingPairs=sorted(map(list,set(map(tuple,hits)))),actualTrianglePairsTested=tested)
  for pid in selected:
   deadline();peer=load(pid);tri=peer.triangles
   active=np.flatnonzero(np.all(tri.max(1)>=low,axis=1)&np.all(tri.min(1)<=high,axis=1)&(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)>0))
   old=cross_against(s,peer,tri,active);new=cross_against(q,peer,tri,active);old_set=set(map(tuple,old['crossingPairs']));new_set=set(map(tuple,new['crossingPairs']))
   contexts.append(dict(peerId=pid,physicalRole='fixed-native-skeleton'if pid in BONES else'fixed-native-eye'if pid in NATIVE_EYE else'active-BioDigital-vessel',defaultHidden=pid in a.get('defaultHidden',[])or bool(parts[pid].get('defaultHidden')),
    allActualOwnerFacesChecked=True,peerTotalFaces=len(tri),peerFacesWithinConservativeSource2mmBounds=len(active),
    sourceAndCandidate=[dict(state='source',**old),dict(state='candidate',**new)],addedCrossingPairs=sorted(new_set-old_set),removedCrossingPairs=sorted(old_set-new_set)))
  report['fixedContext']=dict(allCurrentMetadataBoundsConsidered=len(parts),sourceCropBounds_m=[low.tolist(),high.tolist()],
   scope='Four fixed native bones, native eye peers, and every active BioDigital blood mesh whose bounds overlap the complete source owner plus2mm. Non-BDM vascular guidance and explicitly suppressed/reference alternatives are not active physical obstacles. DefaultHidden alone does not exclude a source.',
   excludedSuppressedOrReferenceBDM=excluded,outsideConservativeBounds=outside,rows=contexts)
  dp=next(row for row in contexts if row['peerId']==DP)
  require('currentDP24PairBaselineConfirmed',len(dp['sourceAndCandidate'][0]['crossingPairs'])==24)
  checks['descendingPalatineProperCrossingsCleared']=not dp['sourceAndCandidate'][1]['crossingPairs']
  checks['noAddedNativeBoneCrossings']=not any(row['addedCrossingPairs']for row in contexts if row['peerId']in BONES)
  checks['noAddedNativeEyeCrossings']=not any(row['addedCrossingPairs']for row in contexts if row['peerId']in NATIVE_EYE)
  checks['noAddedActiveBDMVesselCrossings']=not any(row['addedCrossingPairs']for row in contexts if row['physicalRole']=='active-BioDigital-vessel')
  require('sourceCatalogueStillLiteral',digest_file(SOURCE)==SOURCE_SHA and digest_file(LEDGER)==LEDGER_SHA)
  report['stages'].append('whole-owner-self-and-current-context')
  report.update(complete=True,passed=all(checks.values()),contextAnatomicallyAccepted=False,
   inheritedDefectsRemain=[dict(peerId=row['peerId'],candidateCrossingPairs=row['sourceAndCandidate'][1]['crossingPairs'])for row in contexts if row['sourceAndCandidate'][1]['crossingPairs']])
 except Exception as error:
  report['failure']=dict(type=type(error).__name__,message=str(error),traceback=traceback.format_exc())
 finally:
  report['elapsedSeconds']=time.monotonic()-tick
  atomic_json(args.output,report)
  print(json.dumps(dict(path=str(args.output.resolve()),sha256=digest_file(args.output),complete=report['complete'],passed=report['passed'],failure=report.get('failure',{}).get('message'),elapsedSeconds=report['elapsedSeconds'])))
 return 0 if report['passed']else 1

if __name__=='__main__':raise SystemExit(main())
