"""Independent, read-only verification of the declared Angular/IOV rim overlay."""
import hashlib,json,time
from collections import Counter
from pathlib import Path
import numpy as np
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.distance import triangle_pairs
from geometry_tools.algorithms.triangle_cross_candidates import TriangleCrossIndex
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.shared_source_normals import transport_shared_source_normals
from geometry_tools.source_port_dependencies import check_source_mouth

ROOT=Path(__file__).resolve().parents[1];B=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1';D=B/'right-angular-iov-shared-rim-v2951-agent'
REPORT=D/'report.json';PLAN=ROOT/'server/right-angular-iov-junction-plan.json';OUT=D/'independent-check.json';CP=B/'ocular-right-iov-profile-preview-v2948-install/circulation-atlas.json'
PREFIX='BDM:male:human_20_male_cardiovascular_system-';ANG=PREFIX+'right_angular_vein_ID';IOV=PREFIX+'right_inferior_ophthalmic_vein_ID';OWNERS=[ANG,IOV]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
identity=lambda x:hashlib.sha256(json.dumps(x,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def aliases(v):
 _,first,inv=np.unique(v,axis=0,return_index=True,return_inverse=True);return first[inv]
def loops(m):return boundary_loops(aliases(m.vertices)[m.faces])
def codec(v,o):return (v-o).astype('<f4').astype(float)+o
def area(m):
 t=m.triangles;return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
def whole_pairs(m,p=None):
 contacts,proof=pair_contacts(m,p if p is not None else m,0.,selfcheck=p is None,max_pairs=200000,max_tests=3000000)
 assert proof['complete_within_pair_budget']and not proof['truncated']
 pairs=sorted(set(map(tuple,contacts['pairs'][np.isin(contacts['relation'],['crossing','coplanar-overlap'])].tolist())))
 return dict(crossingPairs=pairs,proof=proof)

def main():
 tick=time.monotonic();assert not OUT.exists();r=json.loads(REPORT.read_text());plan=json.loads(PLAN.read_text());catsha=sha(CP)
 assert catsha==r['sourceCatalogueSha256']==plan['sourceCatalogueSha256']=='948e0b6dcb45985f115bd795ae58a47ea0c860e98ddc6e01389cfdfa993f4b24'
 assert sha(PLAN)==r['sourcePlanSha256']=='14cd0068495fcfb0cca37041f6baf2c637d93a08bd1c33f76c6c74251a393099'
 assert sha(r['sourceCodePath'])==r['sourceCodeSha256'];assert sha(r['normalPolicyPath'])==r['normalPolicySha256']
 atlas=json.loads(CP.read_text());parts={p['id']:p for p in atlas['parts']};assert len(parts)==8779;reader=AtlasReader(max_open_banks=4);cache={}
 def get(pid):
  if pid not in cache:cache[pid]=reader.read_part(atlas,parts[pid],frame='model:local-ta98:assembled-metres',catalogue_path=CP,catalogue_sha=catsha)
  return cache[pid]
 source={pid:get(pid)for pid in OWNERS};rows={x['partId']:x for x in r['rows']};assert set(rows)==set(OWNERS);prows={x['partId']:x for x in plan['rows']};target={pid:Mesh.load(rows[pid]['path'])for pid in OWNERS}
 ledgerPath=Path(r['sourceLedgerPath']);assert sha(ledgerPath)==r['sourceLedgerSha256']==plan['sourceLedgerSha256'];ledger=json.loads(ledgerPath.read_text());assert len(ledger['interfaces'])==542
 link=plan['assembly']['canonicalInterface'];ar=np.array(link['firstRing']);ir=np.array(link['secondRing']);assert link['first']==ANG and link['second']==IOV
 origin=np.asarray(parts[IOV]['geometryOrigin']);sa,si=source[ANG],source[IOV];qa,qi=target[ANG],target[IOV];ideal=sa.vertices.copy();ideal[ar]=si.vertices[ir]
 expected=codec(ideal,origin);np.testing.assert_array_equal(qa.vertices,expected);np.testing.assert_array_equal(qi.vertices,si.vertices);np.testing.assert_array_equal(qa.vertices[ar],qi.vertices[ir]);np.testing.assert_array_equal(qi.vertices[[107,109,111,113,115,117]],si.vertices[[107,109,111,113,115,117]])
 local=ideal-origin;f=local.astype('<f4');fv=f.astype(float);spacing=np.maximum(abs(np.nextafter(f,np.float32(np.inf)).astype(float)-fv),abs(fv-np.nextafter(f,np.float32(-np.inf)).astype(float)));bounds=spacing/2+abs(np.spacing(local))/2+abs(np.spacing(expected))/2;errors=abs(qa.vertices-ideal);assert np.all(errors<=bounds)
 codecproof={'componentError_m':errors.tolist(),'componentBound_m':bounds.tolist(),'maximumActualIdealToStoredError_m':float(np.linalg.norm(errors,axis=1).max()),'maximumAnalyticalEuclideanBound_m':float(np.linalg.norm(bounds,axis=1).max()),'maximumRebaseOnlyDisplacement_m':float(np.linalg.norm(codec(sa.vertices,origin)-sa.vertices,axis=1).max()),'maximumCanonicalPhysicalCorrection_m':float(np.linalg.norm(sa.vertices[ar]-si.vertices[ir],axis=1).max()),'allComponentBoundsPassed':True}
 checks={};arrayrows=[];words={};normalproof={}
 for pid in OWNERS:
  s,q=source[pid],target[pid];row=rows[pid];assert sha(row['path'])==row['sha256'];assert sha(row['sourcePath'])==row['sourceSha256'];snapshot=Mesh.load(row['sourcePath']);np.testing.assert_array_equal(snapshot.vertices,s.vertices);np.testing.assert_array_equal(snapshot.faces,s.faces);assert snapshot.attributes_sha==s.attributes_sha;assert row['sourcePartRecord']==parts[pid]
  assert s.frame==q.frame;np.testing.assert_array_equal(s.faces,q.faces);assert s.source['source_vertex_ids']==q.source['source_vertex_ids'];assert set(s.attributes)==set(q.attributes)
  assert q.attributes['normals'].dtype==s.attributes['normals'].dtype
  for key,val in s.attributes.items():
   if key!='normals':np.testing.assert_array_equal(val,q.attributes[key]);assert val.dtype==q.attributes[key].dtype
  np.testing.assert_array_equal(aliases(s.vertices),aliases(q.vertices));assert sorted(map(lambda x:sorted(x.tolist()),loops(s)))==sorted(map(lambda x:sorted(x.tolist()),loops(q)))
  assert np.isfinite(q.vertices).all();np.testing.assert_array_equal(codec(q.vertices,origin),q.vertices)
  independent,npf=transport_shared_source_normals({pid:s},{pid:q.vertices},[[pid]],maximum_facet_rotation_degrees=65,wall_seconds=5);words[pid]=independent[pid].copy();normalproof[pid]=npf
  if pid==ANG:fixed=np.unique(np.concatenate([x for x in loops(s)if set(x)!=set(ar)]))
  else:fixed=np.setdiff1d(np.arange(len(s.vertices)),ir)
  assert set(fixed)==set(prows[pid]['fixedLiteralNormalRawIds']);words[pid][fixed]=s.attributes['normals'][fixed];np.testing.assert_array_equal(q.attributes['normals'][fixed],s.attributes['normals'][fixed])
  st,qt=topology(s),topology(q)
  for key in ['triangles','used_vertices','components','open_edges','nonmanifold_edges','inconsistent_winding_edges','invalid_vertex_links','duplicate_faces','zero_area_faces']:assert st[key]==qt[key]
  before=whole_pairs(s);after=before if np.array_equal(s.vertices,q.vertices)else whole_pairs(q);added=sorted(set(after['crossingPairs'])-set(before['crossingPairs']));assert not added
  arrayrows.append({'partId':pid,'sourceGeometrySha256':s.geometry_sha,'candidateGeometrySha256':q.geometry_sha,'sourceAttributesSha256':s.attributes_sha,'candidateAttributesSha256':q.attributes_sha,'sourceBankPath':s.source['bank_path'],'sourceBankSha256':s.source['bank_sha256'],'sourceCompactToBankVertexIndices':s.source['source_vertex_ids'],'sourcePath':row['sourcePath'],'sourceSha256':sha(row['sourcePath']),'path':row['path'],'sha256':sha(row['path']),'facesAndNonNormalAttributesLiteral':True,'fixedNormalIndices':fixed.tolist(),'sourceTopology':st,'candidateTopology':qt,'wholeSelf':{'source':before,'candidate':after,'addedCrossingPairs':added}})
 areas={pid:area(target[pid])for pid in OWNERS};finish=[]
 for x,y in zip(ar,ir):
  vector=np.zeros(3)
  for pid,raw in [(ANG,x),(IOV,y)]:
   s=source[pid];eq=np.flatnonzero(np.all(s.vertices==s.vertices[raw],axis=1));incident=np.flatnonzero(np.any(np.isin(s.faces,eq),axis=1));inc=incident[np.isfinite(areas[pid][incident])&(areas[pid][incident]>0)];w=float(areas[pid][inc].sum());sign=int(s.attributes.get('normalSigns',np.ones(len(s.vertices),int))[raw]);effective=words[pid][raw].astype(np.int32)*sign;vector+=w*effective/np.linalg.norm(effective)
  canonical=np.rint(np.clip(vector/np.linalg.norm(vector),-1,1)*32767).astype('<i2')
  for pid,raw in [(ANG,x),(IOV,y)]:
   s=source[pid];eq=np.flatnonzero(np.all(s.vertices==s.vertices[raw],axis=1));sign=s.attributes.get('normalSigns',np.ones(len(s.vertices),int))[eq];words[pid][eq]=(canonical[None].astype(int)*sign[:,None]).astype('<i2')
  finish.append(canonical.tolist())
 for pid in OWNERS:np.testing.assert_array_equal(words[pid],target[pid].attributes['normals'])
 incident=[e for e in ledger['interfaces']if e['first']in OWNERS or e['second']in OWNERS];unchanged=[e for e in ledger['interfaces']if e not in incident];assert len(incident)==11 and len(unchanged)==531;mouthrows=[]
 for e in incident:
  first,second=e['first'],e['second'];a,b=get(first),get(second);aa,bb=target.get(first,a),target.get(second,b);x,y=e['firstRing'],e['secondRing'];sg=np.linalg.norm(a.vertices[x]-b.vertices[y],axis=1);qg=np.linalg.norm(aa.vertices[x]-bb.vertices[y],axis=1);motion=np.linalg.norm(aa.vertices[x]-a.vertices[x],axis=1)+np.linalg.norm(bb.vertices[y]-b.vertices[y],axis=1);assert np.all(qg<=sg+motion+2e-16)
  proof=check_source_mouth(aa,bb,x,y);assert proof['passed']
  if {first,second}==set(OWNERS):assert np.max(qg)==0
  mouthrows.append({'first':first,'second':second,'firstRing':x,'secondRing':y,'sourceGap_m':sg.tolist(),'candidateGap_m':qg.tolist(),'endpointMotionBound_m':motion.tolist(),'candidateMouthProof':proof})
 identities={pid:parts[pid]for e in unchanged for pid in [e['first'],e['second']]};assert not set(identities)&set(OWNERS);unaffectedIdentity=identity({'interfaces':unchanged,'ownerRecords':identities,'sourceCatalogueSha256':catsha})
 old12=[]
 for e in ledger['interfaces']:
  names=[e['first'].removeprefix(PREFIX),e['second'].removeprefix(PREFIX)];vortex=any('vorticose_veins_ID'in n for n in names)and any('ophthalmic_vein_ID'in n for n in names);cav=any('cavernous_sinus_ID'in n for n in names)and any('ophthalmic_vein_ID'in n for n in names)
  if vortex or cav:
   aa=target.get(e['first']);bb=target.get(e['second']);aa=get(e['first'])if aa is None else aa;bb=get(e['second'])if bb is None else bb;gap=np.linalg.norm(aa.vertices[e['firstRing']]-bb.vertices[e['secondRing']],axis=1);assert gap.max()==0;old12.append({'first':e['first'],'second':e['second'],'maximumGap_m':0.})
 assert len(old12)==12
 sourceLocal=whole_pairs(sa,si);targetLocal=whole_pairs(qa,qi);assert len(sourceLocal['crossingPairs'])==7 and not targetLocal['crossingPairs'];joined=[]
 for stage,pair in [('source',[sa,si]),('candidate',[qa,qi])]:
  merged=Mesh(np.vstack([pair[0].vertices,pair[1].vertices]),np.vstack([pair[0].faces,pair[1].faces+len(pair[0].vertices)]),pair[0].frame);joined.append({'state':stage,'topology':topology(merged)})
 assert len(joined[0]['topology']['components'])==2 and len(joined[1]['topology']['components'])==1;assert joined[0]['topology']['open_edges']==88 and joined[1]['topology']['open_edges']==76
 for key in ['nonmanifold_edges','inconsistent_winding_edges','duplicate_faces','zero_area_faces']:assert joined[1]['topology'][key]==0
 result={'schema':'independent-right-angular-iov-junction-check/v1','complete':False,'passed':False,'readOnly':True,'geometryModified':False,'anatomicalAcceptance':False,'wholeLumenValidated':False,'sourceCataloguePath':str(CP),'sourceCatalogueSha256':catsha,'sourceReportPath':str(REPORT),'sourceReportSha256':sha(REPORT),'sourcePlanPath':str(PLAN),'sourcePlanSha256':sha(PLAN),'sourceLedgerPath':str(ledgerPath),'sourceLedgerSha256':sha(ledgerPath),'codecProof':codecproof,'rows':arrayrows,'independentSourceNormalTransports':normalproof,'sharedEffectiveNormalConsensusWords':finish,'incidentMouths':mouthrows,'unaffected531InterfaceIdentity':{'count':531,'sourceIdentitySha256':unaffectedIdentity,'candidateOverlayIdentitySha256':unaffectedIdentity,'candidateOverlayOwnerIds':OWNERS,'noOverlayOwnersInTheseInterfaces':True,'sourceCatalogueUnmodified':True},'allTwelvePriorExactMouths':old12,'AngularIOVWholePair':{'source':sourceLocal,'candidate':targetLocal},'joinedTopology':joined,'physicalContext':{'allCurrentPartBoundsConsidered':8779,'wholeBothOwners':True,'defaultHiddenAloneExcluded':False,'rows':[]},'sourceCodePath':str(Path(__file__).resolve()),'sourceCodeSha256':sha(__file__),'elapsedSeconds':time.monotonic()-tick}
 atomic_json(OUT,result)
 exclusions=Counter();nearExcluded=[];selected=[];hidden=set(atlas.get('defaultHidden',[]));both=np.concatenate([s.triangles for s in source.values()]+[q.triangles for q in target.values()]);lo=both.min((0,1));hi=both.max((0,1));proposals=0
 for part in atlas['parts']:
  pid=part['id'];system=part.get('system');box=np.asarray(part['bounds']);near=bool(np.all(box[1]>=lo)&np.all(box[0]<=hi));reason=None
  if pid in OWNERS:reason='two-candidate-owners-checked-independently'
  elif part.get('suppressed'):reason='suppressed'
  elif part.get('referenceOnly'):reason='reference-only'
  elif system in ['arterial','venous','nervous','lymphatic']and not pid.startswith('BDM:'):reason='retained-ZA-vascular-neural-or-lymph-guidance'
  elif ':female:'in pid.lower()or pid.startswith('TA98FEMALE:')or any('female'in str(part.get(k,'')).lower()for k in ['reviewSource','sourceSex','sex']):reason='foreign-female-alternative'
  elif not near:reason='outside-source-and-candidate-whole-owner-bounds'
  if reason:
   exclusions[reason]+=1
   if near:nearExcluded.append({'partId':pid,'reason':reason,'defaultHidden':pid in hidden or bool(part.get('defaultHidden'))})
  else:selected.append(part)
 def triangle_cross(owner,peer):
  nonlocal proposals
  ot=owner.triangles;pt=peer.triangles;low=ot.min((0,1));high=ot.max((0,1));valid=np.linalg.norm(np.cross(pt[:,1]-pt[:,0],pt[:,2]-pt[:,0]),axis=1)>0;ids=np.flatnonzero(valid&np.all(pt.max(1)>=low,axis=1)&np.all(pt.min(1)<=high,axis=1));hits=[];touch=[];tested=0
  if len(ids):
   for batch in TriangleCrossIndex(pt[ids]).batches(ot,batch_size=128,maximum_proposals=3000000):
    assert time.monotonic()-tick<110,'Finite110-second independent audit';tested+=len(batch);proposals+=len(batch);proof=triangle_pairs(ot[batch[:,0]],pt[ids[batch[:,1]]]);mapped=np.c_[batch[:,0],ids[batch[:,1]]];hits.extend(mapped[np.isin(proof['relation'],['crossing','coplanar-overlap'])].tolist());touch.extend(mapped[proof['relation']=='touching'].tolist())
  return {'crossingPairs':sorted(set(map(tuple,hits))),'touchingPairs':sorted(set(map(tuple,touch))),'actualTrianglePairsTested':tested,'peerActualOverlappingFaces':len(ids)}
 for part in selected:
  assert time.monotonic()-tick<110,'Finite110-second independent audit';pid=part['id'];peer=get(pid);record={'peerId':pid,'system':part.get('system'),'defaultHidden':pid in hidden or bool(part.get('defaultHidden')),'peerGeometrySha256':peer.geometry_sha,'peerAttributesSha256':peer.attributes_sha,'peerSourceBankPath':peer.source['bank_path'],'peerSourceBankSha256':peer.source['bank_sha256'],'peerRecordSha256':identity(part),'actualWholePeerFaces':len(peer.faces),'owners':[]}
  for owner in OWNERS:
   ss=triangle_cross(source[owner],peer);qq=ss if np.array_equal(source[owner].vertices,target[owner].vertices)else triangle_cross(target[owner],peer);before=set(ss['crossingPairs']);after=set(qq['crossingPairs']);added=sorted(after-before);removed=sorted(before-after);rr={'partId':owner,'allActualOwnerFacesChecked':True,'sourceGeometrySha256':source[owner].geometry_sha,'candidateGeometrySha256':target[owner].geometry_sha,'source':ss,'candidate':qq,'addedCrossingPairs':added,'removedCrossingPairs':removed}
   if added:rr['addedWitnessTriangles']=[{'actorFaceId':f,'peerFaceId':g,'sourceActorTriangle_m':source[owner].triangles[f].tolist(),'candidateActorTriangle_m':target[owner].triangles[f].tolist(),'peerTriangle_m':peer.triangles[g].tolist()}for f,g in added]
   record['owners'].append(rr)
  result['physicalContext']['rows'].append(record)
  if any(o['addedCrossingPairs']for o in record['owners']):result['firstFailurePeer']=pid;break
  if len(result['physicalContext']['rows'])%8==0:result['elapsedSeconds']=time.monotonic()-tick;atomic_json(OUT,result,overwrite=True)
  if pid not in {p for e in incident for p in [e['first'],e['second']]}:cache.pop(pid,None)
 context=result['physicalContext'];context.update({'selectedWholeOwnerBoundsPeerCount':len(selected),'actualPeersChecked':len(context['rows']),'exclusionCounts':dict(exclusions),'nearExcludedPeers':nearExcluded,'actualTrianglePairsTested':proposals,'sourceWholeOwnerBounds_m':[lo.tolist(),hi.tolist()],'totalAddedCrossingPairs':sum(len(o['addedCrossingPairs'])for row in context['rows']for o in row['owners']),'totalRemovedCrossingPairs':sum(len(o['removedCrossingPairs'])for row in context['rows']for o in row['owners'])});assert sha(CP)==catsha;result['complete']=len(context['rows'])==len(selected);result['passed']=result['complete']and not context['totalAddedCrossingPairs'];result['elapsedSeconds']=time.monotonic()-tick;atomic_json(OUT,result,overwrite=True)
 print(json.dumps({'path':str(OUT),'sha256':sha(OUT),'complete':result['complete'],'passed':result['passed'],'elapsedSeconds':result['elapsedSeconds'],'checkedPeers':context['actualPeersChecked'],'selectedPeers':len(selected),'addedPairs':context['totalAddedCrossingPairs'],'removedPairs':context['totalRemovedCrossingPairs'],'selfAdded':0,'incidentMouths':len(mouthrows),'unchangedIdentities':531,'rebaseMax_m':codecproof['maximumRebaseOnlyDisplacement_m'],'contextNonzero':[{'peer':row['peerId'],'owners':[{'id':o['partId'],'source':len(o['source']['crossingPairs']),'candidate':len(o['candidate']['crossingPairs']),'added':len(o['addedCrossingPairs'])}for o in row['owners']if o['source']['crossingPairs']or o['candidate']['crossingPairs']]}for row in context['rows']if any(o['source']['crossingPairs']or o['candidate']['crossingPairs']for o in row['owners'])]}))

if __name__=='__main__':main()
