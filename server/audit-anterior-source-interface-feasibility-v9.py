"""Stage 1 only: objective source/recipient correspondence, without fitting."""
from pathlib import Path
import sys,time,json,sqlite3
from collections import defaultdict
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.review import prepare

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';BP=ROOT/'public/models/atlas.json';AUTH=ROOT/'.local-models/ta98-authority.sqlite'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1';OUT=BASE/'anterior-source-interface-feasibility-v9'
FRAME='model:local-ta98:assembled-metres';ROOTT=['A12.2.06.001','A12.2.07.022','A12.2.07.029']

def serial(q):
    if isinstance(q,np.ndarray):return q.tolist()
    if isinstance(q,np.generic):return q.item()
    raise TypeError(type(q).__name__)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);raw=CAT.read_bytes();current_sha=__import__('hashlib').sha256(raw).hexdigest();a=json.loads(raw);b=json.loads(BP.read_text());ap={p['id']:p for p in a['parts']};bp={p['id']:p for p in b['parts']};reader=AtlasReader();meshes={};trees={};sources={};entries={};source_record_bindings={}
    def get(pid):
        if pid not in meshes:
            m=reader.read_part(a,ap[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=current_sha);meshes[pid]=m;source_record_bindings[pid]=identity(ap[pid]);sources[pid]=dict(record=ap[pid],bankSha256=m.source['bank_sha256'],geometrySha256=m.geometry_sha,attributeSha256=m.attributes_sha,unchangedFixedNativeSource=True)
        return meshes[pid]
    def tree(pid):
        if pid not in trees:trees[pid]=BVH(meshes[pid])
        return trees[pid]
    def nearest(first,second):
        m,n=meshes[first],meshes[second];r=tree(first).closest_mesh(tree(second));r.update(firstPartId=first,secondPartId=second,firstTriangle_m=m.triangles[r['first_face_id']],secondTriangle_m=n.triangles[r['second_face_id']]);return r
    def rim_ref(pid,loop,index):
        m=get(pid);v=m.vertices[loop];f=np.c_[np.arange(len(v)),np.roll(np.arange(len(v)),-1),np.arange(len(v))];rid=pid+' actual boundary '+str(index);meshes[rid]=Mesh(v,f,FRAME,{},dict(diagnostic_line_queries_only=True,notATissueMesh=True));return rid,dict(partId=pid,actualSourceBoundaryLoop=index,sourceVertexIds=loop.tolist(),worldPoints_m=v.tolist(),centroid_m=v.mean(0).tolist(),sourceGeometrySha256=m.geometry_sha,pairedWallAnnulusEstablished=False,notAnEmittedDiskOrCap=True)
    db=sqlite3.connect('file:'+str(AUTH)+'?mode=ro',uri=True);db.row_factory=sqlite3.Row;termrows={r['id']:dict(r)for r in db.execute('select * from ta98')};children=defaultdict(list)
    for r in termrows.values():children[r['parent_id']].append(r['id'])
    authority=[];seen=set()
    for root in ROOTT:
        stack=[root]
        while stack:
            t=stack.pop()
            if t in seen:continue
            seen.add(t);authority.append(termrows[t]);stack.extend(children[t])
    census=json.loads((ROOT/'reports/ta98-blood-circulation-inventory-installed-v11.json').read_text());censusby={r['id']:r for r in census['terms']};sourceconcept={c['id']:c for c in b['concepts']}
    inventory=[]
    for t in authority:
        c=sourceconcept.get('FMA'+str(t['fma_id']));row=censusby.get(t['id'],{})
        inventory.append(dict(authority=t,existingCurrentDirectParts=row.get('resolved_direct_part_ids',[]),existingCurrentDescendantOnlyParts=row.get('resolved_descendant_only_part_ids',[]),originalBp3dConcept=c,plannedRepresentation='No construction in stage 1. Classify each authority term and establish actual source/interface identity before any wall or region partition.',missingRepresentations=row.get('missing_requirements',[])))
    atomic_json(OUT/'complete-authority-source-inventory.json',dict(sourceCatalogueSha256=current_sha,authoritySha256=digest_file(AUTH),roots=ROOTT,terms=inventory,allDescendantsInventoried=True,noNewRepresentationsClaimed=True))
    h=a['heightAdjustment'];scale=float(h['scale']);translation=np.array(h['translationMeters']);donors=['FJ1654','FJ1654M','FJ1655','FJ1682','FJ1682M']
    for pid in donors:
        x=reader.read_part(b,bp[pid],frame='source:BodyParts3D-4.0:packed-native-metres-y-up',catalogue_path=BP);v=(x.vertices*scale+translation).astype('<f4').astype(float);key='BP3D-height-replay:'+pid
        meshes[key]=Mesh(v,x.faces,FRAME,x.attributes,dict(x.source,existing_recorded_height_replay=h,no_new_registration=True));sources[key]=dict(originalRecord=bp[pid],sourceGeometrySha256=x.geometry_sha,sourceBankSha256=x.source['bank_sha256'],replayedGeometrySha256=meshes[key].geometry_sha,recordedHeightTransformOnly=True,noAnatomicalRegistrationEstablished=True)
    native=['ZA:Anterior cerebral artery.l','ZA:Anterior cerebral artery.r','ZA:Anterior communicating artery','ZA:Internal carotid artery.l','ZA:Internal carotid artery.r','ZA:Pericallosal artery.l','ZA:Pericallosal artery.r','ZA:Callosomarginal artery.l','ZA:Callosomarginal artery.r','ZA:Optic chiasm.l','ZA:Optic chiasm.r','ZA:Sphenoid bone:Bone-4']
    for pid in native:get(pid)
    interfaces=[]
    for side,dp,di in [('l','FJ1654M','FJ1682M'),('r','FJ1654','FJ1682')]:
        aca='ZA:Anterior cerebral artery.'+side;ica='ZA:Internal carotid artery.'+side;aloops=boundary_loops(get(aca).faces);iloops=boundary_loops(get(ica).faces);ar=[];ir=[]
        assert len(aloops)==2 and len(iloops)==2
        for i,q in enumerate(aloops):
            rid,info=rim_ref(aca,q,i);info['toFixedIca']=nearest(rid,ica);info['toReplayedBpAca']=nearest(rid,'BP3D-height-replay:'+dp);ar.append((rid,info))
        prox=int(np.argmin([x[1]['toFixedIca']['distance_m']for x in ar]));dist=1-prox
        for i,q in enumerate(iloops):
            rid,info=rim_ref(ica,q,i);info['toFixedAca']=nearest(rid,aca);info['toReplayedBpAca']=nearest(rid,'BP3D-height-replay:'+dp);info['toReplayedBpIca']=nearest(rid,'BP3D-height-replay:'+di);ir.append((rid,info))
        terminal=int(np.argmin([x[1]['toFixedAca']['distance_m']for x in ir]));continuations=[]
        for target in ['ZA:Pericallosal artery.'+side,'ZA:Callosomarginal artery.'+side]:
            loops=boundary_loops(get(target).faces);candidate=[]
            for j,q in enumerate(loops):
                rid,info=rim_ref(target,q,j);info['toActualDistalAcaRim']=nearest(rid,ar[dist][0]);info['toReplayedBpAca']=nearest(rid,'BP3D-height-replay:'+dp);candidate.append(info)
            chosen=min(candidate,key=lambda q:q['toActualDistalAcaRim']['distance_m'])if candidate else None
            continuations.append(dict(partId=target,actualOpenBoundaryLoops=len(loops),closestActualRecipientRim=chosen,distalAcaRimToRecipientSurface=nearest(ar[dist][0],target),replayedBpAcaToWholeRecipientSurface=nearest('BP3D-height-replay:'+dp,target),interfaceEstablished=False))
        interfaces.append(dict(side=side,nativeIcaRootRim=ir[terminal][1],nativeIcaOtherProtectedRim=ir[1-terminal][1],nativeAcaProximalRim=ar[prox][1],nativeAcaDistalRim=ar[dist][1],nativeDistalContinuations=continuations,sourceDonorClosedCastHasNoLiteralOpenMouthOrPairedAnnulus=True,fitEstablished=False))
    brain_ids=[p['id']for p in a['parts']if p['id'].startswith('ZA:')and p['id'].rsplit(':',1)[-1]in ['Frontal lobe','Temporal lobe','Parietal lobe','Occipital lobe']];vv=[];ff=[];off=0
    for pid in brain_ids:
        m=get(pid);vv.append(m.vertices);ff.append(m.faces+off);off+=len(m.vertices)
    meshes['Fixed-ZA-cortex-context']=Mesh(np.vstack(vv),np.vstack(ff),FRAME,{},dict(literal_fixed_source_part_ids=brain_ids));context=['ZA:Optic chiasm.l','ZA:Optic chiasm.r','ZA:Sphenoid bone:Bone-4','Fixed-ZA-cortex-context'];context_results=[]
    for donor in ['BP3D-height-replay:FJ1654','BP3D-height-replay:FJ1654M','BP3D-height-replay:FJ1655']:
        for target in context:
            d,proof=pair_contacts(meshes[donor],meshes[target],0.,first_bvh=tree(donor),second_bvh=tree(target),max_pairs=250000,max_tests=2000000);assert proof['complete_within_pair_budget']
            context_results.append(dict(donor=donor,fixedContext=target,closest=nearest(donor,target),actualBoundaryContactSummary=proof,actualCrossingFacePairs=d['pairs'][d['relation']=='crossing'],volumeOccupancyOrAnatomicalClearanceCertified=False))
    contract=dict(schema='anterior-communicating-whole-assembly-stage1-contract/v9',wholeAssembly='Bilateral ACA/ACom communicating assembly and explicit ICA root/distal ACA interfaces',authorityRoots=ROOTT,initializationOnly=dict(originalBp3dCatalogueSha256=digest_file(BP),existingHeightScale=scale,existingHeightTranslation_m=translation.tolist()),donorSourceIds=donors,allOtherAnatomyFixed=True,fixedNativeParts=native+brain_ids,adjoiningInterfaces=interfaces,constructionAuthorizedAtThisStage=False,pairedWallConstructionBlockedUntilActualRootRecipientFit=True,noSourceCapsAreValidOpenBloodMouths=True,noA1A2SubdivisionBeforeActualJunction=True,forearmInstallationIsIndependentAndMustBePreservedOnLaterRebase=True,permittedCurrentAction='read-only correspondence measurements and diagnostic source snapshots',futureAcceptanceRequires=['Explicit source/recipient root correspondence with real paired interfaces','Fixed native brain/skeleton/unrelated artery replay','Whole joined arterial wall and blood-space topology, continuity and calibre','Native rendered whole surface and longitudinal junction cuts'])
    atomic_json(OUT/'assembly-contract.json',json.loads(json.dumps(contract,default=serial)))
    def reviewentry(pid,color,lower=None):
        m=meshes[pid];ids=np.arange(len(m.faces))if lower is None else np.flatnonzero(np.all(m.triangles[:,:,1]>=lower,axis=1));assert len(ids);q=Mesh(m.vertices,m.faces[ids],FRAME,m.attributes,dict(m.source,diagnostic_only=True,source_face_ids=ids.tolist(),no_source_change=True));path=OUT/('review-'+str(len(entries))+'.npz');q.save(path);e=dict(mesh=str(path),name=pid+(' [literal root faces]'if lower is not None else''),color=color);entries[pid]=e;return e
    reviewed=[]
    for pid in ['BP3D-height-replay:'+x for x in donors[:3]]:reviewed.append(reviewentry(pid,[.10,.72,.94,1]))
    for pid in ['BP3D-height-replay:'+x for x in donors[3:]]:reviewed.append(reviewentry(pid,[.10,.48,.70,1],1.700))
    for pid in native[:3]:reviewed.append(reviewentry(pid,[.95,.20,.13,1]))
    for pid in native[3:5]:reviewed.append(reviewentry(pid,[.70,.17,.10,1],1.700))
    for pid in native[9:11]:reviewed.append(reviewentry(pid,[.96,.87,.16,1]))
    spec=dict(frame=FRAME,meshes=reviewed,views=[dict(name='Actual unchanged ZA red vs existing-height-replayed BP3D cyan; inferior',direction=[0,-1,0],up=[0,0,1]),dict(name='Actual unchanged ZA red vs replayed BP3D cyan; anterior',direction=[0,0,1]),dict(name='Actual unchanged ZA red vs replayed BP3D cyan; oblique',direction=[.6,.5,1])],resolution=1200,wall_seconds=60)
    atomic_json(OUT/'comparison-review-spec.json',spec);queue=prepare(None,spec)
    after=json.loads(CAT.read_text());afterparts={p['id']:p for p in after['parts']};fixed_records_unchanged=all(identity(afterparts[pid])==sha for pid,sha in source_record_bindings.items());assert fixed_records_unchanged
    proof=dict(schema='anterior-source-interface-feasibility-stage1/v9',sourceCatalogueSha256=current_sha,currentCatalogueSha256AtFinish=digest_file(CAT),protectedCerebralRecordsLiteralDespiteIndependentInstall=fixed_records_unchanged,authorityTermsInventoried=len(inventory),sources=sources,interfaces=interfaces,actualInitializedDonorContextChecks=context_results,trueRootRecipientFitEstablished=False,pairedWallConstructionPermitted=False,sourceClosedCastsHaveNoEstablishedPairedRootAnnuli=True,initializedSourcesWereNotFittedOrRegistered=True,noGeometryOrLiveMetadataChanged=True,queue=queue,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',json.loads(json.dumps(proof,default=serial)))
    print(json.dumps(dict(output=str(OUT/'proof.json'),sourceCatalogueSha256=current_sha,authorityTerms=len(inventory),interfaces=[dict(side=r['side'],nativeIcaRimToInitializedBpAca_mm=1000*r['nativeIcaRootRim']['toReplayedBpAca']['distance_m'],nativeIcaRimToInitializedBpIca_mm=1000*r['nativeIcaRootRim']['toReplayedBpIca']['distance_m'],nativeAcaDistalRimToInitializedBpAca_mm=1000*r['nativeAcaDistalRim']['toReplayedBpAca']['distance_m'])for r in interfaces],contextCrossingCounts=[dict(donor=r['donor'],context=r['fixedContext'],pairs=r['actualBoundaryContactSummary']['crossing_face_pairs'])for r in context_results],reviewId=queue['request_id'],elapsedSeconds=proof['elapsedSeconds'])))

if __name__=='__main__':main()
