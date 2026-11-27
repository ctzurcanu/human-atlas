"""Read-only coherent closure contract and two bounded native review queues."""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
OUT=BASE/'whole-foot-moving-closure-contract-v30'
FRAME='model:local-ta98:assembled-metres'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    a=json.loads(CAT.read_bytes());ps={p['id']:p for p in a['parts']};reader=AtlasReader();bindings={};stages=[];samples=[];nervecontacts=[]
    interfaces=json.loads((BASE/'whole-foot-native-interfaces-v26/evidence.json').read_bytes());runtime=json.loads((BASE/'whole-foot-runtime-qualification-v28.json').read_bytes());retained={r['id']:r['retained']for r in runtime['rows']}
    for side in ['l','r']:
        for group in ['marginal','dorsal','deep']:
            for mode in ['rays','contacts']:
                p=BASE/f'whole-foot-compartment-v29-{side}-{group}-{mode}/evidence.json';d=json.loads(p.read_bytes());assert d['sourceCatalogueSha256']==EXPECTED and d['stage']==dict(side=side,group=group,mode=mode)
                stages.append(dict(path=str(p),sha256=digest_file(p),stage=d['stage'],elapsedSeconds=d['elapsedSeconds'],complete=True));samples.extend(d['samples']);nervecontacts.extend(d['actualNeighborNerveContacts'])
    # No transform is solved or selected. These samples expose coarse
    # corridor variation and open reference uncertainties for the next method.
    corridors=[]
    for pid in sorted(set(q['sourcePartId']for q in samples)):
        q=[x for x in samples if x['sourcePartId']==pid];iv=[x for x in q if x['coarseCentreTranslationIntervalAlongRay_m']is not None]
        low=max(iv,key=lambda x:x['coarseCentreTranslationIntervalAlongRay_m'][0])if iv else None
        high=min(iv,key=lambda x:x['coarseCentreTranslationIntervalAlongRay_m'][1])if iv else None
        corridors.append(dict(sourcePartId=pid,kind=q[0]['kind'],samples=len(q),samplesWithBothRayBoundaries=len(iv),samplesWithPositiveRayInterval=sum(x['positiveWidthIfBothHits']for x in q),
            conservativeUniformSameAxisInterval_m=[low['coarseCentreTranslationIntervalAlongRay_m'][0],high['coarseCentreTranslationIntervalAlongRay_m'][1]]if iv else None,
            lowWitness=low,highWitness=high,feasible3DRigidTransformEstablished=False,globalInfeasibilityNotClaimed=True,
            rayOnlyLimit='Directional source-ring circumsphere support, current one-sided fascia and open skin; not an exact perpendicular tissue depth or occupied 3D corridor.'))
    closure=[]
    for side in ['l','r']:
        superficial=['ZA:'+p+side for p in ['Dorsal venous arch of foot.','Dorsal metatarsal veins.','Dorsal digital veins of foot.']]
        transition=['ZA:'+p+side for p in ['Great saphenous vein.','Small saphenous vein.']]
        perforators=['ZA:Intercapitular veins of foot.'+side]
        deep=['ZA:'+p+side for p in ['Plantar venous arch.','Plantar metatarsal veins.','Plantar digital veins.','Medial plantar veins.','Lateral plantar veins.']]
        sourcegraphs=[q for q in interfaces['actualSourceSurfaceContactPairs']if q['first'].endswith('.'+side)and q['gate']['returned_triangle_pairs']]
        closure.append(dict(side=side,wholeAssembly='Superficial dorsal and marginal foot venous network with explicitly tethered perforators and proximal saphenous transitions.',
            proposedCoherentMovingSurfaceOwners=superficial,proximalTransitionOwners=transition,
            perforatorBridgeOwners=perforators,protectedDeepPlantarRecipientOwners=deep,
            fixedProximalRecipientOwners=['ZA:Femoral vein.'+side,'ZA:Popliteal vein.'+side,'ZA:Posterior tibial veins.'+side,'ZA:Anterior tibial veins.'+side,'ZA:Fibular veins.'+side],
            actualSourceSurfaceContactGraph=sourcegraphs,
            sourceBoundaryMouthReferences=[q for q in interfaces['nearbyCompleteSourceMouthPairs']if any(q['first'].startswith(x+'#')or q['second'].startswith(x+'#')for x in superficial+transition+perforators)],
            actualHighestProximalSourceEnds=[q for q in interfaces['actualProximalTissueEnds']if q['trunkPartId'].endswith('.'+side)],
            exactFlowOrFiniteWallInterfaceEstablished=False,wholeSourceFamilyIsNotASingleClosedMesh=True,
            motionNotAuthorizedForDeepStructuresByMereContact=True,
            couplingRules=['All three superficial source groups and both marginal continuations must share one coordinated anatomical target corridor; the ankle display boundary does not fix a bad current vessel station.',
                'Proper station transport may bridge from the repositioned superficial family into literal retained proximal GSV/SSV source controls; choose fixed transition support from measured actual calf/thigh context, not an arbitrary foot plane.',
                'Every original metatarsal/digital/arch junction must be explicitly reconstructed with common outer-interface coordinates and a shared normal chart; current overlapping open casts are not annuli.',
                'Intercapitular source cylinders span superficial and deep systems. Their superficial endpoints may follow the dorsal family only through an explicit tethered course; keep deep plantar recipient coordinates fixed until the corresponding actual junction is defined.',
                'Do not move fascia, bones, nerves, arterial courses, neighboring muscle/tendon anatomy, skin or heart to create room. Deep plantar veins remain in their anatomical intermuscular compartment.',
                'Do not infer closed loop containers from source names or a hierarchy parent; every connected source component and every true boundary loop retains its own identity.']))
    # Original-only rendering. Crop triangles for a local camera, without
    # inventing caps or moving any source point, normal, UV or bone.
    entries=[];snapshots=[]
    ids=[r['id']for r in runtime['rows']if r['retained']and r['id']in ps and ps[r['id']].get('bounds')and ps[r['id']]['bounds'][1][0]>.05 and ps[r['id']]['bounds'][0][0]<.18]
    ids=[pid for pid in ids if not pid.endswith('.r')and not any(z in pid for z in ['.r:',':right'])]
    ids+=['ZA:'+p+'.l'+suffix for p,suffix in [('Calcaneus',':Bone-5'),('Talus',':Bone-1'),('First metatarsal bone',':Bone'),('Fifth metatarsal bone',':Bone')]]
    ids=list(dict.fromkeys(ids));assert len(ids)<=100
    low=np.array([.045,.0,-.12]);high=np.array([.18,.145,.16])
    for pid in ids:
        p=ps[pid];m=reader.read_part(a,p,frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED);t=m.triangles;keep=np.flatnonzero(np.all(t.max(1)>=low,axis=1)&np.all(t.min(1)<=high,axis=1))
        if not len(keep):continue
        path=OUT/(identity(pid)[:20]+'.npz');Mesh(m.vertices,m.faces[keep],FRAME,m.attributes,dict(m.source,originalSourceFaceIds=keep.tolist(),literalSourceDiagnosticCrop=True,noCapsAdded=True)).save(path)
        bindings[pid]=dict(recordSha256=identity(p),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'])
        color=[.08,.30,.85,1.]if p['system']=='venous'else([1.,.8,.04,1.]if p['system']=='nervous'else([.8,.62,.5,.18]if p['system']=='integumentary'else([.75,.47,.43,.24]if p['system']=='fascia'else[.65,.65,.60,.34])))
        entries.append(dict(mesh=str(path),name='ORIGINAL '+pid+'; literal uncapped diagnostic crop',color=color))
        snapshots.append(dict(partId=pid,path=str(path),sha256=digest_file(path),sourceFaceIds=keep,sourceGeometrySha256=m.geometry_sha,sourceAttributesSha256=m.attributes_sha))
    reviews=[]
    for name,point,direction in [('whole-superficial-deep', [.107,.043,.015],[1.,.8,1.]),('lateral-ankle-compartment',[.110,.045,-.044],[1.,.25,-.2])]:
        spec=dict(frame=FRAME,meshes=entries,views=[dict(name='ORIGINAL whole left foot venous family, fascia, skin and nerves; '+name,direction=direction),
            dict(name='ORIGINAL true uncapped axial compartment section; '+name,direction=[0,1,0],section=dict(normal=[0,1,0],point=point,clear_outer=True))],resolution=1400,wall_seconds=60)
        path=OUT/(name+'-review-spec.json');atomic_json(path,spec);request=prepare(None,spec);reviews.append(dict(spec=str(path),request=request))
    result=dict(schema='whole-foot-venous-moving-closure-and-compartment-contract/v30',sourceCatalogueSha256=EXPECTED,
        sourceInterfaceInventorySha256=digest_file(BASE/'whole-foot-native-interfaces-v26/evidence.json'),runtimeQualificationSha256=digest_file(BASE/'whole-foot-runtime-qualification-v28.json'),
        completedReadOnlyStages=stages,sourceBindings=bindings,coherentClosures=closure,coarseAvailableCorridors=corridors,
        actualNativeNerveContacts=nervecontacts,reviewSnapshots=snapshots,reviews=reviews,
        primaryReferences=[dict(url='https://pubmed.ncbi.nlm.nih.gov/40176585/',type='50-foot cadaveric study',support='Superficial dorsal drainage into GSV/SSV at malleoli; documented perforating and retinacular variants do not establish arbitrary native intersections as normal.'),
            dict(url='https://pubmed.ncbi.nlm.nih.gov/22847928/',type='400-foot injected cadaveric dissection',support='Deep plantar reservoir between plantar muscles; medial perforators connect deep plantar and superficial medial marginal systems.'),
            dict(url='https://anatomypubs.onlinelibrary.wiley.com/doi/10.1002/ar.21085',type='10-foot cadaveric dissection',support='Separate deep plantar venous arch/interconnections and dorsal connections; anatomical source does not support moving all foot vessels as a superficial group.')],
        singleWholeFootRigidCorrectionEstablished=False,rigidOrProperFieldNotYetSelected=True,
        nextAction='Review actual whole compartment. Then define a source-bound superficial moving family plus explicit tethered proximal/perforator transitions, and a constrained section-corridor graph. Solve one coupled field with all fixed current obstacles, source calibre and common junction/normal constraints. Do not repeat SSV-only lateral offsets or install a medial-only shortcut.',
        excludedClaims=['No accepted anatomical repair','No constructed new geometry','No no-intersection or occupied-volume proof from ray samples','No completed finite wall or blood lumen','No source visibility or atlas change'],
        anatomyAcceptance=False,installed=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'contract.json',json.loads(json.dumps(result,default=clean)));assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(contract=str(OUT/'contract.json'),readOnlyStages=len(stages),samples=len(samples),nervePairs=len(nervecontacts),reviewObjects=len(entries),elapsedSeconds=result['elapsedSeconds'],
        reviews=[q['request']['request_id']for q in reviews],actualNerveCrossingPairs=sum(q['gate']['crossing_face_pairs']for q in nervecontacts))))

if __name__=='__main__':main()
