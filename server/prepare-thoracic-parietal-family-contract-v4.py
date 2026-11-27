"""Executable thoracic family obligations and original-only source review.

The existing parent stem has skeletal crossings. No new branch geometry is
constructed before the parent/interface correction passes the same context.
"""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
OUT=BASE/'whole-family-contract-v4'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    inventory=json.loads((BASE/'source-inventory.json').read_bytes());measure=json.loads((BASE/'ita-actual-interface-landmarks-v3/evidence.json').read_bytes());runtime=json.loads((BASE/'actual-runtime-context-v2.json').read_bytes())
    a=json.loads(CAT.read_bytes());ps={p['id']:p for p in a['parts']};reader=AtlasReader();source={q['id']:q for q in inventory['sources']};bindings={};localization=[];reviews=[]
    for row in measure['rows']:
        side=row['side'];pid=row['parentPartId'];d=np.load(BASE/'ita-actual-interface-landmarks-v3'/(side+'-ita-source-stations.npz'));station=np.full(len(d['physical_vertices']),-1,int)
        for j,rr in enumerate(d['rings']):station[rr]=j
        fs=station[d['physical_faces']]
        for contact in row['wholeStemContextContacts']:
            if not contact['gate']['crossing_face_pairs']:continue
            pairs=np.array(contact['sourceFacePairs'],int);originalSourceFaces=np.unique(pairs[:,0]);interval=np.unique(np.c_[fs[originalSourceFaces].min(1),fs[originalSourceFaces].max(1)],axis=0)
            localization.append(dict(side=side,parentPartId=pid,fixedNeighborPartId=contact['partId'],actualCrossingFacePairs=contact['gate']['crossing_face_pairs'],parentSourceFaces=originalSourceFaces,
                fixedNeighborSourceFaces=np.unique(pairs[:,1]),sourceStationIntervals=interval,physicalDefectNotMerelyMissingTA98Label=True))
    obligations=[]
    for q in inventory['terms']:
        term=q['id'];name=q['name_en'];entry=dict(term=term,parent= q['parent_id'],name=name,currentSources=q['sourcePartIds'],status=q['representationStatus'],classification=q['classification'],
            explicitTakeoffEstablished=False,recipientEntryPatchEstablished=False,newMeshConstructed=False)
        if term.startswith('A12.2.08.'):
            entry.update(wholeParentAssemblies=['ZA:Internal thoracic artery.l','ZA:Internal thoracic artery.r'],originEvidence='Actual native unbranched ITA profiles plus a primary-reference-supported intercostal or terminal level; a new branch takeoff requires an explicitly selected anatomical field and source socket.',
                proposedRecipientCategory=('actual sternal periosteal surface'if term.endswith('.035')else'actual intercostal muscle tissue'if term.endswith('.039')else'current specific recipient tissue, unresolved until named interface is measured'),
                sourcePrerequisite='Repair the existing ITA/sternal and fourth-costal-cartilage crossings, with unchanged parent profile calibre and fixed native terminal/subclavian anatomy. Parent roots/terminals are overlapping casts, not matching annuli.')
        else:
            entry.update(wholeParentAssemblies=['ZA:Thoracic aorta','ZA:Posterior intercostal arteries.l','ZA:Posterior intercostal arteries.r','ZA:Subcostal artery.l','ZA:Subcostal artery.r'],
                originEvidence='Actual source aortic/intercostal component and the referenced vertebral/intercostal level. Primary cadaver evidence varies in branch origins; no equal-spacing or nearest-container heuristic.',
                proposedRecipientCategory='Named thoracic, parietal, cutaneous or spinal recipient; explicit surface/canal/root interface required per branch.',
                sourcePrerequisite='Resolve physical source roots and complete fixed context for the particular intercostal component before adding dorsal, spinal, radicular or collateral branches.')
        if '.11.01' in term and term in ['A12.2.11.011','A12.2.11.012','A12.2.11.013','A12.2.11.014','A12.2.11.015','A12.2.11.016']:
            entry['additionalConstraint']='Spinal canal, cord and root geometry stay fixed. A source nerve-root name does not establish a root vessel mouth, spinal arterial topology or a safe foraminal corridor.'
        obligations.append(entry)
    # One original-source diagnostic queue: side stem and its actual bone,
    # cartilage and recipient muscle context, never a replacement UI scene.
    entries=[];snapshots=[]
    ids=['ZA:Internal thoracic artery.l','ZA:Left subclavian artery','ZA:Musculophrenic artery.l','ZA:Superior epigastric artery.l',
        'ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6','ZA:Transversus thoracis muscle.l']
    ids+=['ZA:Costal cartilage of '+n+' rib.l'for n in ['second','third','fourth','fifth','sixth','seventh']]
    for pid in ids:
        p=ps[pid];m=reader.read_part(a,p,frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED);path=OUT/('source-'+str(len(entries))+'.npz');m.save(path)
        entries.append(dict(mesh=str(path),name='ORIGINAL '+pid,color=[.9,.12,.17,1.]if p['system']=='arterial'else([.78,.60,.43,.65]if p['system']=='muscular'else[.73,.73,.68,1.])))
        binding=dict(partId=pid,path=str(path),sha256=digest_file(path),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256']);snapshots.append(binding);bindings[pid]=binding
    point=np.load(BASE/'ita-actual-interface-landmarks-v3/l-ita-source-stations.npz')['centers'][70].tolist()
    spec=dict(frame=FRAME,meshes=entries,views=[dict(name='ORIGINAL left ITA with actual fixed sternum/costal/transversus context',direction=[.7,.2,-1.]),
        dict(name='ORIGINAL true uncapped axial section through ITA-sternal source crossings',direction=[0,1,0],section=dict(normal=[0,1,0],point=point,clear_outer=True))],resolution=1400,wall_seconds=60)
    atomic_json(OUT/'source-review-spec.json',spec);request=prepare(None,spec);reviews.append(dict(spec=str(OUT/'source-review-spec.json'),request=request))
    result=dict(schema='thoracic-parietal-arterial-whole-family-interface-contract/v4',sourceCatalogueSha256=EXPECTED,
        inventorySha256=digest_file(BASE/'authority-inventory.json'),actualSourceInventorySha256=digest_file(BASE/'source-inventory.json'),
        actualRuntimeContextSha256=digest_file(BASE/'actual-runtime-context-v2.json'),actualInterfaceLandmarksSha256=digest_file(BASE/'ita-actual-interface-landmarks-v3/evidence.json'),
        all35TermObligations=obligations,knownOriginalParentCrossings=localization,sourceBindings=bindings,originalSourceReviewSnapshots=snapshots,reviews=reviews,
        fixedOwners=['All current skeleton, sternum and costal cartilages','All 41 installed heart assembly owners and original coronaries','Current muscles, diaphragm, fascia, nerves, skin and existing vascular recipient roots'],
        bestGroundedFirstFamily=dict(terms=['A12.2.08.035','A12.2.08.039','A12.2.08.036','A12.2.08.037'],
            names=['sternal branches','anterior intercostal branches','perforating branches','medial mammary branches'],
            reason='Four missing macroscopic arterial groups share the current native ITA and physically modeled sternum/intercostal/skin recipient neighborhood; a coordinated short anterior-thoracic parent/branch family avoids inventing a long pericardiophrenic route through an unmodeled recipient.',
            readyForBranchConstruction=False,seriousObstacle='Original bilateral ITA has 89 sternum and 30 fourth-costal-cartilage crossing face pairs per side; new takeoffs cannot be accepted on that parent. Its current terminal and root casts also lack matching annuli.'),
        pericardiacophrenicHold=dict(term='A12.2.08.034',reason='Only explicit LOCAL phrenic-nerve source records are suppressed by the actual normal loader; no explicit pericardium source mesh exists. Exact actual nerve/pericardial recipient frame and diaphragmatic entry remain unresolved.'),
        constructionPrerequisites=['Review actual existing source stem/costal cuts. Identify one native anatomical ITA corridor, using fixed posterior sternal/cartilage surfaces and current intercostal/transversus context.',
            'Use one whole-family proper section transport preserving all source profile areas/diameters and literal roots outside declared support; a positive fit cannot move the skeleton, heart or original coronaries.',
            'Create actual reciprocal common parent/child outer interfaces and shared normals. Source outer casts alone do not establish a finite vascular wall or fluid lumen.',
            'Specify each authored branch takeoff by measured source station/socket and each endpoint by actual recipient tissue patch and depth; no whole-organ centroid, floating line or generic container counts as supply.',
            'Check all 656 current available broad-context records through an appropriate actual candidate AABB, not solely the selected bones; then actual emitted surface/section review and guarded local install.'],
        noBranchCourseOrGeometryConstructed=True,noAtlasChanges=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'contract.json',json.loads(json.dumps(result,default=clean)));assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(contract=str(OUT/'contract.json'),terms=len(obligations),elapsedSeconds=result['elapsedSeconds'],review=request['request_id'],
        actualCrossingIntervals=[dict(side=r['side'],neighbor=r['fixedNeighborPartId'],pairs=r['actualCrossingFacePairs'],intervals=r['sourceStationIntervals'].tolist())for r in localization])))

if __name__=='__main__':main()
