"""Source-only ankle landmarks, held face intervals and native review requests."""
from pathlib import Path
from collections import defaultdict,deque
import sys,json,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
PRIOR=BASE/'summary-v17.json'; OUT=BASE/'boundary-review-v18'; FRAME='model:local-ta98:assembled-metres'
REFERENCES=[
 dict(url='https://pubmed.ncbi.nlm.nih.gov/36691095/',doi='10.1186/s13037-023-00351-2',type='primary cadaveric study, 24 ankles',
      support='Distal GSV is anterior to medial malleolus; study uses actual prominent medial/lateral malleolar landmarks. Native lowest bone tip is only a coarse candidate, not the study landmark.'),
 dict(url='https://pubmed.ncbi.nlm.nih.gov/40176585/',doi='10.7417/CT.2025.5201',type='primary cadaveric study, 50 feet',
      support='Dorsal foot venous drainage continues into GSV/SSV at medial/lateral malleolar regions, with documented variants.'),
 dict(url='https://journals.viamedica.pl/folia_morphologica/article/download/FM.a2022.0053/66806',doi='10.5603/FM.a2022.0053',type='primary cadaveric case report',
      support='Published introduction explicitly describes lateral marginal vein into SSV and medial marginal vein into GSV. Search-index text verified; direct publisher PDF request timed out.'),
 dict(url='https://pmc.ncbi.nlm.nih.gov/articles/PMC4843236/',type='primary cadaveric case report',
      support='Normal SSV commences as lateral marginal continuation posterior to lateral malleolus; study reports a variant proximal course. Search-index text verified; PMC browser challenge prevented full page.')]

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def rings(m,faceids,archpoint):
    used=np.unique(m.faces[faceids]);v,inv=np.unique(m.vertices[used],axis=0,return_inverse=True)
    old=np.full(len(m.vertices),-1,int);old[used]=inv;f=old[m.faces[faceids]]
    loops=boundary_loops(f);assert len(loops)==2 and all(len(q)==12 for q in loops)
    first=min(loops,key=lambda q:np.linalg.norm(v[q].mean(0)-archpoint))
    graph=defaultdict(set)
    for tri in f:
        for a,b in zip(tri,np.roll(tri,-1)):graph[int(a)].add(int(b));graph[int(b)].add(int(a))
    d=np.full(len(v),-1,int);d[first]=0;queue=deque(map(int,first))
    while queue:
        x=queue.popleft()
        for y in graph[x]:
            if d[y]<0:d[y]=d[x]+1;queue.append(y)
    assert np.all(d>=0)and np.all(np.bincount(d)==12)
    return dict(centers=np.array([v[d==i].mean(0)for i in range(d.max()+1)]),faceRingMin=d[f].min(1),faceRingMax=d[f].max(1),
                sourceFaceIds=faceids,oldToPhysicalVertex=old,physicalRing=d,physicalVertices=v)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);prior=json.loads(PRIOR.read_text());sha=digest_file(CAT);assert sha==prior['sourceCatalogueSha256']
    a=json.loads(CAT.read_text());ps={p['id']:p for p in a['parts']};src={p['partId']:p for p in prior['sourceMeshes']};reader=AtlasReader();cache={};rows=[];reviews=[];snapshots=[]
    def get(pid):
        if pid not in cache:cache[pid]=reader.read_part(a,ps[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
        return cache[pid]
    def snap(m,pid,faceids,label,color,entries):
        faceids=np.asarray(faceids,int)
        if not len(faceids):return
        path=OUT/('literal-view-'+str(len(snapshots))+'.npz')
        view=Mesh(m.vertices,m.faces[faceids],FRAME,m.attributes,dict(m.source,sourceFaceIds=faceids.tolist(),sourcePartId=pid,diagnosticLiteralSubset=True,noAnatomicalBoundaryAccepted=True))
        view.save(path);entries.append(dict(mesh=str(path),name=label,color=color))
        snapshots.append(dict(path=str(path),partId=pid,sourceFaceIds=faceids,sourceGeometrySha256=m.geometry_sha,sourceAttributesSha256=m.attributes_sha,sourceRecordSha256=identity(ps[pid]),allPositionsIndicesNormalsUVsLiteral=True))
    for side in ['l','r']:
        entries=[];margins=[];landmarks=[];tips={}
        for kind,bid in [('medial','ZA:Tibia.'+side+':Bone'),('lateral','ZA:Fibula.'+side+':Bone')]:
            b=get(bid);ids=np.flatnonzero(b.vertices[:,1]==b.vertices[:,1].min());point=b.vertices[ids].mean(0);tips[kind]=point
            landmarks.append(dict(kind=kind,sourcePartId=bid,exactDistalMostSourceVertexIds=ids,point_m=point,
                uncertainty='Coarse anatomical tip candidate on distal native bone; prominent malleolar surface and exact named-vein transition need rendered review. No arbitrary nearest whole-bone centroid used.'))
        for term,pid,kind,color in [('A12.3.11.021','ZA:Great saphenous vein.'+side,'medial',[.05,.8,.95,1.]),('A12.3.11.020','ZA:Small saphenous vein.'+side,'lateral',[.9,.12,.63,1.])]:
            m=get(pid);sourceFaces=np.array(src[pid]['components'][0]['sourceFaceIds'],int);arch=next(q for q in prior['sourceMouthReplay']if q['sourcePartId']==pid)
            rr=rings(m,sourceFaces,np.array(arch['firstSourceRim']['center_m']));centers=rr['centers'];tip=tips[kind]
            # This is a coarse landmark-tied display proposal, not a final
            # anatomical assignment. Keep its ±5 mm band explicit for review.
            delta=centers[:,1]-tip[1];cross=np.flatnonzero(delta>=0)
            assert len(cross);boundary=int(cross[0]);band=np.flatnonzero(abs(delta)<=.005)
            below=sourceFaces[rr['faceRingMax']<=boundary];above=sourceFaces[rr['faceRingMin']>=boundary]
            mixed=sourceFaces[~((rr['faceRingMax']<=boundary)|(rr['faceRingMin']>=boundary))]
            assert not len(mixed) and len(below)+len(above)==len(sourceFaces)
            upperDisplay=above[m.triangles[above].mean(1)[:,1]<.22]
            snap(m,pid,below,side+' '+kind+' proposed foot interval (HELD); literal source faces',color,entries)
            snap(m,pid,upperDisplay,side+' retained '+ps[pid]['name']+' ankle course; not renamed',[.21,.31,.70,1.],entries)
            inherited=next((q for q in prior['additionalNativePackingChecks']if q['first']==pid),None)
            crossings=[]
            if inherited:
                pairs=np.asarray(inherited['sourceFacePairs'],int)
                for portion,ids in [('proposedFoot',below),('retainedOtherSource',np.setdiff1d(np.arange(len(m.faces)),below))]:
                    subset=pairs[np.isin(pairs[:,0],ids)];tri=m.triangles[np.unique(subset[:,0])]if len(subset)else np.empty((0,3,3))
                    crossings.append(dict(portion=portion,pairs=len(subset),sourceFacePairs=subset,
                        sourceTriangleExtent_m=[tri.min((0,1)),tri.max((0,1))]if len(tri)else None,
                        lowestCrossingAboveCoarseMalleolarTip_m=float(tri[:,:,1].min()-tip[1])if len(tri)else None))
            margins.append(dict(term=term,kind=kind,sourcePartId=pid,sourceComponent=0,sourceBankSha256=m.source['bank_sha256'],
                sourceGeometrySha256=m.geometry_sha,sourceAttributesSha256=m.attributes_sha,sourceCourseRingsFromArch=centers,
                coarseLandmarkPoint_m=tip,coarseBoundaryRing=boundary,boundaryRingCenter_m=centers[boundary],
                uncertaintyBandRingIds=band,uncertaintyBand_m=.005,proposedLiteralFootSourceFaceIds=below,
                retainedWholeSourceControl=True,retainedPrincipalSourceFaceIds=above,allOtherSourceComponentsRetained=True,
                noFaceDuplicatedOrRemoved=True,nonOwningReferenceOnly=True,geometryChanges=0,
                fullMouthReplay=arch,inheritedTibiaCrossingLocalization=crossings,
                anatomyStatus='HELD for native source render: malleolus-tip level is a coarse proposal, not an accepted histological or sharp anatomical boundary.'))
        for pid in ['ZA:Dorsal venous arch of foot.'+side,'ZA:Dorsal metatarsal veins.'+side,'ZA:Dorsal digital veins of foot.'+side]:
            m=get(pid);snap(m,pid,np.arange(len(m.faces)),side+' '+ps[pid]['name'],[.24,.48,.91,1.],entries)
        for pid in ['ZA:Tibia.'+side+':Bone','ZA:Fibula.'+side+':Bone','ZA:Talus.'+side+':Bone-1','ZA:Calcaneus.'+side+':Bone','ZA:First metatarsal bone.'+side+':Bone','ZA:Fifth metatarsal bone.'+side+':Bone']:
            if pid not in ps:continue
            b=get(pid);ids=np.flatnonzero(b.triangles.mean(1)[:,1]<.22)
            snap(b,pid,ids,side+' fixed native '+ps[pid]['name']+'; distal display crop only',[.72,.72,.68,1.],entries)
        spec=dict(frame=FRAME,meshes=entries,views=[dict(name=side+' ankle-foot SOURCE; dorsal',direction=[0,1,.4],up=[0,0,1]),
            dict(name=side+' ankle-foot SOURCE; medial',direction=[-1 if side=='l'else 1,.1,.1]),
            dict(name=side+' ankle-foot SOURCE; lateral',direction=[1 if side=='l'else-1,.1,.1]),
            dict(name=side+' ankle-foot SOURCE; coronal diagnostic section at medial arch/marginal mouths',direction=[0,0,1],section=dict(normal=[0,0,1],point=margins[0]['fullMouthReplay']['firstSourceRim']['center_m'],clear_outer=True))],resolution=1600,wall_seconds=60)
        atomic_json(OUT/(side+'-review-spec.json'),spec);request=prepare(None,spec);reviews.append(dict(side=side,review=request,spec=str(OUT/(side+'-review-spec.json'))))
        rows.append(dict(side=side,actualCoarseDistalBoneLandmarks=landmarks,marginProposals=margins))
    assert digest_file(CAT)==sha
    result=dict(schema='source-bound-marginal-vein-boundary-proposal/v18',sourceCatalogueSha256=sha,priorSummarySha256=digest_file(PRIOR),
        fullApplicableInventory=str(PRIOR),authorityParentCaveat='Both marginal-vein terms are listed beneath small saphenous vein in authoritative TA98. Preserve that terminology parent; physical medial continuation is GSV and is a separate anatomical relation.',
        primaryReferences=REFERENCES,rows=rows,diagnosticSnapshots=snapshots,reviews=reviews,
        status='Held for actual native render review. No physical or live catalogue changes. Reference-region acceptance is separate from flow-interface acceptance.',
        noNewAnatomicalGeometry=True,noCapOrPlug=True,noSourcePositionNormalUVChange=True,noInstall=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proposal.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(output=str(OUT/'proposal.json'),reviews=[dict(side=q['side'],id=q['review']['request_id'])for q in reviews],
        boundaries=[dict(side=q['side'],parts=[dict(term=p['term'],faces=len(p['proposedLiteralFootSourceFaceIds']),boundaryRing=p['coarseBoundaryRing'],crossings=[dict(portion=c['portion'],pairs=c['pairs'],lowestAboveTip_mm=1000*c['lowestCrossingAboveCoarseMalleolarTip_m']if c['lowestCrossingAboveCoarseMalleolarTip_m']is not None else None)for c in p['inheritedTibiaCrossingLocalization']])for p in q['marginProposals']])for q in rows],elapsedSeconds=result['elapsedSeconds'])))

if __name__=='__main__':main()
