"""Read-only native ACA region controls; exact surface views, no new owners."""
from pathlib import Path
from collections import deque,defaultdict
import sys,time,json,copy
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,atomic_bytes,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

B=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1';OUT=B/'native-aca-bounded-references-v3';CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';EXPECTED='b8fdbe1c874df6ff62bf985d15708a4e4804eacf95f2475c0d3a6d07e82cec69';FRAME='model:local-ta98:assembled-metres'

def levels(mesh,loop):
    graph=defaultdict(set)
    for a,b,c in mesh.faces:
        for u,v in [(a,b),(b,c),(c,a)]:graph[int(u)].add(int(v));graph[int(v)].add(int(u))
    result=np.full(len(mesh.vertices),-1,int);result[loop]=0;q=deque(map(int,loop))
    while q:
        a=q.popleft()
        for b in graph[a]:
            if result[b]<0:result[b]=result[a]+1;q.append(b)
    counts=np.bincount(result);assert np.all(result>=0)and np.all(counts==len(loop)),'Native course is not an unbranched corresponding-ring tube'
    return result,np.array([mesh.vertices[result==j].mean(0)for j in range(result.max()+1)])

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;assert not OUT.exists();OUT.mkdir();a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};reader=AtlasReader();census=json.loads((ROOT/'reports/ta98-blood-circulation-inventory-installed-v11.json').read_text());terms={r['id']:r for r in census['terms']};children=defaultdict(list)
    for r in terms.values():children[r['parent_id']].append(r['id'])
    stack=['A12.2.07.022'];desc=[]
    while stack:
        c=stack.pop();desc.append(c);stack.extend(children[c])
    inventory=dict(schema='native-aca-authority-and-source-inventory/v1',root='A12.2.07.022',sourceCatalogueSha256=EXPECTED,censusSha256=digest_file(ROOT/'reports/ta98-blood-circulation-inventory-installed-v11.json'),authoritySha256=census['authority_cache_sha256'],terms=[terms[c]for c in sorted(desc)],wholeNativeAssembly=['ZA:Anterior cerebral artery.l','ZA:Anterior cerebral artery.r'],plannedMissingTermControls=['A12.2.07.023','A12.2.07.033'],noAdditionalUnseenBranchesClaimed=True,noGeometryReconstruction=True)
    atomic_json(OUT/'authority-inventory.json',inventory)
    atomic_json(ROOT/'reports/ta98-anterior-cerebral-native-region-inventory-v3.json',inventory)
    comm=reader.read_part(a,parts['ZA:Anterior communicating artery'],frame=FRAME,catalogue_path=CAT);results=[];records={};data=bytearray();ancestry={};termgroups=defaultdict(list)
    for side in ['l','r']:
        pid='ZA:Anterior cerebral artery.'+side;p=parts[pid];m=reader.read_part(a,p,frame=FRAME,catalogue_path=CAT);ica=reader.read_part(a,parts['ZA:Internal carotid artery.'+side],frame=FRAME,catalogue_path=CAT);loops=boundary_loops(m.faces);assert len(loops)==2 and all(len(q)==12 for q in loops)
        ends=np.array([m.vertices[q].mean(0)for q in loops]);hit,_=BVH(ica).closest_points(ends);proximal=int(hit['distance_m'].argmin());ring,centers=levels(m,loops[proximal]);assert hit['distance_m'][1-proximal]>hit['distance_m'][proximal]+.002,'Proximal ICA end is ambiguous'
        contacts,check=pair_contacts(m,comm,0.,max_pairs=250000,max_tests=2000000);assert check['complete_within_pair_budget'];atomic_npz(OUT/('native-ACo-contact-'+side+'.npz'),**contacts)
        actual=contacts['pairs'];assert len(actual)>0,'No actual native ACom crossing/contact footprint; cannot assign A1/A2 by a coordinate guess'
        cf=np.unique(actual[:,0]);low=int(ring[m.faces[cf]].min());high=int(ring[m.faces[cf]].max());split=(low+high)//2;assert 2<=split<=ring.max()-2,'Actual communicating footprint is at an endpoint, not a usable A1/A2 boundary'
        selected=[np.flatnonzero(ring[m.faces].max(1)<=split),np.flatnonzero(ring[m.faces].min(1)>=split)];assert len(np.intersect1d(*selected))==0 and np.array_equal(np.sort(np.r_[*selected]),np.arange(len(m.faces)))
        raw=reader.asset(a['chunks'][p['chunk']]['url']).read_bytes();n=p['vertexCount'];rawfaces=reader.view(np.frombuffer(raw,dtype='u1'),p['indices'],p['indexCount'],'<u4').reshape(-1,3);assert len(rawfaces)==len(m.faces)
        buffers={}
        for key,width in [('positions',12),('normals',6),('normalSigns',1),('uvs',8),('colors',3)]:
            if key not in p:continue
            data.extend(b'\0'*((-len(data))%4));buffers[key]=len(data);data.extend(raw[p[key]:p[key]+n*width])
        for term,ids in zip(['A12.2.07.023','A12.2.07.033'],selected):
            rid='TA98ARTREF:male:'+term+':'+side;data.extend(b'\0'*((-len(data))%4));offset=len(data);data.extend(rawfaces[ids].astype('<u4').tobytes());record=copy.deepcopy(p);record.update(buffers);record.update(id=rid,conceptId=rid,name=terms[term]['name_en'].capitalize()+(' (left)'if side=='l'else' (right)'),indices=offset,indexCount=int(len(ids)*3),triangles=len(ids),chunk=0,ta98Term=term,ta98Kind='exact',nativeVascularRegionOf=pid,nonOwningSurfaceReference=True,pairedTissue=False,sectionSolid=False,defaultHidden=True,referenceRole='bounded-native-arterial-course-surface',sourceId=pid,sourceFaceIds=ids.tolist(),bounds=[m.triangles[ids].min((0,1)).tolist(),m.triangles[ids].max((0,1)).tolist()],coverageLimitation='Non-owning bounded surface reference on the unchanged native ACA course. A1/A2 boundary uses its actual native ACom contact footprint and one existing ring; finite wall layers, an unobstructed communicating lumen, and junction perfusion have not been established.')
            records[rid]=record;ancestry[rid]=dict(sourcePartId=pid,sourceBankSha256=m.source['bank_sha256'],sourceGeometrySha256=m.geometry_sha,sourceAttributeSha256=m.attributes_sha,sourceFaceIds=ids.tolist(),sourceContactFootprintRingRange=[low,high],boundarySourceRing=split,boundaryPoints=m.vertices[ring==split].tolist(),representation='non-owning-bounded-tissue-surface-reference',actualAuthorityEntityType=terms[term]['authority_type'],notAnOpeningOrNewTubeSolid=True);termgroups[term].append(rid)
        results.append(dict(sourcePartId=pid,sourceGeometrySha256=m.geometry_sha,sourceBankSha256=m.source['bank_sha256'],sourceLicense=p['license'],sourceTriangles=len(m.faces),sourceRings=len(centers),proximalIcaSurfaceDistances_m=hit['distance_m'].tolist(),actualAComContact=check,communicatingFootprintRingRange=[low,high],boundarySourceRing=split,boundaryWithinActualContactFootprint=True,boundaryCoordinate_m=centers[split].tolist(),boundedRegionFaceCounts=[len(q)for q in selected],fullSourceFaceCoverNoOverlap=True,sourceCourseAndCalibreUnchanged=True,noNewPhysicalOwners=True))
    bank=OUT/'bounded-reference-views.bin';atomic_bytes(bank,bytes(data));cat=copy.deepcopy(a);cat['chunks']=[dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=len(data))];cat['parts']=list(records.values())
    for rid,p in records.items():
        new=reader.read_part(cat,p,frame=FRAME);sourceid=ancestry[rid]['sourcePartId'];old=reader.read_part(a,parts[sourceid],frame=FRAME,catalogue_path=CAT);ids=np.array(ancestry[rid]['sourceFaceIds']);assert np.array_equal(new.triangles,old.triangles[ids]);ow=old.attributes['normals'][old.faces[ids]];nw=new.attributes['normals'][new.faces];assert np.array_equal(ow,nw)
    atomic_json(OUT/'records.json',records);atomic_json(OUT/'source-face-ancestry.json',ancestry);atomic_json(OUT/'term-controls.json',dict(schema='native-arterial-bounded-reference-controls/v1',sourceCatalogueSha256=EXPECTED,bankSha256=digest_file(bank),physicalOwnersAdded=0,terms=[dict(term=c,parent=terms[c]['parent_id'],representation='tissue-region-surface-reference',partIds=ids,nonOwning=True)for c,ids in sorted(termgroups.items())]))
    proof=dict(schema='native-aca-bounded-reference-proof/v3',sourceCatalogueSha256=EXPECTED,bank=str(bank),bankSha256=digest_file(bank),bankBytes=len(data),terms=2,referenceControls=4,physicalOwnersAdded=0,allCopiedPositionNormalBuffersLiteral=True,allReferenceTrianglesAndCornerWordsLiteralSourceSubsets=True,actualLandmarkSource='ZA:Anterior communicating artery with actual ICA proximal-end witnesses',wholeNativeArterySurfacesRemainUnchanged=True,notFiniteWallOrLumenAcceptance=True,sourceInterfacesStillNeedPatentBloodSpaceReconstruction=True,rows=results,installed=False,elapsedSeconds=time.monotonic()-began);atomic_json(OUT/'proof.json',proof);assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(terms=2,controls=4,physicalOwnersAdded=0,bankBytes=len(data),bankSha256=proof['bankSha256'],sourceBoundaryRings=[dict(id=r['sourcePartId'],ring=r['boundarySourceRing'],contact=r['communicatingFootprintRingRange'])for r in results],elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
