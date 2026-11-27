"""Read-only current-frame branch retrieval; no carrier aliases or construction."""
from pathlib import Path
from collections import defaultdict,deque
import sys,time,json,sqlite3
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import patch_components,pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';CENSUS=ROOT/'reports/ta98-blood-circulation-inventory-installed-v12.json';AUTH=ROOT/'.local-models/ta98-authority.sqlite';OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/upper-limb-native-network-opportunities-v13';FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def ring_course(m,loops):
    if len(loops)!=2 or len(loops[0])!=len(loops[1]):return dict(regularUnbranchedCorrespondingRingCourse=False)
    g=defaultdict(set)
    for row in m.faces:
        for x,y in zip(row,np.roll(row,-1)):g[int(x)].add(int(y));g[int(y)].add(int(x))
    d=np.full(len(m.vertices),-1,int);d[loops[0]]=0;queue=deque(map(int,loops[0]))
    while queue:
        x=queue.popleft()
        for y in g[x]:
            if d[y]<0:d[y]=d[x]+1;queue.append(y)
    regular=bool(np.all(d>=0)and np.all(np.bincount(d)==len(loops[0])))
    if not regular:return dict(regularUnbranchedCorrespondingRingCourse=False)
    centers=np.array([m.vertices[d==i].mean(0)for i in range(d.max()+1)]);return dict(regularUnbranchedCorrespondingRingCourse=True,rings=len(centers),controlsPerRing=len(loops[0]),sourceCenterlineDerivedFromRingMeans_m=centers,approximateCourseLength_m=float(np.linalg.norm(np.diff(centers,axis=0),axis=1).sum()),noInternalSourceBranchJunction=True)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);source_sha=digest_file(CAT);a=json.loads(CAT.read_text());parts={p['id']:p for p in a['parts']};census=json.loads(CENSUS.read_text());assert census['source_catalogue_sha256']==source_sha;status={r['id']:r for r in census['terms']};db=sqlite3.connect('file:'+str(AUTH)+'?mode=ro',uri=True);db.row_factory=sqlite3.Row;terms={r['id']:dict(r)for r in db.execute('select * from ta98')};children=defaultdict(list)
    for r in terms.values():children[r['parent_id']].append(r['id'])
    roots=['A12.2.09.005','A12.2.09.013','A12.3.08.014'];seen=set();inventory=[]
    for root in roots:
        todo=[root]
        while todo:
            c=todo.pop()
            if c in seen:continue
            seen.add(c);r=status[c];inventory.append(dict(authority=terms[c],currentDirectPartIds=r['resolved_direct_part_ids'],currentDescendantOnlyPartIds=r['resolved_descendant_only_part_ids'],activeExactSourcePartIds=r['active_exact_source_part_ids'],suppressedExactSourcePartIds=r['suppressed_exact_source_part_ids'],currentMissingRequirements=r['missing_requirements'],representationStatus=r['representation_status']));todo.extend(children[c])
    atomic_json(OUT/'complete-subtree-inventory.json',dict(sourceCatalogueSha256=source_sha,censusSha256=digest_file(CENSUS),authoritySha256=digest_file(AUTH),roots=roots,terms=inventory))
    reader=AtlasReader();meshes={};trees={};bindings={};rows=[]
    def get(pid):
        if pid not in meshes:meshes[pid]=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=source_sha);bindings[pid]=identity(parts[pid])
        return meshes[pid]
    def tree(pid):
        if pid not in trees:trees[pid]=BVH(meshes[pid])
        return trees[pid]
    def closest(x,y):
        r=tree(x).closest_mesh(tree(y));r.update(firstPartId=x,secondPartId=y,firstTriangle_m=meshes[x].triangles[r['first_face_id']],secondTriangle_m=meshes[y].triangles[r['second_face_id']]);return r
    source_prefixes=['ZA:Thoraco-acromial artery.','ZA:Pectoral branches of thoraco-acromial artery.','O3M:Thoracoacromial artery Acromial br.','O3M:Thoracoacromial artery Deltoid br.','O3M:Thoracoacromial artery Pectoral br.','ZA:Subscapular artery.','ZA:Thoracodorsal artery.','ZA:Circumflex scapular artery.','ZA:Cephalic vein.','ZA:Basilic vein.','O3M:Superficial veins of upper limb.']
    for side in ['l','r']:
        landmark_ids=['ZA:Clavicle.'+side+':Bone','ZA:Scapula.'+side+':Bone','ZA:Clavicular part of deltoid muscle.'+side+':Abductor','ZA:Acromial part of deltoid muscle.'+side+':Abductor','ZA:Clavicular head of pectoralis major muscle.'+side+':Adductor','ZA:Pectoralis minor muscle.'+side]
        for pid in landmark_ids:get(pid)
        sourceids=[p+side for p in source_prefixes if p+side in parts]
        for pid in sourceids:
            m=get(pid);loops=boundary_loops(m.faces);components=patch_components(m,np.arange(len(m.faces)));comp=[]
            for i,ids in enumerate(components):
                ids=np.array(ids,int);q=Mesh(m.vertices,m.faces[ids],FRAME,m.attributes,m.source);key=pid+'#component-'+str(i);meshes[key]=q
                info=dict(component=i,sourceFaceIds=ids,triangles=len(ids),bounds_m=q.bounds,disjointLiteralFaceComponent=True)
                if pid.startswith('O3M:Superficial'):
                    trunk='ZA:Thoraco-acromial artery.'+side;get(trunk);info['toActualThoracoAcromialTrunk']=closest(key,trunk)
                comp.append(info)
            entry=dict(partId=pid,name=parts[pid]['name'],licence=parts[pid].get('license'),record=parts[pid],bankSha256=m.source['bank_sha256'],geometrySha256=m.geometry_sha,attributeSha256=m.attributes_sha,topology=topology(m),actualOpenBoundaryLoops=[dict(sourceVertexIds=q,points_m=m.vertices[q],centroid_m=m.vertices[q].mean(0))for q in loops],components=comp,sourceCourse=ring_course(m,loops))
            if ('Thoraco' in pid or 'Subscapular artery' in pid)and not pid.startswith('O3M:Superficial'):
                entry['actualLandmarkClosestTriangles']=[closest(pid,target)for target in landmark_ids]
            rows.append(entry)
        # The actual named parent/child surfaces are measured, never welded.
    junctions=[]
    for side in ['l','r']:
        for first,second in [('ZA:Thoraco-acromial artery.','ZA:Pectoral branches of thoraco-acromial artery.'),('ZA:Thoraco-acromial artery.','O3M:Thoracoacromial artery Acromial br.'),('ZA:Thoraco-acromial artery.','O3M:Thoracoacromial artery Deltoid br.'),('ZA:Subscapular artery.','ZA:Thoracodorsal artery.'),('ZA:Subscapular artery.','ZA:Circumflex scapular artery.')]:
            x,y=first+side,second+side
            if x not in meshes or y not in meshes:continue
            d,gate=pair_contacts(meshes[x],meshes[y],0.,first_bvh=tree(x),second_bvh=tree(y),max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget'];junctions.append(dict(firstPartId=x,secondPartId=y,closest=closest(x,y),actualContactGate=gate,actualZeroDistanceFacePairs=d['pairs'],anatomicalOrLumenContinuityEstablished=False))
    missing=[r['authority']['id']for r in inventory if not r['currentDirectPartIds']];after=json.loads(CAT.read_text());currentparts={p['id']:p for p in after['parts']};assert all(identity(currentparts[p])==h for p,h in bindings.items())
    result=dict(schema='upper-limb-native-source-network-opportunities/v13',sourceCatalogueSha256=source_sha,censusSha256=digest_file(CENSUS),authorityInventory=str(OUT/'complete-subtree-inventory.json'),inventoryTerms=len(inventory),missingDirectTermsInScope=missing,sourceMeshes=rows,actualNamedSourceJunctions=junctions,landmarkLimitation='Current catalogue has real clavicle/scapula/deltoid/pectoral meshes, but no independently identified current coracoid-process or acromion source-face region. Whole-scapula closest triangles are not silently called those landmarks.',candidateControlsPrepared=0,noBranchesInferredFromWholeCarrierNames=True,noPrivateDenseGeometryReadOrCopied=True,noGeometryOrLiveChanges=True,sourceRecordsStillLiteral=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(output=str(OUT/'proof.json'),inventoryTerms=len(inventory),missing=missing,sourceSummaries=[dict(id=r['partId'],triangles=r['topology']['triangles'],components=len(r['components']),openLoops=len(r['actualOpenBoundaryLoops']),regularUnbranched=r['sourceCourse'].get('regularUnbranchedCorrespondingRingCourse'))for r in rows],junctions=[dict(first=r['firstPartId'],second=r['secondPartId'],gap_mm=1000*r['closest']['distance_m'],crossings=r['actualContactGate']['crossing_face_pairs'])for r in junctions],elapsedSeconds=result['elapsedSeconds'])))
if __name__=='__main__':main()
