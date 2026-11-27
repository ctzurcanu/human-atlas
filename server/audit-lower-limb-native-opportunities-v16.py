"""Bounded source-only lower-limb branch retrieval; no geometry writes or aliases."""
from pathlib import Path
import sys, time, json, sqlite3
from collections import defaultdict, deque
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file, identity
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import patch_components, pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
CENSUS=ROOT/'reports/ta98-blood-circulation-inventory-installed-v12.json'
AUTH=ROOT/'.local-models/ta98-authority.sqlite'
OUT=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16/summary.json'
FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray): return x.tolist()
    if isinstance(x,np.generic): return x.item()
    raise TypeError(type(x).__name__)

def physical(m):
    v, inv=np.unique(m.vertices,axis=0,return_inverse=True)
    return Mesh(v,inv[m.faces],m.frame,{},dict(diagnosticExactWeldingOnly=True))

def course(m,loops):
    if len(loops)!=2 or len(loops[0])!=len(loops[1]): return None
    graph=defaultdict(set)
    for row in m.faces:
        for x,y in zip(row,np.roll(row,-1)): graph[int(x)].add(int(y)); graph[int(y)].add(int(x))
    d=np.full(len(m.vertices),-1,int); d[loops[0]]=0; queue=deque(map(int,loops[0]))
    while queue:
        x=queue.popleft()
        for y in graph[x]:
            if d[y]<0: d[y]=d[x]+1; queue.append(y)
    if np.any(d<0) or not np.all(np.bincount(d)==len(loops[0])): return None
    centers=np.array([m.vertices[d==i].mean(0) for i in range(d.max()+1)])
    return dict(rings=len(centers),verticesPerRing=len(loops[0]),literalRingCenters_m=centers,
                approximateCourseLength_m=float(np.linalg.norm(np.diff(centers,axis=0),axis=1).sum()),
                unbranchedRegularSourceCourse=True)

def main():
    began=time.monotonic(); sha=digest_file(CAT); a=json.loads(CAT.read_text()); ps={p['id']:p for p in a['parts']}
    census=json.loads(CENSUS.read_text()); assert census['source_catalogue_sha256']==sha
    statuses={t['id']:t for t in census['terms']}
    db=sqlite3.connect('file:'+str(AUTH)+'?mode=ro',uri=True); db.row_factory=sqlite3.Row
    authority=[dict(r) for r in db.execute("select * from ta98 where id like 'A12.2.16.%' or id like 'A12.3.11.%' order by id")]
    inventory=[]
    for t in authority:
        s=statuses[t['id']]
        inventory.append(dict(authority=t,directPartIds=s['resolved_direct_part_ids'],
            activeExactPartIds=s['active_exact_source_part_ids'],suppressedExactPartIds=s['suppressed_exact_source_part_ids'],
            representationStatus=s['representation_status'],missingRequirements=s['missing_requirements']))
    reader=AtlasReader(); meshes={}; trees={}; bindings={}; sources=[]; tests=[]
    def get(pid):
        if pid not in meshes:
            meshes[pid]=reader.read_part(a,ps[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
            bindings[pid]=identity(ps[pid])
        return meshes[pid]
    def tree(pid):
        if pid not in trees: trees[pid]=BVH(meshes[pid])
        return trees[pid]
    def test(x,y,contact=True):
        get(x); get(y); c=tree(x).closest_mesh(tree(y))
        c.update(firstTriangle_m=meshes[x].triangles[c['first_face_id']],secondTriangle_m=meshes[y].triangles[c['second_face_id']])
        r=dict(first=x,second=y,closest=c)
        if contact:
            d,g=pair_contacts(meshes[x],meshes[y],0.,first_bvh=tree(x),second_bvh=tree(y),max_pairs=250000,max_tests=2000000)
            assert g['complete_within_pair_budget']; r.update(actualContactGate=g,actualContactSourceFacePairs=d['pairs'])
        tests.append(r); return r
    for side in ['l','r']:
        prefixes=['ZA:Great saphenous vein.','ZA:Small saphenous vein.','ZA:Dorsal venous arch of foot.',
                  'ZA:Dorsal metatarsal veins.','ZA:Dorsal digital veins of foot.','O3M:Dorsal venous netwo.']
        for prefix in prefixes:
            pid=prefix+side; m=get(pid); welded=physical(m); loops=boundary_loops(welded.faces)
            components=patch_components(m,np.arange(len(m.faces)))
            sources.append(dict(partId=pid,name=ps[pid]['name'],sourceRecord=ps[pid],recordSha256=bindings[pid],
                bankSha256=m.source['bank_sha256'],geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,
                topology=topology(m),bounds_m=m.bounds,exactCoordinateBoundaryLoops=[dict(vertexIds=q,points_m=welded.vertices[q],center_m=welded.vertices[q].mean(0))for q in loops],
                components=[dict(component=i,sourceFaceIds=ids,triangles=len(ids),bounds_m=Mesh(m.vertices,m.faces[ids],FRAME).bounds)for i,ids in enumerate(components)],
                sourceCourse=course(welded,loops),notPairedWallGeometry=True))
        for x,y in [('ZA:Great saphenous vein.','ZA:Dorsal venous arch of foot.'),
                    ('ZA:Small saphenous vein.','ZA:Dorsal venous arch of foot.'),
                    ('ZA:Dorsal metatarsal veins.','ZA:Dorsal venous arch of foot.'),
                    ('ZA:Dorsal digital veins of foot.','ZA:Dorsal metatarsal veins.'),
                    ('O3M:Dorsal venous netwo.','ZA:Great saphenous vein.'),
                    ('O3M:Dorsal venous netwo.','ZA:Small saphenous vein.'),
                    ('O3M:Dorsal venous netwo.','ZA:Dorsal venous arch of foot.')]: test(x+side,y+side)
        for prefix in ['ZA:Great saphenous vein.','ZA:Small saphenous vein.','O3M:Dorsal venous netwo.']:
            for landmark in ['ZA:Tibia.'+side+':Bone','ZA:Fibula.'+side+':Bone','ZA:Talus.'+side+':Bone-1',
                             'ZA:First metatarsal bone.'+side+':Bone','ZA:Fifth metatarsal bone.'+side+':Bone']:
                test(prefix+side,landmark,False)
    assert digest_file(CAT)==sha
    result=dict(schema='lower-limb-native-opportunities/v16',sourceCatalogue=str(CAT),sourceCatalogueSha256=sha,
        censusSha256=digest_file(CENSUS),authoritySha256=digest_file(AUTH),fullApplicableLowerLimbInventory=inventory,
        targetMissingTerms=['A12.3.11.020','A12.3.11.021','A12.3.11.011'],sourceMeshes=sources,
        actualRecipientAndLandmarkTests=tests,candidateControlsPrepared=0,noGeometryWrites=True,noLiveChanges=True,
        noPrivateGeometryRead=True,identificationStatus='Pending source-course and recipient review. No carrier alias accepted.',
        elapsedSeconds=time.monotonic()-began)
    OUT.parent.mkdir(parents=True,exist_ok=True); atomic_json(OUT,json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(output=str(OUT),inventoryTerms=len(inventory),sources=[dict(id=r['partId'],faces=r['topology']['triangles'],components=len(r['components']),loops=len(r['exactCoordinateBoundaryLoops']),bounds_m=r['bounds_m'].tolist(),course=r['sourceCourse'])for r in sources],
        recipientTests=[dict(first=r['first'],second=r['second'],gap_mm=1000*r['closest']['distance_m'],crossings=r.get('actualContactGate',{}).get('crossing_face_pairs'))for r in tests if 'actualContactGate'in r],elapsedSeconds=result['elapsedSeconds']),default=clean))

if __name__=='__main__': main()
