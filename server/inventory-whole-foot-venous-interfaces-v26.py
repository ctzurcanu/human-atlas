"""Read-only complete native foot venous source interfaces and tissue ends.

Boundary proximity is retained as diagnostic evidence, never promoted to an
annulus. Source casts are outer surfaces, not finite walls or blood volumes.
"""
from pathlib import Path
import json,sys,time,itertools
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
OUT=BASE/'whole-foot-native-interfaces-v26'
FRAME='model:local-ta98:assembled-metres'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def area(v):
    q=v-v.mean(0);n=np.cross(q,np.roll(q,-1,axis=0)).sum(0)/2
    return np.linalg.norm(n),n/np.linalg.norm(n)

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    a=json.loads(CAT.read_bytes());ps={p['id']:p for p in a['parts']};reader=AtlasReader();cache={};trees={};sources=[];loops=[];pairs=[];ends=[]
    # All literal native venous sources reaching the foot/low ankle, plus
    # the actual proximal recipients for the two superficial trunks.
    pids=[p['id']for p in a['parts']if p['id'].startswith('ZA:')and p.get('system')=='venous'and p.get('bounds')and p['bounds'][0][1]<.14 and p['bounds'][1][1]>.01]
    pids+=['ZA:Femoral vein.'+s for s in ['l','r']]+['ZA:Popliteal vein.'+s for s in ['l','r']]
    def get(pid):
        if pid not in cache:cache[pid]=reader.read_part(a,ps[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=EXPECTED)
        return cache[pid]
    def tree(pid):
        if pid not in trees:trees[pid]=BVH(get(pid))
        return trees[pid]
    for pid in pids:
        m=get(pid);v,first,inv=np.unique(m.vertices,axis=0,return_index=True,return_inverse=True);f=inv[m.faces];ll=boundary_loops(f);path=OUT/(identity(pid)[:20]+'.npz');m.save(path)
        row=dict(id=pid,sourceId=ps[pid]['sourceId'],mesh=str(path),meshSha256=digest_file(path),recordSha256=identity(ps[pid]),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,bankSha256=m.source['bank_sha256'],
            bounds_m=m.bounds,topology=topology(m),outerSurfaceOnly=True,finiteWallAndPatencyUnestablished=True,boundaryLoopIds=[])
        for j,l in enumerate(ll):
            q=v[l];ar,no=area(q);compact=first[l];normal=m.attributes['normals'][compact].astype(float);normal/=np.linalg.norm(normal,axis=1)[:,None]
            loop=dict(id=pid+'#mouth-'+str(j),partId=pid,componentBoundaryOrdinal=j,physicalVertexIds=l,compactVertexIds=compact,points_m=q,center_m=q.mean(0),area_m2=ar,directedNormal=no,
                packedNormals=m.attributes['normals'][compact],unitNormals=normal,diameter_m=float(np.linalg.norm(q[:,None]-q[None],axis=2).max()),
                trueOrientedSourceBoundary=True,notFiniteAnnulus=True)
            loops.append(loop);row['boundaryLoopIds'].append(loop['id'])
        sources.append(row)
    # Full12-point correspondences for every spatially adjacent original
    # mouth pair, including every source component within a grouped item.
    matches=[]
    for x,y in itertools.combinations(loops,2):
        if x['partId'].endswith('.l')!=y['partId'].endswith('.l'):continue
        centerGap=float(np.linalg.norm(x['center_m']-y['center_m']))
        if centerGap>.015 or len(x['points_m'])!=len(y['points_m']):continue
        choices=[]
        for shift in range(len(y['points_m'])):
            # Opposed source boundary winding is required for a real join.
            q=np.roll(y['points_m'][::-1],shift,axis=0);dist=np.linalg.norm(q-x['points_m'],axis=1)
            choices.append((float(dist.max()),float(np.sqrt((dist*dist).mean())),shift,q))
        maximum,rms,shift,target=min(choices,key=lambda q:(q[0],q[1],q[2]))
        nn=np.roll(y['unitNormals'][::-1],shift,axis=0);cos=np.sum(x['unitNormals']*nn,axis=1)
        matches.append(dict(first=x['id'],second=y['id'],centerGap_m=centerGap,opposedCyclicShift=shift,
            completeMouthMaximumCoordinateDifference_m=maximum,rms_m=rms,relativeAreaDifference=float(abs(x['area_m2']/y['area_m2']-1)),
            areaRatio=float(x['area_m2']/y['area_m2']),maximumNormalChartAngle_deg=float(np.degrees(np.arccos(np.clip(cos,-1,1))).max()),
            exactSharedCoordinates=maximum==0.,sharedNormalWords=False,
            directedPlaneNormalAngle_deg=float(np.degrees(np.arccos(np.clip(np.dot(x['directedNormal'],-y['directedNormal']),-1,1)))),
            relation='Diagnostic adjacent original outer-source mouths; this proximity does not establish a joined anatomical interface.'))
    for x,y in itertools.combinations(pids,2):
        if x.endswith('.l')!=y.endswith('.l'):continue
        m,n=get(x),get(y);bb,cc=m.bounds,n.bounds
        if np.any(bb[1]<cc[0])or np.any(bb[0]>cc[1]):continue
        detail,gate=pair_contacts(m,n,0.,first_bvh=tree(x),second_bvh=tree(y),max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget']
        pairs.append(dict(first=x,second=y,gate=gate,pairs=detail['pairs'],relations=detail['relation']))
    for s in ['l','r']:
        for trunk,recipient in [('ZA:Great saphenous vein.','ZA:Femoral vein.'),('ZA:Small saphenous vein.','ZA:Popliteal vein.')]:
            pid=trunk+s;rid=recipient+s;own=[l for l in loops if l['partId']==pid];proximal=max(own,key=lambda q:q['center_m'][1]);q=proximal['points_m'];ring=Mesh(q,np.c_[np.arange(len(q)),np.roll(np.arange(len(q)),-1),np.arange(len(q))],FRAME)
            closest=BVH(ring).closest_mesh(tree(rid));rc=[l for l in loops if l['partId']==rid];near=min(rc,key=lambda z:np.linalg.norm(z['center_m']-proximal['center_m']))
            ends.append(dict(trunkPartId=pid,proximalBoundaryLoopId=proximal['id'],recipientPartId=rid,recipientClosestActualTriangle=closest,
                nearestRecipientSourceMouth=near['id'],recipientMouthCenterDistance_m=float(np.linalg.norm(near['center_m']-proximal['center_m'])),
                proposedFixedReference='Actual highest source boundary, not a source bounding-box centre or invented closed container.',
                literalMatchedAnnulusEstablished=any(z['exactSharedCoordinates']and {z['first'],z['second']}=={proximal['id'],near['id']}for z in matches),
                completeTissueWallOrFlowJoinEstablished=False))
    assert digest_file(CAT)==EXPECTED
    result=dict(schema='whole-native-foot-venous-source-interface-inventory/v26',sourceCatalogueSha256=EXPECTED,sources=sources,sourceBoundaryLoops=loops,
        nearbyCompleteSourceMouthPairs=matches,actualSourceSurfaceContactPairs=pairs,actualProximalTissueEnds=ends,
        sourceFamilyIsNotAutomaticallyOneConnectedTree=True,noPhysicalSourceChanges=True,noNewSourceParentContainers=True,
        defaultRuntimeSourceQualification='Named licensed ZA native source meshes are measured literally; no private or alternative mesh is substituted.',
        anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),sources=len(sources),sourceMouths=len(loops),nearbyMouthPairs=len(matches),exactSharedMouths=sum(q['exactSharedCoordinates']for q in matches),
        fullSurfacePairs=len(pairs),elapsedSeconds=result['elapsedSeconds'],proximalEnds=[dict(trunk=r['trunkPartId'],recipient=r['recipientPartId'],surfaceGap_mm=1000*r['recipientClosestActualTriangle']['distance_m'],recipientMouthCenterGap_mm=1000*r['recipientMouthCenterDistance_m'],literalAnnulus=r['literalMatchedAnnulusEstablished'])for r in ends])))

if __name__=='__main__':main()
