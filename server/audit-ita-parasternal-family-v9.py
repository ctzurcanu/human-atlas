"""Whole changed-support, actual-runtime nonregression for one ITA candidate."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts,self_contact_summary
from geometry_tools.spatial import BVH

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
CANDIDATE=BASE/'one-parasternal-family-candidate-v8'

def intersects(a,b):return bool(np.all(a[1]>=b[0])and np.all(b[1]>=a[0]))

def contacts(a,b,at,bt):
    if not intersects(a.bounds,b.bounds):return dict(crossing_face_pairs=0,coplanar_overlap_face_pairs=0,touching_face_pairs=0,complete_within_pair_budget=True,broadPhaseDisjoint=True),[]
    detail,gate=pair_contacts(a,b,0.,first_bvh=at,second_bvh=bt,max_pairs=250000,max_tests=2000000)
    wrong=np.isin(detail['relation'],['crossing','coplanar-overlap'])
    return gate,detail['pairs'][wrong].tolist()

def ap_enclosure(mesh,points):
    # Actual AP source-face rays only. This is explicit open-surface slice
    # evidence, not a closed-solid containment claim or a blood lumen test.
    tri=mesh.triangles;lo=points[:,:2].min(0);hi=points[:,:2].max(0)
    ids=np.flatnonzero(np.all(tri[:,:,:2].max(1)>=lo,axis=1)&np.all(tri[:,:,:2].min(1)<=hi,axis=1));tt=tri[ids]
    a=tt[:,1,:2]-tt[:,0,:2];b=tt[:,2,:2]-tt[:,0,:2];det=a[:,0]*b[:,1]-a[:,1]*b[:,0];valid=abs(det)>1e-20
    counts=[];examples=[];ambiguous=[]
    for start in range(0,len(points),256):
        q=points[start:start+256];d=q[:,None,:2]-tt[None,:,0,:2];w1=np.zeros(d.shape[:2]);w2=w1.copy()
        np.divide(d[:,:,0]*b[None,:,1]-d[:,:,1]*b[None,:,0],det[None],out=w1,where=valid[None]);np.divide(a[None,:,0]*d[:,:,1]-a[None,:,1]*d[:,:,0],det[None],out=w2,where=valid[None])
        good=valid[None]&(w1>=0)&(w2>=0)&(w1+w2<=1);z=tt[None,:,0,2]+w1*(tt[None,:,1,2]-tt[None,:,0,2])+w2*(tt[None,:,2,2]-tt[None,:,0,2]);back=good&(z<q[:,None,2]);front=good&(z>q[:,None,2]);inside=back.any(1)&front.any(1);counts.extend((start+np.flatnonzero(inside)).tolist())
        for k in np.flatnonzero(inside)[:max(0,12-len(examples))]:
            hit=np.flatnonzero(good[k]);examples.append(dict(pointIndex=int(start+k),point_m=q[k].tolist(),sourceFaces=ids[hit].tolist(),hitZ_m=z[k,hit].tolist(),barycentric=np.c_[1-w1[k,hit]-w2[k,hit],w1[k,hit],w2[k,hit]].tolist()))
    return dict(queryPoints=len(points),actualTrianglesUsed=len(ids),apBracketedPoints=len(counts),pointIds=counts,examples=examples,
        unmodeledClosedSolidContainmentNotInferred=True,parallelTriangleProjectionsIgnored=int((~valid).sum()))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--side',choices=['l','r'],required=True);args=parser.parse_args();side=args.side
    began=time.monotonic();OUT=BASE/('complete-current-context-candidate-v9-'+side);OUT.mkdir(exist_ok=False)
    proposal=json.loads((CANDIDATE/'evidence.json').read_bytes());runtime=json.loads((BASE/'ita-distal-closure-context-v6/evidence.json').read_bytes());v3=json.loads((BASE/'ita-actual-interface-landmarks-v3/evidence.json').read_bytes())
    moving={q['partId']:q for q in proposal['movingOwners']};selected={k:q for k,q in moving.items()if k.endswith('.'+side)};supports={q['partId']:q for q in proposal['changedSupport']};meshes={};trees={};records=[]
    for pid,row in selected.items():
        assert digest_file(row['candidate'])==row['candidateSha256']
        source=Mesh.load(row['source']['path']);candidate=Mesh.load(row['candidate']);faces=np.array(supports[pid]['sourceFaceIds'],int)
        src=Mesh(source.vertices,source.faces[faces],source.frame);dst=Mesh(candidate.vertices,candidate.faces[faces],candidate.frame)
        meshes[pid]=(source,candidate,src,dst,faces);trees[pid]=(BVH(src),BVH(dst))
    count=0
    for context in runtime['rows']:
        pid=context['partId'];count+=1
        if pid in moving:
            records.append(dict(partId=pid,role='moving-family-owner',handledByCompleteMutualFamilyChecks=True));continue
        if not context['retainedTriangles']:
            records.append(dict(partId=pid,role='fixed',candidateAndSourceChangedSupportsDisjointFromWholeSourceTriangleAABBs=True));continue
        assert digest_file(context['mesh'])==context['sha256']; obstacle=Mesh.load(context['mesh']);relevant=[k for k,(whole,cm,s,d,f)in meshes.items()if intersects(s.bounds,obstacle.bounds)or intersects(d.bounds,obstacle.bounds)]
        if not relevant:
            records.append(dict(partId=pid,role='fixed',candidateAndSourceChangedSupportAABBsDisjoint=True));continue
        tree=BVH(obstacle);checks=[]
        for key in relevant:
            whole,cm,src,dst,faces=meshes[key];before,pbefore=contacts(src,obstacle,trees[key][0],tree);after,pafter=contacts(dst,obstacle,trees[key][1],tree)
            checks.append(dict(movingPartId=key,sourceGate=before,candidateGate=after,sourceDefectSourceFacePairs=[[int(faces[a]),int(context['sourceFaceIds'][b])]for a,b in pbefore],
                candidateDefectSourceFacePairs=[[int(faces[a]),int(context['sourceFaceIds'][b])]for a,b in pafter],
                crossingsNonRegressing=after['crossing_face_pairs']<=before['crossing_face_pairs'],
                overlapsNonRegressing=after['coplanar_overlap_face_pairs']<=before['coplanar_overlap_face_pairs']))
        records.append(dict(partId=pid,role='fixed',system=context['system'],checks=checks))
        if count%100==0:atomic_json(OUT/('progress-'+str(count)+'.json'),dict(processed=count,elapsedSeconds=time.monotonic()-began))
    mutual=[];selfChecks=[];keys=list(selected)
    for key in keys:
        src,candidate,*_=meshes[key];before,_=self_contact_summary(src,max_pairs=250000,max_tests=2000000);after=selected[key]['wholeSelfContact']
        selfChecks.append(dict(partId=key,sourceGate=before,candidateGate=after,crossingsNonRegressing=after['crossing_face_pairs']<=before['crossing_face_pairs']))
    for i,key in enumerate(keys):
        for other in keys[i+1:]:
            s,c,*_=meshes[key];ss,cc,*_=meshes[other];before,pb=contacts(s,ss,BVH(s),BVH(ss));after,pa=contacts(c,cc,BVH(c),BVH(cc))
            mutual.append(dict(first=key,second=other,sourceGate=before,candidateGate=after,sourcePairs=pb,candidatePairs=pa,crossingsNonRegressing=after['crossing_face_pairs']<=before['crossing_face_pairs']))
    # Direct full-source AP tests use v3's whole native meshes, not the
    # conservative clipped context views used for surface broad phase.
    boneIds=['ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6']+['ZA:Costal cartilage of '+n+' rib.'+side for n in ['first','second','third','fourth','fifth','sixth','seventh']]
    skeletal=[]
    for key,(src,dst,ss,dd,faces)in meshes.items():
        old=np.unique(np.r_[src.vertices[np.unique(src.faces[faces])],src.triangles[faces].mean(1)],axis=0);new=np.unique(np.r_[dst.vertices[np.unique(dst.faces[faces])],dst.triangles[faces].mean(1)],axis=0)
        for bid in boneIds:
            obstacle=Mesh.load(v3['sourceBindings'][bid]['path']);before=ap_enclosure(obstacle,old);after=ap_enclosure(obstacle,new)
            skeletal.append(dict(movingPartId=key,fixedPartId=bid,source=before,candidate=after,sourceAndCandidateWholeMeshGeometrySha256=obstacle.geometry_sha))
    complete=all(z[g]['complete_within_pair_budget']for q in records for z in q.get('checks',[])for g in ['sourceGate','candidateGate'])and all(q[g]['complete_within_pair_budget']for q in mutual for g in ['sourceGate','candidateGate'])
    regressions=[dict(partId=q['partId'],moving=z['movingPartId'],before=z['sourceGate']['crossing_face_pairs'],after=z['candidateGate']['crossing_face_pairs'],oldOverlap=z['sourceGate']['coplanar_overlap_face_pairs'],newOverlap=z['candidateGate']['coplanar_overlap_face_pairs'])for q in records for z in q.get('checks',[])if not(z['crossingsNonRegressing']and z['overlapsNonRegressing'])]
    result=dict(schema='ita-one-candidate-all656-context-nonregression/v9',side=side,candidateEvidenceSha256=digest_file(CANDIDATE/'evidence.json'),runtimeContextSha256=digest_file(BASE/'ita-distal-closure-context-v6/evidence.json'),
        all656FixedOrMovingOwnersAccountedFor=len(records)==656,completeWithinPairBudgets=complete,rows=records,wholeSelfChecks=selfChecks,wholeMutualFamilyChecks=mutual,fullSourceSkeletalAPEnclosure=skeletal,
        regressions=regressions,allOtherAnatomyFixed=True,acceptedForInstallation=False,finiteBloodLumenNotEstablished=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',result)
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],complete=complete,regressions=regressions,
        candidateAPBracketed=[dict(owner=q['movingPartId'],bone=q['fixedPartId'],before=q['source']['apBracketedPoints'],after=q['candidate']['apBracketedPoints'])for q in skeletal if q['source']['apBracketedPoints']or q['candidate']['apBracketedPoints']])))

if __name__=='__main__':main()
