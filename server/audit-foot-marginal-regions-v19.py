"""Literal foot-vein source regions: complete skeletal screen and mouth cuts.

No anatomical coordinate, face, UV or normal is changed. Diagnostic source
triangle crops are kept open. This audit does not establish a blood lumen.
"""
from pathlib import Path
import json, sys, time
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file, identity
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology, winding_contains
from geometry_tools.contact import pair_contacts, patch_components, self_contact_summary
from geometry_tools.spatial import BVH
from geometry_tools.review import prepare

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
PROPOSAL=BASE/'boundary-review-v18/proposal.json'
OUT=BASE/'skeletal-and-mouth-review-v19'
FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    start=time.monotonic();assert not OUT.exists();OUT.mkdir()
    proposal=json.loads(PROPOSAL.read_bytes());sha=digest_file(CAT)
    assert sha==proposal['sourceCatalogueSha256']
    a=json.loads(CAT.read_bytes());parts={p['id']:p for p in a['parts']}
    reader=AtlasReader();cache={};bindings={};trees={};reviews=[];rows=[]
    def get(pid):
        if pid not in cache:
            cache[pid]=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
            bindings[pid]=dict(recordSha256=identity(parts[pid]),geometrySha256=cache[pid].geometry_sha,
                attributesSha256=cache[pid].attributes_sha,bankSha256=cache[pid].source['bank_sha256'])
        return cache[pid]
    def tree(pid):
        if pid not in trees:trees[pid]=BVH(get(pid))
        return trees[pid]
    for side in proposal['rows']:
        for r in side['marginProposals']:
            pid=r['sourcePartId'];source=get(pid);faceids=np.array(r['proposedLiteralFootSourceFaceIds'],int)
            region=Mesh(source.vertices,source.faces[faceids],FRAME,source.attributes,
                dict(source.source,sourceFaceIds=faceids.tolist(),literalSubset=True))
            rbvh=BVH(region);bounds=region.bounds
            skeletal=[p for p in a['parts'] if p.get('system')=='skeletal' and p.get('bounds')]
            candidate=[];excluded=[]
            for p in skeletal:
                box=np.array(p['bounds'],float)
                # Recorded catalogue bounds are broad-phase retrieval only;
                # actual source arrays govern every surviving pair.
                if np.any(box[1]<bounds[0])or np.any(box[0]>bounds[1]):excluded.append(p['id']);continue
                b=get(p['id']);actual=b.bounds
                if np.any(actual[1]<bounds[0])or np.any(actual[0]>bounds[1]):excluded.append(p['id']);continue
                detail,gate=pair_contacts(region,b,0.,first_bvh=rbvh,second_bvh=tree(p['id']),max_pairs=250000,max_tests=2000000)
                assert gate['complete_within_pair_budget']
                top=topology(b);contained=None
                if top['closed_oriented_manifold'] and gate['crossing_face_pairs']==0:
                    # One literal vertex per connected surface component,
                    # combined with complete triangle-crossing exclusion.
                    comps=patch_components(region,np.arange(len(region.faces)))
                    probes=np.array([region.vertices[region.faces[q[0],0]]for q in comps])
                    contained=winding_contains(b,probes)
                closest=rbvh.closest_mesh(tree(p['id']))
                candidate.append(dict(partId=p['id'],gate=gate,sourceFacePairs=detail['pairs'],
                    closestDistance_m=closest['distance_m'],boneTopology=top,closedBoneComponentContainment=contained,
                    defaultHidden=bool(p.get('defaultHidden')),referenceOrDerivativeNotOmitted=True))
            selfgate,_=self_contact_summary(region,max_pairs=250000,max_tests=2000000)
            assert selfgate['complete_within_pair_budget']
            # Focus the actual original vessel/arch mouths, retaining all
            # triangles whose boxes meet a 12 mm diagnostic neighborhood.
            mouth=r['fullMouthReplay'];center=np.array(mouth['firstSourceRim']['center_m'])
            entries=[];snapshots=[]
            for qid,label,color in [(pid,'original '+r['kind']+' course',[.05,.8,.95,1.]),
                    (mouth['recipientPartId'],'original dorsal venous arch',[.55,.22,.88,1.])]:
                q=get(qid);t=q.triangles;keep=np.flatnonzero(np.all(t.max(1)>=center-.012,axis=1)&np.all(t.min(1)<=center+.012,axis=1))
                assert len(keep)
                view=Mesh(q.vertices,q.faces[keep],FRAME,q.attributes,dict(q.source,sourceFaceIds=keep.tolist(),diagnosticTriangleCrop=True,uncapped=True))
                path=OUT/(side['side']+'-'+r['kind']+'-'+str(len(entries))+'.npz');view.save(path)
                entries.append(dict(mesh=str(path),name=label+'; literal uncapped diagnostic crop',color=color))
                snapshots.append(dict(partId=qid,path=str(path),sha256=digest_file(path),sourceFaceIds=keep,
                    sourceGeometrySha256=q.geometry_sha,sourceAttributesSha256=q.attributes_sha))
            points=np.array(mouth['firstSourceRim']['points_m']);_,_,vh=np.linalg.svd(points-center,full_matrices=False)
            tangent=vh[-1];ref=np.eye(3)[np.argmin(abs(tangent))];n1=np.cross(tangent,ref);n1/=np.linalg.norm(n1);n2=np.cross(tangent,n1)
            spec=dict(frame=FRAME,meshes=entries,views=[dict(name=side['side']+' '+r['kind']+' actual source mouths; end-on',direction=tangent.tolist()),
                dict(name=side['side']+' '+r['kind']+' actual source mouths; first uncapped longitudinal cut',direction=n1.tolist(),section=dict(normal=n1.tolist(),point=center.tolist(),clear_outer=True)),
                dict(name=side['side']+' '+r['kind']+' actual source mouths; orthogonal uncapped longitudinal cut',direction=n2.tolist(),section=dict(normal=n2.tolist(),point=center.tolist(),clear_outer=True))],resolution=1400,wall_seconds=45)
            specpath=OUT/(side['side']+'-'+r['kind']+'-spec.json');atomic_json(specpath,spec);request=prepare(None,spec)
            reviews.append(dict(side=side['side'],kind=r['kind'],spec=str(specpath),request=request,snapshots=snapshots))
            rows.append(dict(side=side['side'],term=r['term'],sourcePartId=pid,sourceFaceIds=faceids,
                sourceRegionBounds_m=bounds,wholeSourceFaces=len(source.faces),regionFaces=len(faceids),
                skeletalSourcesConsidered=len(skeletal),skeletalBroadPhaseExcludedIds=excluded,
                completeCandidateSkeletalPairs=candidate,regionTopology=topology(region),regionSelfContact=selfgate,
                bloodLumenEstablished=False,nativeMouthCoordinateMismatch_m=mouth['maxVertexDifference_m']))
    assert digest_file(CAT)==sha
    result=dict(schema='literal-foot-marginal-source-skeletal-mouth-audit/v19',sourceCatalogueSha256=sha,
        proposalSha256=digest_file(PROPOSAL),sourceBindings=bindings,regions=rows,reviews=reviews,
        noCoordinateOrAttributeChange=True,noAnatomicalFacesChanged=True,noCapAdded=True,noInstalledChanges=True,
        elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],
        regions=[dict(term=x['term'],side=x['side'],skeletalPairs=len(x['completeCandidateSkeletalPairs']),
            crossings=sum(q['gate']['crossing_face_pairs']for q in x['completeCandidateSkeletalPairs']),
            selfCrossings=x['regionSelfContact']['crossing_face_pairs'])for x in rows],
        reviews=[dict(side=q['side'],kind=q['kind'],id=q['request']['request_id'])for q in reviews])))

if __name__=='__main__':main()
