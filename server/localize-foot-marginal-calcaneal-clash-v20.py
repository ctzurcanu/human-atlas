"""Actual inherited lateral source/bone crossings, no physical edits."""
from pathlib import Path
import json, sys, time, importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.review import prepare

CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/lower-limb-native-network-opportunities-v16'
OUT=BASE/'lateral-calcaneal-localization-v20'
FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def main():
    start=time.monotonic();OUT.mkdir(exist_ok=False)
    proposal=json.loads((BASE/'boundary-review-v18/proposal.json').read_bytes())
    audit=json.loads((BASE/'skeletal-and-mouth-review-v19/evidence.json').read_bytes())
    summary=json.loads((BASE/'summary-v17.json').read_bytes())
    sha=digest_file(CAT);assert sha==audit['sourceCatalogueSha256']==proposal['sourceCatalogueSha256']
    a=json.loads(CAT.read_bytes());parts={p['id']:p for p in a['parts']};reader=AtlasReader()
    spec=importlib.util.spec_from_file_location('footbound',ROOT/'server/prepare-marginal-vein-boundaries-v18.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[];reviews=[]
    for row in audit['regions']:
        if row['term']!='A12.3.11.020':continue
        side=row['side'];pid=row['sourcePartId'];boneid='ZA:Calcaneus.'+side+':Bone-5'
        contact=next(q for q in row['completeCandidateSkeletalPairs']if q['partId']==boneid)
        pairs=np.array(contact['sourceFacePairs'],int);assert contact['gate']['crossing_face_pairs']==158==len(pairs)
        m=reader.read_part(a,parts[pid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
        b=reader.read_part(a,parts[boneid],frame=FRAME,catalogue_path=CAT,catalogue_sha=sha)
        sourceids=np.array(row['sourceFaceIds'],int);sourcefaces=sourceids[pairs[:,0]];tris=m.triangles[np.unique(sourcefaces)];center=tris.mean((0,1));low=tris.min((0,1));high=tris.max((0,1))
        s=next(q for q in summary['sourceMeshes']if q['partId']==pid);p=next(q for rr in proposal['rows']for q in rr['marginProposals']if q['sourcePartId']==pid)
        sourcecomponent=np.array(s['components'][0]['sourceFaceIds'],int);rr=module.rings(m,sourcecomponent,np.array(p['fullMouthReplay']['firstSourceRim']['center_m']))
        ix=np.searchsorted(sourcecomponent,sourcefaces);rmin=rr['faceRingMin'][ix];rmax=rr['faceRingMax'][ix]
        entries=[];snaps=[]
        for q,qid,color in [(m,pid,[.9,.12,.63,1.]),(b,boneid,[.72,.72,.68,1.])]:
            t=q.triangles;keep=np.flatnonzero(np.all(t.max(1)>=low-.010,axis=1)&np.all(t.min(1)<=high+.010,axis=1))
            path=OUT/(side+'-'+str(len(entries))+'.npz');Mesh(q.vertices,q.faces[keep],FRAME,q.attributes,dict(q.source,sourceFaceIds=keep.tolist(),diagnosticOnly=True,uncapped=True)).save(path)
            entries.append(dict(mesh=str(path),name=qid+'; actual inherited crossing crop, not repaired',color=color))
            snaps.append(dict(partId=qid,path=str(path),sha256=digest_file(path),sourceFaceIds=keep,sourceGeometrySha256=q.geometry_sha,sourceAttributesSha256=q.attributes_sha))
        normal=np.array([1.,0.,0.]);spec=dict(frame=FRAME,meshes=entries,views=[dict(name=side+' actual SSV foot / calcaneus inherited crossing; lateral',direction=[1 if side=='l'else-1,0,0]),
            dict(name=side+' true sagittal uncapped crossing cut',direction=normal.tolist(),section=dict(normal=normal.tolist(),point=center.tolist(),clear_outer=True))],resolution=1400,wall_seconds=40)
        path=OUT/(side+'-spec.json');atomic_json(path,spec);request=prepare(None,spec)
        reviews.append(dict(side=side,spec=str(path),request=request,snapshots=snaps))
        rows.append(dict(side=side,sourcePartId=pid,bonePartId=boneid,crossingPairs=len(pairs),sourceFacePairs=np.c_[sourcefaces,pairs[:,1]],
            crossedSourceRingIntervals=np.unique(np.c_[rmin,rmax],axis=0),triangleExtent_m=[low,high],center_m=center,
            proposedFootBoundarySourceRing=p['coarseBoundaryRing'],sourceMouthStillUnrepaired=True,
            anatomicalAcceptance=False,allCoordinatesFacesAttributesUnchanged=True))
    assert digest_file(CAT)==sha
    result=dict(schema='literal-lateral-foot-calcaneus-crossing-localization/v20',sourceCatalogueSha256=sha,
        priorSkeletalAuditSha256=digest_file(BASE/'skeletal-and-mouth-review-v19/evidence.json'),rows=rows,reviews=reviews,
        decision='Lateral marginal source interval held. True pre-existing calcaneal crossings cannot be renamed as a correct independent reference.',
        noPhysicalOrCatalogueChange=True,elapsedSeconds=time.monotonic()-start)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],
        rows=[dict(side=q['side'],crossings=q['crossingPairs'],ringIntervals=q['crossedSourceRingIntervals'].tolist())for q in rows],
        requests=[dict(side=q['side'],id=q['request']['request_id'])for q in reviews])))

if __name__=='__main__':main()
