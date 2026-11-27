"""Replay the one rejected source-preserving candidate and queue comparison."""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.review import prepare
from geometry_tools.algorithms.tube_sections import section_metrics

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
OUT=BASE/'candidate-replay-and-comparison-v10'
CANDIDATE=BASE/'one-parasternal-family-candidate-v8'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    p=json.loads((CANDIDATE/'evidence.json').read_bytes()); replay=[]
    for row in p['movingOwners']:
        assert digest_file(row['candidate'])==row['candidateSha256']and digest_file(row['field'])==row['fieldSha256']
        src=Mesh.load(row['source']['path']);candidate=Mesh.load(row['candidate']);d=np.load(row['field']);v=d['physical_vertices'];inv=d['source_to_physical'];rings=d['rings'];c=v[rings].mean(1);R=d['rotation'];target=d['target_centers'];changed=d['changed_stations'];exact=np.einsum('sij,svj->svi',R,v[rings]-c[:,None])+target[:,None]
        encoded=exact.astype('<f4').astype(float);encoded[~changed]=v[rings][~changed];assert np.array_equal(encoded,d['encoded_sections'])
        raw=v.copy();raw[rings]=encoded;raw=raw[inv];assert np.array_equal(raw,candidate.vertices)
        st=np.full(len(v),-1,int)
        for j,r in enumerate(rings):st[r]=j
        expected=src.attributes['normals'].copy()
        for j in np.flatnonzero(changed):
            ids=np.flatnonzero(st[inv]==j);expected[ids]=np.clip(np.rint((src.attributes['normals'][ids].astype(float)/32767@R[j].T)*32767),-32767,32767).astype('<i2')
        assert np.array_equal(expected,candidate.attributes['normals'])and np.array_equal(src.faces,candidate.faces)
        for key in src.attributes:
            if key!='normals':assert np.array_equal(src.attributes[key],candidate.attributes[key])
        before=section_metrics(v[rings][None]);after=section_metrics(encoded[None]);assert np.all(abs(after['area'][0]-before['area'][0])<=d['area_roundoff_bound'])
        protected=np.ones(len(v),bool);protected[rings[changed]]=False;protected=protected[inv]
        assert np.array_equal(src.vertices[protected],candidate.vertices[protected])and np.array_equal(src.attributes['normals'][protected],candidate.attributes['normals'][protected])
        replay.append(dict(partId=row['partId'],geometryAndAllNormalsReplayedExactly=True,sourceFacesUVAndOtherAttributesLiteral=True,unchangedControlsLiteral=True,
            properRotationOrthogonalityMaximum=float(abs(R.transpose(0,2,1)@R-np.eye(3)).max()),minimumProperDeterminant=float(np.linalg.det(R).min()),
            relativeSectionAreaError=float((abs(after['area']-before['area'])/before['area']).max()),sourceGeometrySha256=src.geometry_sha,candidateGeometrySha256=candidate.geometry_sha))
    audit=[json.loads((BASE/('complete-current-context-candidate-v9-'+s)/'evidence.json').read_bytes())for s in ['l','r']]
    assert all(q['completeWithinPairBudgets']and q['all656FixedOrMovingOwnersAccountedFor']for q in audit)
    atlas=json.loads(CAT.read_bytes());parts={q['id']:q for q in atlas['parts']};reader=AtlasReader();fixed=[]
    ids=['ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6','ZA:Pleura','ZA:Transversus thoracis muscle.l','ZA:Transversus thoracis muscle.r','ZA:Anterior segmental artery of left lung']
    ids+=['ZA:Costal cartilage of '+n+' rib.'+s for s in ['l','r']for n in ['first','second','third','fourth','fifth','sixth','seventh']]
    for i,pid in enumerate(ids):
        mesh=reader.read_part(atlas,parts[pid],frame='model:local-ta98:assembled-metres',catalogue_path=CAT,catalogue_sha=EXPECTED);path=OUT/('fixed-'+str(i)+'.npz');mesh.save(path)
        color=[.5,.6,.9,.08]if pid=='ZA:Pleura'else[.85,.45,.16,.35]if'Transversus'in pid else[.75,.06,.85,1.]if'artery'in pid else[.74,.74,.69,.7]
        fixed.append(dict(mesh=str(path),name='FIXED '+pid,color=color))
    center=np.load(BASE/'ita-actual-interface-landmarks-v3/l-ita-source-stations.npz')['centers'][70].tolist();requests=[]
    for kind in ['source','candidate']:
        entries=[]
        for row in p['movingOwners']:
            path=row['source']['path']if kind=='source'else row['candidate'];entries.append(dict(mesh=path,name=kind.upper()+' '+row['partId'],color=[.92,.12,.15,1.]if kind=='source'else[.05,.8,.95,1.]))
        spec=dict(frame='model:local-ta98:assembled-metres',meshes=entries+fixed,
            views=[dict(name=kind.upper()+' whole arterial family with unchanged skeletal and soft-tissue context',direction=[.7,.15,-1.]),
                dict(name=kind.upper()+' true uncapped axial cut at original ITA station70',direction=[0,1,0],section=dict(normal=[0,1,0],point=center,clear_outer=True))],resolution=1400,wall_seconds=60)
        path=OUT/(kind+'-review-spec.json');atomic_json(path,spec);request=prepare(None,spec);requests.append(dict(kind=kind,spec=str(path),request=request))
    result=dict(schema='one-ita-candidate-independent-replay-comparison/v10',sourceCatalogueSha256=EXPECTED,
        candidateEvidenceSha256=digest_file(CANDIDATE/'evidence.json'),allSixSourceArrayReplays=replay,
        leftAuditSha256=digest_file(BASE/'complete-current-context-candidate-v9-l/evidence.json'),rightAuditSha256=digest_file(BASE/'complete-current-context-candidate-v9-r/evidence.json'),
        sameFixedContextAndViewDirectionsInBothRequests=True,requests=requests,disposition='REJECTED: new arterial, pleural and transversus contacts. Not installed.',
        nextRequiredConstraint='A two-sided native parasternal compartment bounded anteriorly by posterior cartilage/sternal surfaces and posteriorly by anterior transversus/pleural surfaces, with actual complete-ring source calibre and all current-context nonregression.',
        noSecondOffsetTrial=True,noAtlasChanges=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',result);assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],replayedOwners=len(replay),requests=[dict(kind=q['kind'],request_id=q['request']['request_id'])for q in requests])))

if __name__=='__main__':main()
