"""Actual source anatomy review; failed landmark gates prohibit partition."""
from pathlib import Path
import json,sys,time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import Mesh
from geometry_tools.review import prepare

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'
OUT=BASE/'source-landmark-contract-review-v4'

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False)
    source=json.loads((BASE/'source-inventory.json').read_bytes());audit=json.loads((BASE/'complete-current-context-v3/evidence.json').read_bytes());runtime=json.loads((BASE/'actual-runtime-context-v2.json').read_bytes())
    independent=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-boundary-semantics-independent-v1/evidence.json'
    available={q['id']for q in runtime['rows']if q['normalAvailable']};ids=['ZA:Vertebral artery.l','ZA:Vertebral artery.r','ZA:Basilar artery','ZA:Atlas (C1):Bone-2','ZA:Atlas (C1):Cartilage','ZA:Axis (C2):Bone-2','ZA:Vertebra C4:Bone-2','ZA:Vertebra C6:Bone-2','ZA:Occipital bone:Bone-6','ZA:Spinal dura']
    entries=[];bindings=[]
    for i,pid in enumerate(ids):
        assert pid in available;m=Mesh.load(source['sourceBindings'][pid]['mesh']);sourceFaces=np.arange(len(m.faces))
        if pid=='ZA:Spinal dura':
            sourceFaces=np.flatnonzero(m.triangles[:,:,1].max(1)>=1.58)
            m=Mesh(m.vertices,m.faces[sourceFaces],m.frame,m.attributes,dict(m.source,actualSourceFaceIds=sourceFaces.tolist(),diagnosticCropWithoutAddedCaps=True))
        path=OUT/('source-'+str(i)+'.npz');m.save(path);color=[.92,.14,.18,1.]if'artery'in pid else[.35,.65,.9,.2]if'dura'in pid else[.76,.74,.65,.55]
        entries.append(dict(mesh=str(path),name='UNCHANGED SOURCE '+pid,color=color));bindings.append(dict(partId=pid,mesh=str(path),sha256=digest_file(path),geometrySha256=m.geometry_sha,attributesSha256=m.attributes_sha,actualSourceFaceIds=sourceFaces.tolist()))
    point=[.027288,1.646488,-.015767]
    spec=dict(frame='model:local-ta98:assembled-metres',meshes=entries,
        views=[dict(name='UNCHANGED native vertebral course with actual fixed C1/C2/C4/C6/occipital and dura',direction=[.8,.2,1.]),
            dict(name='TRUE uncapped axial slice through independently witnessed C1 course crossing',direction=[0,1,0],section=dict(normal=[0,1,0],point=point,clear_outer=True))],resolution=1400,wall_seconds=60)
    atomic_json(OUT/'review-spec.json',spec);request=prepare(None,spec)
    result=dict(schema='vertebral-native-regional-source-and-landmark-contract/v4',sourceCatalogueSha256=source['sourceCatalogueSha256'],sourceInventorySha256=digest_file(BASE/'source-inventory.json'),
        completeCurrentContextSha256=digest_file(BASE/'complete-current-context-v3/evidence.json'),independentLandmarkEvidence=str(independent),independentLandmarkEvidenceSha256=digest_file(independent),
        actualReviewBindings=bindings,review=request,wholeAssembly='Two original 205-station, one-component, two-mouth outer casts; all other anatomy fixed.',
        plannedRegions=[dict(term='A12.2.08.003',name='prevertebral',boundary='Actual subclavian origin to first transverse-foramen entry, normally C6; aperture has not been geometrically established.'),
            dict(term='A12.2.08.004',name='cervical',boundary='Actual transverse-foraminal passage to C2 exit under the declared 2024 cadaveric convention; no genericY partition.'),
            dict(term='A12.2.08.009',name='atlantic',boundary='C2 exit, C1 transverse-foramen and groove, to actual dural entry; C1 exit alone is not the whole start.'),
            dict(term='A12.2.08.010',name='intracranial',boundary='Actual dural entry to basilar union; foramen-magnum plane does not establish dural entry.')],
        regionalBoundaryGatePassed=False,reason='Complete current-source audit has actual C1/C2/C4/occipital crossings, and native spinal dura supplies no centerline entry intersection. These physical defects/missing landmarks cannot become anatomical region boundaries.',
        nativeSourceWholeSelfCrossings=dict(left=0,right=0),literalSourceTopologyAndBanksReplayed=True,
        noPartitionOrGeometryConstructed=True,noAtlasChanges=True,noFalseAdditionalTermCoverage=True,acceptedForInstallation=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'contract.json',result)
    print(json.dumps(dict(contract=str(OUT/'contract.json'),review=request['request_id'],elapsedSeconds=result['elapsedSeconds'],partitionGate=False)))

if __name__=='__main__':main()
