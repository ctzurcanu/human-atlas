"""Restore actual source perforator contours after common GSV course transport.

Explicitly owned four cylindrical strips per source family; existing vertices,
faces and source ancestry survive. This builds a diagnostic, never installs it.
"""
import hashlib
import importlib.util
import argparse
import json
import shutil
import sys
import time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
_spec=importlib.util.spec_from_file_location('perforator_ring_topology',ROOT/'server/audit-biodigital-perforator-rings.py')
_module=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(_module)
rings=_module.rings
from geometry_tools.algorithms.tube_sections import redirect_local_sections, section_metrics, bend_diagnostics
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs


def main():
    started=time.perf_counter();w=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'
    previous=w/'saphenous-transport-v2/candidate-v1';base=w/'lower-limb-joints-v4/candidate-v1'
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',default='saphenous-sections-v3')
    args=parser.parse_args();assert Path(args.stage).name==args.stage and args.stage not in ('.','..')
    stage=w/args.stage;stage.mkdir(exist_ok=False);out=stage/'candidate-v1';out.mkdir()
    (out/'checkpoint.json').write_text(json.dumps({'status':'source-section-transport-preparing'})+'\n')
    assert json.loads((previous/'checkpoint.json').read_text())['status']=='candidate-numerical-gates-passed-anatomy-unverified'
    for p in previous.iterdir():
        if p.suffix=='.npz' or p.name=='instance-ownership.json':shutil.copyfile(p,out/p.name)
    contract=json.loads((w/'saphenous-transport-v2/contract.json').read_text())
    contract['stagePurpose']='Common GSV finite transport followed by explicit source perforator-strip rigid contour transport. Shared boundaries/mouth contacts remain subject to audit.'
    contract['source_section_transport']={'baseline':'lower-limb-joints-v4/candidate-v1','sourceFamilies':['Left femoral canal perforator veins','Right femoral canal perforator veins'],'sourceTopology':'four exact open cylindrical strips each; two literal mouths and every cyclic contour recovered by edge adjacency','method':'proper minimum rotations between source and common-transport centreline tangents for interior contours; both common-mapped mouths retained exactly','commonMappedCentreLinesPreserved':True,'wholeBodyAcceptance':False}
    (stage/'contract.json').write_text(json.dumps(contract,indent=2)+'\n');(out/'input-contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    ownership=json.loads((base/'instance-ownership.json').read_text())['instances']
    with np.load(previous/'venous.npz')as z:arrays={k:z[k]for k in z.files}
    with np.load(base/'venous.npz')as z:source=z['vertices'];source_faces=z['faces']
    assert np.array_equal(arrays['faces'],source_faces)
    vertices=arrays['vertices'].copy();original=vertices.copy();transports=[]
    rotations=np.broadcast_to(np.eye(3),(len(vertices),3,3)).copy()
    for part in ownership:
        if part['name'] not in contract['source_section_transport']['sourceFamilies']:continue
        start=part['candidateVertexStart'];stop=start+part['candidateVertexCount'];fs=part['candidateFaceStart'];fe=fs+part['candidateFaceCount']
        f=source_faces[fs:fe]-start;groups=rings(source[start:stop],f)
        for ci,group in enumerate(groups):
            ids=start+np.asarray(group)
            p=source[ids][None];current=vertices[ids][None]
            if np.max(np.linalg.norm(p-current,axis=3))<1e-9:continue
            centers=current[0].mean(1)
            new,rotation=redirect_local_sections(p,centers,dtype=np.float64)
            # A mouth is the adjoining interface. Its existing common map
            # must match the mother surface; it must not be reoriented by
            # an interior centreline tangent estimate.
            new[:,[0,-1]]=current[:,[0,-1]]
            a,b=section_metrics(p),section_metrics(new)
            pair_error=float(np.max(np.abs(a['shape_distances'][:,1:-1]-b['shape_distances'][:,1:-1])))
            assert pair_error<1e-12
            vertices[ids]=new[0]
            rotations[ids]=rotation[:,None]
            bend=bend_diagnostics(centers,new)
            transports.append(dict(sourceId=part['sourceInstanceId'],component=ci,sourceVertexIds=(ids-start).tolist(),
                maximumRingPairDistanceError_m=pair_error,minimumAreaRatio=float(np.min(b['area']/a['area'])),
                minimumInteriorAreaRatio=float(np.min(b['area'][:,1:-1]/a['area'][:,1:-1])),
                commonMappedMouthCoordinatesExactlyRetained=bool(np.array_equal(new[:,[0,-1]],current[:,[0,-1]])),
                maximumAreaRatio=float(np.max(b['area']/a['area'])),minimumRotationDeterminant=float(np.linalg.det(rotation).min()),
                centersBefore_m=p[0].mean(1).tolist(),commonMappedCenters_m=centers.tolist(),
                estimatedMinimumSweepJacobian=bend['minimum_sampled_sweep_jacobian'],
                estimatedTightBendStations=bend['failed_stations'].tolist()))
    moved=np.linalg.norm(vertices-original,axis=1)>1e-9;faces=arrays['faces'];changed=moved[faces].any(1)
    f=faces[changed];tri=vertices[f];before=original[f]
    cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1)/2
    assert np.all(area>1e-18)
    # Recompute owned smooth vertex normals from literal repaired triangles.
    all_cross=np.cross(vertices[faces[:,1]]-vertices[faces[:,0]],vertices[faces[:,2]]-vertices[faces[:,0]])
    raw=np.zeros_like(vertices)
    for k in range(3):np.add.at(raw,faces[:,k],all_cross)
    norm=np.linalg.norm(raw,axis=1);unit=np.divide(raw,norm[:,None],out=np.zeros_like(raw),where=norm[:,None]>0)
    key='attribute_normals';assert key in arrays
    affected_normals=np.unique(faces[changed])
    arrays[key]=arrays[key].copy();arrays[key][affected_normals]=unit[affected_normals];arrays['vertices']=vertices
    np.savez(out/'venous.npz',**arrays)
    test=Mesh(tri.reshape(-1,3),np.arange(len(tri)*3).reshape(-1,3),'active-male-metres');bvh=BVH(test)
    crossed=[];tested=0
    # Same-source adjacent faces are expected to touch; only transverse
    # nonadjacent repaired source triangles can establish self penetration.
    for batch in bvh.candidate_batches(bvh,0,batch_size=2048):
        batch=batch[batch[:,0]<batch[:,1]]
        if not len(batch):continue
        shared=(f[batch[:,0],:,None]==f[batch[:,1],None,:]).any(axis=(1,2));batch=batch[~shared]
        if not len(batch):continue
        result=triangle_pairs(tri[batch[:,0]],tri[batch[:,1]]);tested+=len(batch)
        crossed.extend(batch[result['relation']=='crossing'].tolist())
    source_tri=source[f]
    source_cross=np.cross(source_tri[:,1]-source_tri[:,0],source_tri[:,2]-source_tri[:,0])
    expected=np.einsum('nij,nj->ni',rotations[f].mean(1),source_cross)
    cosine=np.einsum('ni,ni->n',cross,expected)/(np.linalg.norm(cross,axis=1)*np.linalg.norm(expected,axis=1))
    inversions=int((cosine<=0).sum())
    report={'schema':'explicit-source-perforator-section-transport/v1','transports':transports,
        'sourceCatalogueGeometryUnchanged':True,'otherVesselVerticesAndNormalsExactlyUnchanged':True,
        'verticesChangedBeyondCommonField':int(moved.sum()),'sourceFacesChangedBeyondCommonField':int(changed.sum()),
        'newZeroAreaFaces':int((area<=1e-18).sum()),'repairedSelfTrianglePairsTested':tested,
        'invertedAgainstTransportedSourceFrame':inversions,'minimumOrientationCosine':float(cosine.min()),
        'repairedSelfTransverseCrossings':len(crossed),'selfCrossingPairs':crossed,
        'anatomicalAcceptance':False,'minimumClearanceAndSharedMouthChecksPending':True,'elapsedSeconds':time.perf_counter()-started}
    (stage/'section-transport-report.json').write_text(json.dumps(report,indent=2)+'\n')
    status='candidate-numerical-gates-passed-anatomy-unverified' if not crossed and not inversions else 'rejected-repaired-self-crossings-or-orientation'
    previous_report=json.loads((previous/'candidate-report.json').read_text())
    for key in ['analyticProof','minimumVertexJacobianDeterminant','minimumFaceCentroidJacobianDeterminant','minimumFaceOrientationCosineAgainstCentroidJacobian']:
        previous_report.pop(key,None)
    previous_report.update(status=status,commonMapReport=str(previous/'candidate-report.json'),
        explicitSourceSectionOverride=report,contractSha256=hashlib.sha256((stage/'contract.json').read_bytes()).hexdigest(),
        analyticProofScope='Common transported field only; explicit rigid source contours additionally preserve pairwise distances, while side sheets and interfaces require separate audit.')
    (out/'candidate-report.json').write_text(json.dumps(previous_report,indent=2)+'\n')
    (out/'checkpoint.json').write_text(json.dumps({'status':status,'elapsedSeconds':report['elapsedSeconds']})+'\n')
    print(json.dumps({k:report[k]for k in ['verticesChangedBeyondCommonField','sourceFacesChangedBeyondCommonField','repairedSelfTransverseCrossings','elapsedSeconds']}))
    print(json.dumps([dict(source=t['sourceId'],component=t['component'],areaRatio=t['minimumAreaRatio'],sweepJacobian=t['estimatedMinimumSweepJacobian'])for t in transports]))


if __name__=='__main__':main()
