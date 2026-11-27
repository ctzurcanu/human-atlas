"""One local FM aperture probe anchored to actual midsagittal cortical lips."""
from pathlib import Path
import sys,json,time,importlib.util
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.distance import barycentric
BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'

def main():
    began=time.monotonic();out=BASE/'foramen-magnum-source-plane-independent-v4';out.mkdir(exist_ok=False)
    inv=json.loads((BASE/'source-inventory.json').read_bytes());vs=[];fs=[];owners=[]
    for pid,row in inv['sourceBindings'].items():
        if not pid.startswith('ZA:Occipital bone:'):continue
        assert digest_file(row['mesh'])==row['sha256'];z=np.load(row['mesh']);fs.append(z['faces']+sum(len(q)for q in vs));vs.append(z['vertices']);owners += [dict(partId=pid,sourceFace=i)for i in range(len(z['faces']))]
    v,ix=np.unique(np.vstack(vs),axis=0,return_inverse=True);f=ix[np.vstack(fs)]
    spec=importlib.util.spec_from_file_location('ap',Path(__file__).with_name('identify-native-vertebral-apertures-independent-v1.py'));mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    contours,diagnostic=mod.literal_section(v[:,[1,0,2]],f,0.);assert len(contours)==2 and not diagnostic['noncycleSectionComponents']
    front=max(contours,key=lambda q:q['center_m'][2]);back=min(contours,key=lambda q:q['center_m'][2]);b=np.array(front['worldPoints_m'])[:,[1,0,2]];o=np.array(back['worldPoints_m'])[:,[1,0,2]]
    bi=int(b[:,2].argmin());oi=int(o[:,2].argmax());basion=b[bi];opisthion=o[oi]
    def witness(point,row,i):
        ids=np.array([row['sourceFaceIds'][i-1],row['sourceFaceIds'][i]]);ww,valid=barycentric(np.broadcast_to(point,(len(ids),3)),v[f[ids]]);k=int(np.argmax(np.min(ww,axis=1)));return dict(**owners[int(ids[k])],analysisFace=int(ids[k]),point_m=point,barycentric=ww[k],classification='Actual cortical lip candidate; anatomical identification requires rendered review.')
    n=np.cross(np.array([1.,0,0]),opisthion-basion);n/=np.linalg.norm(n)
    if n[1]<0:n=-n
    x=np.array([1.,0,0]);z=np.cross(x,n);basis=np.c_[x,n,z];center=(basion+opisthion)/2;offset=.00025;probe=center+n*offset
    cc,diag=mod.literal_section(v@basis,f,float(probe@n));result=[]
    for i,c in enumerate(cc):
        if not c['holeByNestedContourParity']or not c['centroidInside']:continue
        centerWorld=np.array(c['center_m'])@basis.T
        if abs(centerWorld[0])>.008:continue
        xyz=np.array(c['worldPoints_m'])@basis.T;result.append(dict(contourId=i,worldPoints_m=xyz,center_m=centerWorld,area_m2=c['area_m2'],sourcePlaneCentroidClearance_m=c['centroidClearance_m'],sourceFaceProvenance=[owners[k]for k in c['sourceFaceIds']],
            localClosedSectionCurve=True,sourcePlaneNormal=n,containingSourceContours=c['containingContourIds']))
    atomic_json(out/'evidence.json',dict(schema='native-foramen-magnum-local-source-plane/v4',sourceCatalogueSha256=inv['sourceCatalogueSha256'],sourceInventorySha256=digest_file(BASE/'source-inventory.json'),
        midsagittalSourceLipCandidates=dict(anteriorBasion=witness(basion,front,bi),posteriorOpisthion=witness(opisthion,back,oi)),lipChordLength_m=float(np.linalg.norm(opisthion-basion)),
        apertureProbePlane=dict(point_m=probe,normal=n,cranialOffsetFromLipChord_m=offset,notADuralEntryBoundary=True),selectedLocalApertures=result,
        allProbeContours=[dict(center_m=np.array(c['center_m'])@basis.T,area_m2=c['area_m2'],hole=c['holeByNestedContourParity'],centroidInside=c['centroidInside'])for c in cc],sectionDiagnostics=diag,
        wholeOccipitalSourceDefectsRetained=dict(boundaryEdges=40,nonmanifoldEdges=39),noGeometryOrPartitionChanged=True,closedSolidNotAsserted=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began))
    print(json.dumps(dict(evidence=str(out/'evidence.json'),lipCandidates_mm=[(basion*1000).tolist(),(opisthion*1000).tolist()],normal=n.tolist(),apertures=[dict(center_mm=(r['center_m']*1000).tolist(),area_mm2=r['area_m2']*1e6,clearance_mm=r['sourcePlaneCentroidClearance_m']*1000,points=len(r['worldPoints_m']))for r in result],allContours=len(cc),badComponents=len(diag['noncycleSectionComponents']),elapsedSeconds=time.monotonic()-began)))

if __name__=='__main__':main()
