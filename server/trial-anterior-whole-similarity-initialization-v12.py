"""One reversible whole-source initialization, using complete junction unions."""
from pathlib import Path
import sys,time,json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.fit import landmark_transform
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.audit import topology,winding_contains

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/arterial-source-opportunities-v1';OUT=BASE/'anterior-whole-similarity-initialization-v12';CORR=BASE/'anterior-donor-counterpart-curves-v10/proof.json';RECIP=BASE/'anterior-source-interface-feasibility-v9/proof.json'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';BP=ROOT/'public/models/atlas.json';EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de';FRAME='model:local-ta98:assembled-metres'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def moments(curves):
    aa=np.vstack(curves);bb=np.vstack([np.roll(x,-1,axis=0)for x in curves]);length=np.linalg.norm(bb-aa,axis=1);total=length.sum();origin=aa.mean(0);a,b=aa-origin,bb-origin
    mu=(length[:,None]*(a+b)/2).sum(0)/total
    second=sum(w*(2*np.outer(x,x)+np.outer(x,y)+np.outer(y,x)+2*np.outer(y,y))/6 for w,x,y in zip(length,a,b))/total;cov=second-np.outer(mu,mu);values,basis=np.linalg.eigh(cov);normal=basis[:,0];uv=basis[:,1:];weighted=[]
    for c in curves:
        xy=(c-origin)@uv;cross=xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1);signed=cross.sum()/2
        if signed!=0:
            cen=((xy+np.roll(xy,-1,axis=0))*cross[:,None]).sum(0)/(6*signed);weighted.append((abs(signed),origin+uv@cen))
    area=sum(x[0]for x in weighted);acenter=sum(w*q for w,q in weighted)/area if area else origin+mu
    return dict(lineIntegralCentroid_m=origin+mu,lineIntegralCovariance_m2=cov,lineIntegralPerimeter_m=float(total),geometricRmsSpread_m=float(np.sqrt(max(0,np.trace(cov)))),maximumBoundaryRadiusFromLineCentroid_m=float(np.linalg.norm(aa-origin-mu,axis=1).max()),principalRmsSpreads_m=np.sqrt(np.maximum(values,0)),diagnosticPlaneNormal=normal,componentCount=len(curves),sumOfAbsoluteProjectedComponentAreas_m2=float(area),projectedAreaIntegralCentroid_m=acenter,lineVsProjectedAreaCentroidSeparation_m=float(np.linalg.norm(acenter-origin-mu)),areaMomentsAreDiagnosticNotAnAnnulus=True,geometricSpreadIsNotStatisticalConfidenceOrAnatomicalLocalizationError=True)

def distance_to_curves(points,curves):
    a=np.vstack(curves);b=np.vstack([np.roll(x,-1,axis=0)for x in curves]);e=b-a;den=(e*e).sum(1);delta=points[:,None,:]-a[None,:,:];t=np.clip(np.einsum('nki,ki->nk',delta,e)/den,0,1);q=a[None,:,:]+t[:,:,None]*e[None,:,:];return np.linalg.norm(points[:,None,:]-q,axis=2).min(1)

def curve_residual(first,second,step=.000025):
    def direction(src,dst):
        samples=[];largest=0.;quad=[];weights=[];g=np.array([.5-np.sqrt(3/5)/2,.5,.5+np.sqrt(3/5)/2]);gw=np.array([5/18,4/9,5/18])
        for c in src:
            for a,b in zip(c,np.roll(c,-1,axis=0)):
                length=float(np.linalg.norm(b-a));n=max(1,int(np.ceil(length/step)));largest=max(largest,length/n);samples.extend(a+np.arange(n)[:,None]/n*(b-a));quad.extend(a+g[:,None]*(b-a));weights.extend(length*gw)
        sd=distance_to_curves(np.array(samples),dst);qd=distance_to_curves(np.array(quad),dst);w=np.array(weights);roundoff=1e-12
        return dict(sampledMaximumLower_m=float(sd.max()),continuousMaximumUpper_m=float(sd.max()+largest/2+roundoff),fullCurveCoverageMaximumSampleInterval_m=largest,samples=len(samples),arcLengthGauss3Rms_m=float(np.sqrt((w*qd*qd).sum()/w.sum())),continuousBoundUsesDistanceOneLipschitz=True)
    return dict(sourceUnionToActualRecipientRim=direction(first,second),actualRecipientRimToSourceUnion=direction(second,first))

def main():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED;resume='--resume-identical-matrix' in sys.argv
    if resume:
        assert OUT.exists()and not (OUT/'proof.json').exists();saved_matrix=np.array(Mesh.load(OUT/'FJ1654.npz').source['proper_whole_similarity_matrix'])
    else:OUT.mkdir(exist_ok=False)
    candidate_dir=OUT/'alias-preserving-replay'if resume else OUT
    if resume:candidate_dir.mkdir(exist_ok=False)
    a=json.loads(CAT.read_text());b=json.loads(BP.read_text());ap={p['id']:p for p in a['parts']};bp={p['id']:p for p in b['parts']};corr=json.loads(CORR.read_text());recip=json.loads(RECIP.read_text());reader=AtlasReader();interfaces=[];xx=[];yy=[]
    for source,target in zip(corr['rows'],recip['interfaces']):
        assert source['side']==target['side']
        for role,components,rim in [('ICA-root',source['sourceRootIntersection']['components'],target['nativeIcaRootRim']),('distal-ACA-continuation',source['sourceDistalContinuationIntersection']['components'],target['nativeAcaDistalRim'])]:
            curves=[np.array(c['sourceCurvePoints_m'])for c in components];native=[np.array(rim['worldPoints_m'])];sm,tm=moments(curves),moments(native);xx.append(sm['lineIntegralCentroid_m']);yy.append(tm['lineIntegralCentroid_m']);interfaces.append(dict(side=source['side'],role=role,sourceCompleteCurves=curves,targetActualRim=native,sourceMoments=sm,targetMoments=tm,sourcePartId=source['donorAca'],recipientPartId=rim['partId'],recipientSourceBoundaryLoop=rim['actualSourceBoundaryLoop'],actualSourceFacePairs=[pair for c in components for pair in c['actualSourceFacePairs']],coarseRegionLandmarkOnly=True,exactAnnulusCorrespondenceClaimed=False))
    spec=dict(mode='similarity',moving=['FJ1654','FJ1654M','FJ1655','FJ1682','FJ1682M'],sourceFrame=corr['sourceFrame'],targetFrame=FRAME,source_landmarks=np.array(xx),target_landmarks=np.array(yy),weights=np.ones(4),scale_bounds=[.5,2.],oneTrialOnly=True,landmarks='Line-integral centroids of complete named source junction boundary unions, including all three distal loops; actual native recipient boundary line-integral centroids',geometricScaleBoundsAreMathematicalInitializationBoundsNotAnatomicalApproval=True,noNativeAnatomyMoves=True,noDeformationOrWallConstruction=True)
    atomic_json(OUT/'initialization-spec.json',json.loads(json.dumps(spec,default=clean)))
    if resume:
        matrix=saved_matrix;rr=np.array(xx)@matrix[:3,:3].T+matrix[:3,3]-yy
        fit=dict(scale=float(np.linalg.norm(matrix[:3,0])),rotation_determinant=float(np.linalg.det(matrix[:3,:3]/np.linalg.norm(matrix[:3,0]))),landmark_rms_m=float(np.sqrt(np.mean(np.sum(rr*rr,axis=1)))),maximum_landmark_error_m=float(np.linalg.norm(rr,axis=1).max()),resumedOriginalMatrixWithoutAnotherFit=True)
    else:matrix,fit=landmark_transform(spec['source_landmarks'],spec['target_landmarks'],weights=spec['weights'],similarity=True,scale_bounds=spec['scale_bounds'])
    linear=matrix[:3,:3];scale=float(np.linalg.norm(linear[:,0]));rotation=linear/scale;translation=matrix[:3,3];assert np.linalg.det(rotation)>0
    # Independent centered affine replay and similarity invariants; no call
    # back into the transform implementation creates the candidate arrays.
    center=np.mean(xx,axis=0);center_after=linear@center+translation;gram=linear.T@linear;assert np.max(abs(gram-scale*scale*np.eye(3)))<1e-12
    residuals=[]
    for row in interfaces:
        source=[c@linear.T+translation for c in row['sourceCompleteCurves']];res=curve_residual(source,row['targetActualRim']);sm=row['sourceMoments'];tm=row['targetMoments'];normal=rotation@sm['diagnosticPlaneNormal'];angle=float(np.degrees(np.arccos(np.clip(abs(normal@tm['diagnosticPlaneNormal']),-1,1))));center_error=float(np.linalg.norm(linear@sm['lineIntegralCentroid_m']+translation-tm['lineIntegralCentroid_m']));residuals.append(dict(side=row['side'],role=row['role'],coarseLandmarkResidual_m=center_error,sourceJunctionGeometricRmsSpreadAfterScale_m=scale*sm['geometricRmsSpread_m'],sourceJunctionMaximumRadiusAfterScale_m=scale*sm['maximumBoundaryRadiusFromLineCentroid_m'],diagnosticPlaneAngle_degrees=angle,fullCurveRimResidual=res,coarseLandmarkDoesNotEstablishAnnularInterface=True))
    candidates={};snapshots=[];replay=[]
    for pid in spec['moving']:
        old=reader.read_part(b,bp[pid],frame=corr['sourceFrame'],catalogue_path=BP);unique,alias=np.unique(old.vertices,axis=0,return_inverse=True);positions=(unique@linear.T+translation)[alias];independent=(old.vertices-center)@linear.T+center_after;error=float(np.linalg.norm(positions-independent,axis=1).max());assert error<1e-12;attrs={k:q.copy()for k,q in old.attributes.items()};attrs['normals']=old.attributes['normals'].astype(float)@rotation.T
        candidate=Mesh(positions,old.faces,FRAME,attrs,dict(old.source,proper_whole_similarity_matrix=matrix.tolist(),diagnostic_initialization_only=True,noAnatomicalAcceptance=True));path=candidate_dir/(pid+'.npz');candidate.save(path);saved=Mesh.load(path);assert np.array_equal(saved.vertices,positions)and np.array_equal(saved.faces,old.faces)and np.array_equal(saved.attributes['normals'],attrs['normals']);assert topology(saved)['closed_oriented_manifold'];candidates[pid]=candidate;snapshots.append(dict(id=pid,path=str(path),sha256=digest_file(path),sourceGeometrySha256=old.geometry_sha,candidateGeometrySha256=candidate.geometry_sha,sourceBankSha256=old.source['bank_sha256']));replay.append(dict(id=pid,maximumIndependentTransformError_m=error,facesLiteral=True,everyLiteralSourceCoordinateEquivalenceClassTransformedOnce=True,exactCoordinateTopologyRemainsClosed=True,sourceUvColorNormalSignsLiteral=all(np.array_equal(saved.attributes[k],old.attributes[k])for k in old.attributes if k!='normals'),normalFieldIsExactlyProperRotationOfSourceWords=True,noIntegerNormalRequantization=True))
    fixed={};bindings={};brainids=[p['id']for p in a['parts']if p['id'].startswith('ZA:')and p['id'].rsplit(':',1)[-1]in ['Frontal lobe','Temporal lobe','Parietal lobe','Occipital lobe']];vv=[];ff=[];off=0
    for pid in ['ZA:Optic chiasm.l','ZA:Optic chiasm.r','ZA:Sphenoid bone:Bone-4']+brainids:
        m=reader.read_part(a,ap[pid],frame=FRAME,catalogue_path=CAT);bindings[pid]=identity(ap[pid])
        if pid in brainids:vv.append(m.vertices);ff.append(m.faces+off);off+=len(m.vertices)
        else:fixed[pid]=m
    fixed['Fixed-ZA-cortex-context']=Mesh(np.vstack(vv),np.vstack(ff),FRAME,{},dict(literalFixedNativeBrainParts=brainids));trees={pid:BVH(m)for pid,m in fixed.items()};topos={pid:topology(m)for pid,m in fixed.items()};packing=[]
    for pid,m in candidates.items():
        tree=BVH(m)
        for target,host in fixed.items():
            data,gate=pair_contacts(m,host,0.,first_bvh=tree,second_bvh=trees[target],max_pairs=250000,max_tests=2000000);assert gate['complete_within_pair_budget'];contain=None
            if not gate['crossing_face_pairs']and not gate['coplanar_overlap_face_pairs']and topos[target]['closed_oriented_manifold']:
                c1=winding_contains(host,m.vertices[:1]);c2=winding_contains(m,host.vertices[:1]);contain=dict(candidateFirstControlInHost=c1,hostFirstControlInCandidate=c2,bothOneConnectedClosedSurface=True)
            packing.append(dict(candidatePartId=pid,fixedPartId=target,actualGlobalBoundaryGate=gate,actualCrossingSourceFacePairs=data['pairs'][data['relation']=='crossing'],hostClosedOrientedVolume=topos[target]['closed_oriented_manifold'],containmentWhenBoundaryDisjoint=contain,clearanceAccepted=False))
    # Strict rim-residual screening is deliberately separate from the coarse
    # centroid objective. The moment spread cannot waive an interface gap.
    limits=dict(curveRimRms_m=.00025,curveRimMaximumUpper_m=.0005,diagnosticPlaneAngle_degrees=20)
    rim_pass=all(r['diagnosticPlaneAngle_degrees']<=limits['diagnosticPlaneAngle_degrees']and all(q['arcLengthGauss3Rms_m']<=limits['curveRimRms_m']and q['continuousMaximumUpper_m']<=limits['curveRimMaximumUpper_m']for q in r['fullCurveRimResidual'].values())for r in residuals)
    crossings=sum(p['actualGlobalBoundaryGate']['crossing_face_pairs']+p['actualGlobalBoundaryGate']['coplanar_overlap_face_pairs']for p in packing)
    after=json.loads(CAT.read_text());now={p['id']:p for p in after['parts']};assert all(identity(now[pid])==sha for pid,sha in bindings.items())and digest_file(CAT)==EXPECTED
    proof=dict(schema='anterior-whole-similarity-coarse-initialization/v12',sourceCatalogueSha256=digest_file(BP),recipientCatalogueSha256=EXPECTED,sourceCurvesProofSha256=digest_file(CORR),recipientRimsProofSha256=digest_file(RECIP),movingClosure=spec['moving'],wholeProperSimilarityMatrix=matrix,fitMetrics=fit,uniformScale=scale,rotationDeterminant=float(np.linalg.det(rotation)),rotationAngle_degrees=float(np.degrees(np.arccos(np.clip((np.trace(rotation)-1)/2,-1,1)))),properSimilarityGramMaximumError=float(abs(gram-scale*scale*np.eye(3)).max()),interfaceMomentsAndAncestry=interfaces,perInterfaceResiduals=residuals,meaningfulRimScreeningLimits=limits,rimResidualScreenPassed=rim_pass,actualFixedPackingChecks=packing,actualCrossingOrCoplanarPairs=crossings,allFixedContextAndUnrelatedArteriesUnchanged=True,independentSavedTransformReplay=replay,candidateSnapshots=snapshots,fitAcceptedForAnatomicalPlacement=False,rejected=bool(not rim_pass or crossings),rejectReason='Coarse whole-source initialization is not an actual matched arterial interface; any failed full-curve residual or fixed-context crossing rejects placement.',wallConstructionOrInstallPermitted=False,sourceClosedCastsRemainClosedCasts=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'proof.json',json.loads(json.dumps(proof,default=clean)))
    print(json.dumps(dict(output=str(OUT/'proof.json'),scale=scale,rotationDegrees=proof['rotationAngle_degrees'],centroidRms_mm=1000*fit['landmark_rms_m'],residuals=[dict(side=r['side'],role=r['role'],centroid_mm=1000*r['coarseLandmarkResidual_m'],sourceToRimMaxUpper_mm=1000*r['fullCurveRimResidual']['sourceUnionToActualRecipientRim']['continuousMaximumUpper_m'],rimToSourceMaxUpper_mm=1000*r['fullCurveRimResidual']['actualRecipientRimToSourceUnion']['continuousMaximumUpper_m'],planeAngle=r['diagnosticPlaneAngle_degrees'])for r in residuals],crossingPairs=crossings,rejected=proof['rejected'],elapsedSeconds=proof['elapsedSeconds'])))
if __name__=='__main__':main()
