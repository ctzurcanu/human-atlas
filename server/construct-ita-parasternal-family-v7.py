"""One coupled C2 parasternal native-source arterial transport diagnostic.

No artery calibre, other anatomy, face ownership, UV, root or distant outlet
is reduced/changed. Existing unmatched outer casts remain unmatched casts.
"""
from pathlib import Path
import json, sys, time
import numpy as np
from scipy.interpolate import BSpline, CubicSpline
from scipy.optimize import minimize, LinearConstraint
from scipy.spatial.transform import Rotation
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, atomic_npz, digest_file, identity
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings, align_station_correspondence
from geometry_tools.algorithms.tube_sections import redirect_local_sections, section_metrics, bend_diagnostics
from geometry_tools.audit import topology
from geometry_tools.contact import self_contact_summary

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/thoracic-parietal-vascular-family-v1'
OUT=BASE/'one-parasternal-family-candidate-v8'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def smooth(x):
    t=np.clip(x,0,1);return 6*t**5-15*t**4+10*t**3

def clip_xy(triangle,lo,hi):
    p=list(triangle)
    for dim,bound,sign in [(0,lo[0],1),(0,hi[0],-1),(1,lo[1],1),(1,hi[1],-1)]:
        out=[]
        if not len(p):return np.empty((0,3))
        for a,b in zip(p,p[1:]+p[:1]):
            da=sign*(a[dim]-bound);db=sign*(b[dim]-bound)
            if da>=0:out.append(a)
            if (da>=0)!=(db>=0):out.append(a+(b-a)*(da/(da-db)))
        p=out
    return np.array(p)

def station_data(m,reference=None):
    v,first,inv=np.unique(m.vertices,axis=0,return_index=True,return_inverse=True);f=inv[m.faces]
    courses=recover_station_rings(v,f); result=[]
    for course in courses:
        rr=list(course['rings'])
        if reference is None:
            if v[rr[-1]].mean(0)[1]>v[rr[0]].mean(0)[1]:rr=rr[::-1]
        elif np.linalg.norm(v[rr[-1]].mean(0)-reference)<np.linalg.norm(v[rr[0]].mean(0)-reference):rr=rr[::-1]
        rr,phase=align_station_correspondence(v,f[course['face_ids']],rr); rr=np.array(rr)
        result.append(dict(rings=rr,centers=v[rr].mean(1)))
    return v,inv,f,result

def fit_field(s,window,desired,lower,active):
    # One convex constrained least-squares fit; degree-five B-spline is C4
    # inside and the explicit quintic support window is C2 at its boundary.
    t=np.r_[[s[0]]*6,np.arange(25,s[-1],25),[s[-1]]*6]
    spline=BSpline(t,np.eye(len(t)-6),5); basis=spline(s)*window[:,None]
    dense=np.linspace(s[0],s[-1],500); wd=smooth((dense-s[12])/(s[26]-s[12]))
    values=spline(dense)*wd[:,None]; h=dense[1]-dense[0]; second=np.diff(values,2,axis=0)/h**2
    H=basis.T@basis+2500*second.T@second+1e-8*np.eye(basis.shape[1]); g=basis.T@desired
    fun=lambda c:.5*c@H@c-g@c
    solve=minimize(fun,np.full(basis.shape[1],max(float(lower[active].max()),0.)),jac=lambda c:H@c-g,
        constraints=[LinearConstraint(basis[active],lower[active],np.inf)],method='SLSQP',options=dict(maxiter=100,ftol=1e-10))
    if not solve.success:raise RuntimeError('The one constrained C2 fit failed: '+solve.message)
    field=basis@solve.x
    if (field[active]<lower[active]-1e-7).any():raise RuntimeError('Declared source landmark constraints not met')
    return field,dict(success=bool(solve.success),iterations=solve.nit,message=solve.message,coefficients_mm=solve.x,knots_mm=t,degree=5,
        supportWindow='quintic smoothstep station12 arclength to station26 arclength',landmarkConstraintResidual_mm=field[active]-lower[active])

def apply_course(source,v,inv,rings,target,rotations,changed):
    profiles=v[rings]; moved=np.einsum('sij,svj->svi',rotations,profiles-profiles.mean(1)[:,None])+target[:,None]
    encoded=moved.astype('<f4').astype(float); encoded[~changed]=profiles[~changed]
    output=v.copy();output[rings]=encoded; raw=output[inv]
    station=np.full(len(v),-1,int)
    for j,r in enumerate(rings):station[r]=j
    normal=np.array(source.attributes['normals'],copy=True)
    for j in np.flatnonzero(changed):
        used=np.flatnonzero(station[inv]==j); n=source.attributes['normals'][used].astype(float)/32767
        normal[used]=np.clip(np.rint((n@rotations[j].T)*32767),-32767,32767).astype('<i2')
    attrs={k:z.copy()for k,z in source.attributes.items()};attrs['normals']=normal
    before=section_metrics(profiles[None]);after=section_metrics(encoded[None]);epsilon=float(np.linalg.norm(encoded-moved,axis=2).max())
    perimeter=before['perimeter'][0]; bound=2*epsilon*perimeter+12*epsilon**2+256*np.finfo(float).eps*before['area'][0]
    assert np.all(abs(after['area'][0]-before['area'][0])<=bound)
    return raw,attrs,dict(maximumRelativeAreaError=float((abs(after['area']-before['area'])/before['area']).max()),
        maximumDiameterChange_m=float(abs(after['diameter']-before['diameter']).max()),maximumPointRoundoff_m=epsilon,
        maximumProperOrthogonalityError=float(abs(rotations.transpose(0,2,1)@rotations-np.eye(3)).max()),minimumProperDeterminant=float(np.linalg.det(rotations).min()),
        changedStations=np.flatnonzero(changed),sampledBends=bend_diagnostics(encoded.mean(1),encoded[None])),encoded,bound

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False)
    evidence=json.loads((BASE/'ita-actual-interface-landmarks-v3/evidence.json').read_bytes())
    contract=json.loads((BASE/'ita-distal-closure-context-v6/evidence.json').read_bytes());assert contract['sourceCatalogueSha256']==EXPECTED
    rays=json.loads((BASE/'ita-actual-interface-landmarks-v3/independent-parasternal-rays.json').read_bytes());assert digest_file(BASE/'ita-actual-interface-landmarks-v3/independent-parasternal-rays.json')==contract['sourceRaysSha256']
    bindings=evidence['sourceBindings'];lo,hi=np.array(contract['conservativeCorridorBounds_m']);results=[];candidateBindings={}
    for side in ['l','r']:
        pid='ZA:Internal thoracic artery.'+side; source=Mesh.load(bindings[pid]['path'])
        d=np.load(BASE/'ita-actual-interface-landmarks-v3'/(side+'-ita-source-stations.npz'));v=d['physical_vertices'];inv=d['source_to_physical'];rings=d['rings'];centers=d['centers'];profiles=v[rings]
        s=np.r_[0,np.cumsum(np.linalg.norm(np.diff(centers,axis=0),axis=1))]*1000;window=smooth((s-s[12])/(s[26]-s[12]));sign=1 if side=='l'else -1
        desired=np.array([max(q.get('parasternalX_m',float(sign*centers[j,0]))-sign*centers[j,0],0)*1000 for j,q in enumerate(rays['rows'])])
        # Saved left/right source geometry is independently checked below; the
        # sternal source is bilateral, while recipient cartilage is side-specific.
        if side=='r':assert np.max(abs(centers[:,[1,2]]-np.load(BASE/'ita-actual-interface-landmarks-v3/l-ita-source-stations.npz')['centers'][:,[1,2]]))<1e-9
        active=np.array([j>=24 and 'parasternalX_m'in q for j,q in enumerate(rays['rows'])]); lower=desired-1.
        dx,xfit=fit_field(s,window,desired,lower,active)
        target=centers.copy();target[:,0]+=sign*dx/1000
        cartilage=[('ZA:Costal cartilage of '+n+' rib.'+side,Mesh.load(bindings['ZA:Costal cartilage of '+n+' rib.'+side]['path'])) for n in ['first','second','third','fourth','fifth','sixth','seventh']]
        required=np.zeros(len(rings));witness=[]
        for j,c in enumerate(target):
            radius=np.linalg.norm(profiles[j]-centers[j],axis=1).max(); band=radius+.0005; best=None
            for rid,m in cartilage:
                tt=m.triangles; ids=np.flatnonzero(np.all(tt[:,:,:2].max(1)>=c[:2]-band,axis=1)&np.all(tt[:,:,:2].min(1)<=c[:2]+band,axis=1))
                for fid in ids:
                    polygon=clip_xy(tt[fid],c[:2]-band,c[:2]+band)
                    if len(polygon):
                        point=polygon[np.argmin(polygon[:,2])]
                        if best is None or point[2]<best['point_m'][2]:best=dict(partId=rid,sourceFace=int(fid),point_m=point)
            if best:
                required[j]=max((c[2]+radius+.0005-best['point_m'][2])*1000,0.)
                witness.append(dict(station=j,projectionHalfWidth_m=band,posteriorEnvelope=best,requiredPosteriorOffset_mm=required[j]))
        zactive=(np.arange(len(rings))>=24)&(required>0); zdesired=np.maximum(required,0)
        dz,zfit=fit_field(s,window,zdesired,required,zactive)
        target[:,2]-=dz/1000; target[:13]=centers[:13]
        moved,rot=redirect_local_sections(profiles[None],target,dtype=float);changed=np.any(target!=centers,axis=1)
        # The finite-difference rotation at a protected source boundary is
        # exact identity because the field is zero through station12, and
        # station13 also has zero displacement by the C2 window endpoint.
        changed[:13]=False;rot[:13]=np.eye(3)
        raw,attrs,metrics,encoded,bound=apply_course(source,v,inv,rings,target,rot,changed)
        mesh=Mesh(raw,source.faces,source.frame,attrs,dict(source.source,oneC2ParasternalCandidate=True,sourceGeometrySha256=source.geometry_sha,anatomicalAcceptance=False))
        path=OUT/(side+'-ita.npz');mesh.save(path)
        if side=='l':
            prior=Mesh.load(BASE/'one-parasternal-family-candidate-v7/l-ita.npz');assert mesh.geometry_sha==prior.geometry_sha and mesh.attributes_sha==prior.attributes_sha
        assert np.array_equal(mesh.faces,source.faces)and np.array_equal(attrs['uvs'],source.attributes['uvs'])
        protected=np.isin(inv,np.ravel(rings[:13]));assert np.array_equal(mesh.vertices[protected],source.vertices[protected]);assert np.array_equal(attrs['normals'][protected],source.attributes['normals'][protected])
        selfgate,_=self_contact_summary(mesh,max_pairs=250000,max_tests=2000000);assert selfgate['complete_within_pair_budget']
        fieldpath=OUT/(side+'-ita-field.npz');atomic_npz(fieldpath,physical_vertices=v,source_to_physical=inv,rings=rings,source_centers=centers,target_centers=target,rotation=rot,
            encoded_sections=encoded,changed_stations=changed,area_roundoff_bound=bound,source_faces=source.faces,dx_mm=dx,dz_mm=dz)
        row=dict(partId=pid,candidate=str(path),candidateSha256=digest_file(path),field=str(fieldpath),fieldSha256=digest_file(fieldpath),source=bindings[pid],metrics=metrics,
            xFit=xfit,zFit=zfit,cartilageEnvelopeWitnesses=witness,wholeTopology=topology(mesh),wholeSelfContact=selfgate,sourceOuterDiameterRange_m=[float(d['diameters'].min()),float(d['diameters'].max())],originFirst13StationsLiteral=True)
        candidateBindings[pid]=row;results.append(row)
        # Existing branch-mouth cast centroids agree, but the complete rims
        # do not. Apply ONE common rigid pose to all source mouths, then a C2
        # proper-profile transition to literal source at45 mm physical arc.
        pivot=centers[-1];D=target[-1]-pivot;R=rot[-1]
        for name in ['Musculophrenic artery','Superior epigastric artery']:
            rid='ZA:'+name+'.'+side;src=Mesh.load(bindings[rid]['path']);vv,ii,ff,courses=station_data(src,pivot)
            selected=min(courses,key=lambda z:np.linalg.norm(z['centers'][0]-pivot)); rr=selected['rings'];cc=selected['centers'];arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(cc,axis=0),axis=1))]
            weight=1-smooth((arc-.003)/.042);common=(cc-pivot)@R.T+pivot+D;tt=cc+weight[:,None]*(common-cc);rrr=redirect_local_sections(vv[rr][None],tt,dtype=float)[1]
            # Blend endpoint roll by an actual proper SO(3) rotation. The
            # literal source cast relationship remains the same within3 mm.
            relative=Rotation.from_matrix(R@rrr.transpose(0,2,1)).as_rotvec();blend=Rotation.from_rotvec(weight[:,None]*relative).as_matrix();rrr=blend@rrr
            same=weight==0;rrr[same]=np.eye(3);tt[same]=cc[same];ch=~same
            rraw,aat,mm,ee,bb=apply_course(src,vv,ii,rr,tt,rrr,ch)
            cm=Mesh(rraw,src.faces,src.frame,aat,dict(src.source,oneC2ParasternalCandidate=True,sourceGeometrySha256=src.geometry_sha,anatomicalAcceptance=False))
            cpath=OUT/(side+'-'+('musculophrenic'if name.startswith('Musculo')else'superior-epigastric')+'.npz');cm.save(cpath)
            if side=='l'and name.startswith('Musculo'):
                prior=Mesh.load(BASE/'one-parasternal-family-candidate-v7/l-musculophrenic.npz');assert cm.geometry_sha==prior.geometry_sha and cm.attributes_sha==prior.attributes_sha
            fpath=cpath.with_name(cpath.stem+'-field.npz');atomic_npz(fpath,physical_vertices=vv,source_to_physical=ii,rings=rr,source_centers=cc,target_centers=tt,rotation=rrr,
                encoded_sections=ee,changed_stations=ch,area_roundoff_bound=bb,source_faces=src.faces,common_pivot=pivot,common_translation=D,common_rotation=R,weight=weight)
            untouched=np.ones(len(vv),bool);untouched[rr[ch]]=False;keep=untouched[ii];assert np.array_equal(cm.vertices[keep],src.vertices[keep]);assert np.array_equal(aat['normals'][keep],src.attributes['normals'][keep]);assert np.array_equal(cm.faces,src.faces)
            sg,_=self_contact_summary(cm,max_pairs=250000,max_tests=2000000);assert sg['complete_within_pair_budget']
            crow=dict(partId=rid,candidate=str(cpath),candidateSha256=digest_file(cpath),field=str(fpath),fieldSha256=digest_file(fpath),source=bindings[rid],metrics=mm,wholeTopology=topology(cm),wholeSelfContact=sg,
                unchangedControlsBeyond45mmAndDisconnectedComponentsLiteral=True,completeRootCommonRigidPose=True,unmatchedSourceRimStillNotAnAnnulus=True)
            candidateBindings[rid]=crow;results.append(crow)
        atomic_json(OUT/(side+'-stage-complete.json'),dict(side=side,elapsedSeconds=time.monotonic()-began))
    # Every potentially changed triangle must lie inside the complete-context
    # corridor. Exact untouched controls outside it remain source-owned.
    changedBounds=[]
    for row in results:
        src=Mesh.load(row['source']['path']);dst=Mesh.load(row['candidate']);changed=np.any(src.vertices!=dst.vertices,axis=1);faces=np.flatnonzero(changed[src.faces].any(1));pts=np.vstack([src.vertices[src.faces[faces]].reshape(-1,3),dst.vertices[dst.faces[faces]].reshape(-1,3)])
        bounds=np.array([pts.min(0),pts.max(0)]);within=bool(np.all(bounds[0]>=lo)&np.all(bounds[1]<=hi));changedBounds.append(dict(partId=row['partId'],sourceFaceIds=faces,sourceAndCandidateChangedSupportBounds_m=bounds,insideCompleteContextCorridor=within));assert within
    result=dict(schema='one-coupled-native-ita-parasternal-source-transport/v8',sourceCatalogueSha256=EXPECTED,
        wholeFamilyContractSha256=digest_file(BASE/'whole-family-contract-v4/contract.json'),distalClosureContractSha256=digest_file(BASE/'ita-distal-closure-context-v6/evidence.json'),
        sourceRaysSha256=digest_file(BASE/'ita-actual-interface-landmarks-v3/independent-parasternal-rays.json'),movingOwners=results,changedSupport=changedBounds,
        allOtherAnatomyFixed=True,allFacesUVAndCalibreRetainedWithinEncodingRounding=True,sourceSubclavianAndTerminalUnmatchedOuterInterfacesExplicitlyRetained=True,
        finiteWallOrBloodLumenEstablished=False,all656NeighborAuditPending=True,acceptedForInstallation=False,
        previousIncompleteStageSha256=digest_file(BASE/'one-parasternal-family-candidate-v7/failed-stage.json'),
        identicalPreviouslySavedGeometryAndNormalsReplayed=True,noSecondCourseTrial=True,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)))
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),owners=len(results),elapsedSeconds=result['elapsedSeconds'],maximumAreaError=max(q['metrics']['maximumRelativeAreaError']for q in results),selfCrossings=[(q['partId'],q['wholeSelfContact']['crossing_face_pairs'])for q in results])))

if __name__=='__main__':main()
