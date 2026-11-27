"""One native-profile, skeleton-fixed vertebrobasilar C2 diagnostic.

Landmarks are literal foraminal source guides.  Upper-course constraints are
local source cortical/brainstem tangent planes, not a global solid proof.
No missing branch, finite wall, annulus or blood lumen is invented.
"""
from pathlib import Path
import importlib.util,json,sys,time
import numpy as np
from scipy.interpolate import BSpline,BPoly,CubicSpline
from scipy.spatial.transform import Rotation
from scipy.optimize import minimize
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT,atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_sections import redirect_local_sections,arclength,section_metrics,bend_diagnostics
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
from geometry_tools.algorithms.section_course_qp import project_weighted_section_halfspaces
from geometry_tools.audit import topology
from geometry_tools.contact import self_contact_summary

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/vertebral-native-regional-opportunity-v1'
OUT=BASE/'one-coupled-foraminal-candidate-v9'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'

def clean(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x).__name__)

def smooth(t):
    t=np.clip(t,0,1);return 6*t**5-15*t**4+10*t**3

def station_data(mesh,pivot):
    v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];courses=[]
    for q in recover_station_rings(v,f):
        rings=list(q['rings'])
        if np.linalg.norm(v[rings[-1]].mean(0)-pivot)<np.linalg.norm(v[rings[0]].mean(0)-pivot):rings=rings[::-1]
        rings,_=align_station_correspondence(v,f[q['face_ids']],rings);rings=np.array(rings)
        courses.append(dict(rings=rings,centers=v[rings].mean(1),sourceFaceIds=q['face_ids']))
    return v,inv,courses

def first_segment_hit(tri,p,q):
    d=q-p;e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];h=np.cross(np.broadcast_to(d,e2.shape),e2);det=np.einsum('ij,ij->i',e1,h);ok=abs(det)>1e-20
    f=np.zeros_like(det);f[ok]=1/det[ok];s=p-tri[:,0];u=f*np.einsum('ij,ij->i',s,h);z=np.cross(s,e1);v=f*(z@d);t=f*np.einsum('ij,ij->i',e2,z)
    ids=np.flatnonzero(ok&(u>=0)&(v>=0)&(u+v<=1)&(t>1e-8)&(t<=1))
    return None if not len(ids)else dict(face=int(ids[np.argmin(t[ids])]),fraction=float(t[ids].min()))

def fm_target(c,guide,bone,radius):
    p=np.array(guide['target_m']);axis=np.array(guide['axis']);projected=c-axis*((c-p)@axis);tree=BVH(bone)
    hit=first_segment_hit(bone.triangles,p,projected);end=hit['fraction']if hit else 1.;delta=projected-p
    # Start at the actual central gap and retain the last segment position
    # before the first literal cortical crossing with full profile clearance.
    desired=radius+.00015;lo=0.;hi=end
    if tree.closest_point(p+hi*delta)['distance_m']>=desired:lo=hi
    else:
        for _ in range(32):
            mid=(lo+hi)/2
            if tree.closest_point(p+mid*delta)['distance_m']>=desired:lo=mid
            else:hi=mid
    target=p+lo*delta;nearest=tree.closest_point(target)
    assert nearest['distance_m']>=desired-1e-10 and abs((target-p)@axis)<1e-12
    return target,dict(sourceProjectedPoint_m=projected,actualCentralGap_m=p,firstActualCorticalSegmentHit=hit,retainedSegmentFraction=lo,
        actualSourceClosestWitness=nearest,requiredProfileSphereClearance_m=desired,notDuralEntry=True,notSharedBilateralWaypoint=True)

def lower_course(v,rings,guides,side,fm,bone):
    c=v[rings].mean(1);s=arclength(c);native=CubicSpline(s,c);rows={q['level']:q for q in guides if q['side']==side};keys={8:(c[8],None)};proof=[]
    for level,j in [('C6',25),('C4',50),('C2',70),('C1',90)]:
        g=rows[level];axis=np.array(g['sourceCycleBestFitAxis']);axis*=np.sign(axis@native(s[j],1));target=np.array(g['target_m'])
        keys[j]=(target,axis);proof.append(dict(level=level,station=j,target_m=target,axis=axis,sourceEvidenceSha256=g['sourceEvidenceSha256']))
    # Two neighboring source guide-axis points constrain actual entry course,
    # not only the nearest point of an unrelated bone surface.
    for level,j in [('C2',77),('C1',82)]:
        g=rows[level];axis=np.array(g['sourceCycleBestFitAxis']);axis*=np.sign(axis@native(s[j],1));p=np.array(g['target_m']);target=p+axis*((c[j,1]-p[1])/axis[1]);keys[j]=(target,axis)
        proof.append(dict(level=level+' local source-cycle axis',station=j,target_m=target,axis=axis,initializationNotApertureLip=True))
    target,witness=fm_target(c[136],fm,bone,float(np.linalg.norm(v[rings[136]]-c[136],axis=1).max()));keys[136]=(target,None);keys[155]=(c[155],None)
    ordered=sorted(keys);derivatives=[]
    for j in ordered:
        target,axis=keys[j];d1=np.zeros(3)if axis is None else axis*np.linalg.norm(native(s[j],1))-native(s[j],1)
        d2=np.zeros(3)if axis is None else-native(s[j],2)
        derivatives.append([target-c[j],d1,d2])
    field=BPoly.from_derivatives(s[ordered],derivatives);result=c.copy();take=(np.arange(len(c))>=8)&(np.arange(len(c))<=155);result[take]+=field(s[take]);result[:9]=c[:9];result[155:]=c[155:]
    return result,dict(method='piecewise quintic Hermite displacement, continuous through second derivative',literalSourceSupports=[list(range(9)),list(range(155,171))],guideStations=proof,
        foramenMagnumStation=136,foramenMagnumWitness=witness,knots_m=s[ordered],displacementDerivatives=derivatives)

def upper_coupled(data,clivus,brain_right):
    # Three fields share one literal root translation.  Their local basis
    # vanishes through second derivative at protected support boundaries.
    keys=['l','r','b'];slices={};bases={};count=1
    for key,n,start,end in [('l',205,170,204),('r',205,170,204),('b',37,0,19)]:
        c=data[key]['centers'];s=arclength(c);u=(s-s[start])/(s[end]-s[start]);knots=np.r_[[0.]*6,[.25,.5,.75],[1.]*6];sp=BSpline(knots,np.eye(len(knots)-6),5,extrapolate=False);raw=np.nan_to_num(sp(np.clip(u,0,1)))
        if key=='b':
            w=smooth((s-s[1])/(s[4]-s[1]))*(1-smooth((s-s[15])/(s[19]-s[15])));common=1-smooth((s-s[2])/(s[19]-s[2]));w[np.arange(n)>=19]=0;common[np.arange(n)>=19]=0
        else:
            w=smooth((s-s[170])/(s[175]-s[170]))*(1-smooth((s-s[198])/(s[202]-s[198])));common=smooth((s-s[170])/(s[202]-s[170]));w[:171]=0;common[:171]=0;w[202:]=0;common[202:]=1
        bases[key]=(raw*w[:,None],common,s);slices[key]=slice(count,count+raw.shape[1]);count+=raw.shape[1]
    matrix={}
    for key,(local,common,s)in bases.items():
        B=np.zeros((len(s),count));B[:,0]=common;B[:,slices[key]]=local;matrix[key]=B
    allB=np.vstack(list(matrix.values()));H=allB.T@allB+1e-8*np.eye(count)
    for key,B in matrix.items():
        # Dimensionless index curvature regularization; analytic basis and
        # quintic windows establish C2, not these discrete samples.
        D=np.diff(B,2,axis=0);H+=500*D.T@D
    K=np.linalg.inv(H);weights=[];normals=[];offsets=[];witnesses=[]
    for row in clivus['rows']:
        key='l'if row['owner']=='leftVA'else'b';j=row['station'];c=data[key]['centers'][j];radius=float(np.linalg.norm(data[key]['profiles'][j]-c,axis=1).max())
        for label in ['posteriorClivalCortex','anteriorBrainstemSurface']:
            q=row.get(label)
            if not q:continue
            n=np.array(q['outwardFaceUnitNormal']);p=np.array(q['point_m']);weights.append(matrix[key][j]);normals.append(n);offsets.append(radius+.0001-n@(c-p));witnesses.append(dict(owner=key,station=j,role=label,sourceWitness=q,sphereRadius_m=radius,addedClearance_m=.0001))
            if key=='l':
                # Reflect the query only, then recover its literal right
                # source triangle and normal before applying a constraint.
                p=p*np.array([-1,1,1]);n=n*np.array([-1,1,1]);rid=q['partId'].replace('.l','.r');m=brain_right.get(rid)
                assert m is not None;closest=BVH(m).closest_point(p);p=np.array(closest['point']);tt=m.triangles[closest['face_id']];actual=np.cross(tt[1]-tt[0],tt[2]-tt[0]);actual/=np.linalg.norm(actual);assert actual@n>0.99;n=actual
                cc=data['r']['centers'][j];weights.append(matrix['r'][j]);normals.append(n);offsets.append(radius+.0001-n@(cc-p));witnesses.append(dict(owner='r',station=j,role=label,actualRightSourceClosestWitness=closest,mirrorUsedAsQueryOnly=True,point_m=p,normal=n,sphereRadius_m=radius))
    # Bilateral source junction remains on the actual midsagittal plane.
    for sign in [-1,1]:weights.append(np.eye(count)[0]);normals.append(np.array([sign,0,0]));offsets.append(0.)
    atomic_npz(OUT/'actual-coupled-constraints-before-solve.npz',metric=K,weights=np.array(weights),normals=np.array(normals),offsets_m=np.array(offsets),**{'basis_'+key:B for key,B in matrix.items()})
    atomic_json(OUT/'actual-coupled-witnesses-before-solve.json',json.loads(json.dumps(dict(rows=witnesses,extraRows=['shared junction Dx >= 0','shared junction Dx <= 0'],method='Source local tangent planes; whole-source collision gate still required.',notAProofThatAll3DCoursesAreImpossible=True),default=clean)))
    # Translation-only was rejected with a saved Farkas certificate.  The
    # actual source junction is now one proper rigid pose, not an arbitrary
    # three-ring frozen translation.  Its profile support is measured from
    # each actual rotated ring rather than replaced by an isotropic sphere.
    L=np.linalg.cholesky(K);pivot=data['b']['centers'][0];cache={}
    def evaluate(x):
        if cache.get('x')is not None and np.array_equal(cache['x'],x):return cache['value']
        coef=(L@x[:-3].reshape(count,3))*.001;R=Rotation.from_rotvec(x[-3:]).as_matrix();targets={};rotations={};profiles={}
        for key,B in matrix.items():
            old=data[key]['centers'];common=B[:,0];target=old+B@coef+common[:,None]*((old-pivot)@(R-np.eye(3)).T)
            if key!='b':target[:171]=old[:171]
            else:target[19:]=old[19:]
            rr=redirect_local_sections(data[key]['profiles'][None],target,dtype=float)[1];relative=Rotation.from_matrix(R@rr.transpose(0,2,1)).as_rotvec();rr=Rotation.from_rotvec(common[:,None]*relative).as_matrix()@rr
            if key!='b':rr[:170]=np.eye(3);rr[203:]=R
            else:rr[:2]=R;rr[20:]=np.eye(3)
            targets[key]=target;rotations[key]=rr;profiles[key]=np.einsum('sij,svj->svi',rr,data[key]['profiles']-old[:,None])+target[:,None]
        value=(targets,rotations,profiles,coef,R);cache.update(x=x.copy(),value=value);return value
    planes=[]
    for q in witnesses:
        if 'sourceWitness'in q:n=np.array(q['sourceWitness']['outwardFaceUnitNormal']);p=np.array(q['sourceWitness']['point_m'])
        else:n=np.array(q['normal']);p=np.array(q['point_m'])
        planes.append((q['owner'],q['station'],n,p))
    def constraints(x):
        _,_,pp,_,_=evaluate(x)
        return np.array([((pp[key][j]-p)@n).min()-.0001 for key,j,n,p in planes])*1000
    def objective(x):return .5*x[:-3]@x[:-3]+12.5*(x[-3:]@x[-3:])
    def gradient(x):return np.r_[x[:-3],25*x[-3:]]
    start=time.monotonic();sol=minimize(objective,np.zeros(3*count+3),jac=gradient,constraints=[dict(type='ineq',fun=constraints)],method='SLSQP',options=dict(maxiter=100,ftol=1e-10))
    targets,rotations,_,coefficients,R=evaluate(sol.x);residual=float(max(0,-constraints(sol.x).min()/1000));report=dict(method='single proper-rigid-root and C2 course SLSQP; actual section support against source planes',success=bool(sol.success),iterations=int(sol.nit),message=str(sol.message),maximumResidual_m=residual,elapsedSeconds=time.monotonic()-start,rootProperRotation=R,rootRotationVector_rad=sol.x[-3:],objective=float(sol.fun))
    atomic_json(OUT/'proper-root-solve-receipt.json',json.loads(json.dumps(report,default=clean)))
    atomic_npz(OUT/'proper-root-solve-fields.npz',**{'target_'+k:z for k,z in targets.items()},**{'rotation_'+k:z for k,z in rotations.items()},coefficients_m=coefficients,root_rotation=R)
    if not sol.success or residual>1e-9:raise ValueError('The single proper-root C2 solve failed: '+sol.message+'; residual '+str(residual)+'m')
    return targets,dict(method='one full-3D coupled proper-profile C2 field with one common proper root pose',solver=report,coefficients_m=coefficients,bases=matrix,rootRotation=R,rootPivot_m=pivot,rotations=rotations,
        constraints=witnesses,unmeasuredMidlineBrainstemAndBasilar12to19ExplicitlyUnknown=True,localPlaneConstraintsDoNotCertifyCurvedWholeSurface=True)

def apply(source,v,inv,updates):
    result=v.copy();attrs={k:q.copy()for k,q in source.attributes.items()};fields=[];metrics=[]
    for k,(rings,target,rot,changed)in enumerate(updates):
        profiles=v[rings];centers=profiles.mean(1);moved=np.einsum('sij,svj->svi',rot,profiles-centers[:,None])+target[:,None];encoded=moved.astype('<f4').astype(float);encoded[~changed]=profiles[~changed];result[rings]=encoded
        st=np.full(len(v),-1,int)
        for j,r in enumerate(rings):st[r]=j
        for j in np.flatnonzero(changed):
            ids=np.flatnonzero(st[inv]==j);attrs['normals'][ids]=np.clip(np.rint((source.attributes['normals'][ids].astype(float)@rot[j].T)),-32767,32767).astype('<i2')
        before=section_metrics(profiles[None]);after=section_metrics(encoded[None]);eps=np.linalg.norm(encoded-moved,axis=2).max();bound=2*eps*before['perimeter'][0]+12*eps**2+256*np.finfo(float).eps*before['area'][0];assert np.all(abs(after['area'][0]-before['area'][0])<=bound)
        metrics.append(dict(component=k,relativeAreaError=float((abs(after['area']-before['area'])/before['area']).max()),maximumPointRoundoff_m=float(eps),maximumIsometryError=float(abs(rot.transpose(0,2,1)@rot-np.eye(3)).max()),minimumProperDeterminant=float(np.linalg.det(rot).min()),bends=bend_diagnostics(encoded.mean(1),encoded[None])))
        fields.append(dict(rings=rings,target_centers=target,rotation=rot,encoded_sections=encoded,changed_stations=changed,area_roundoff_bound=bound))
    raw=result[inv];mesh=Mesh(raw,source.faces,source.frame,attrs,dict(source.source,sourceGeometrySha256=source.geometry_sha,oneVertebrobasilarC2Diagnostic=True,anatomicalAcceptance=False));changed=np.any(raw!=source.vertices,axis=1)
    assert np.array_equal(mesh.faces,source.faces)
    for key in attrs:
        if key!='normals':assert np.array_equal(attrs[key],source.attributes[key])
    assert np.array_equal(attrs['normals'][~changed],source.attributes['normals'][~changed])
    return mesh,fields,metrics

def main():
    began=time.monotonic();OUT.mkdir(exist_ok=False);assert digest_file(CAT)==EXPECTED
    inventory=json.loads((BASE/'source-inventory.json').read_bytes());closure=json.loads((BASE/'necessary-basilar-source-closure-v5/contract.json').read_bytes());guides=json.loads((BASE/'source-derived-throat-guides-independent-v5.json').read_bytes());clivus=json.loads((BASE/'posterior-clivus-patch-qualified-rays-independent-v1.json').read_bytes());bindings={**inventory['sourceBindings'],**closure['sourceBindings']};data={}
    for key,name in [('l','ZA:Vertebral artery.l'),('r','ZA:Vertebral artery.r'),('b','ZA:Basilar artery')]:
        z=np.load(BASE/(key+'-literal-source-stations.npz')if key!='b'else BASE/'necessary-basilar-source-closure-v5/basilar-literal-source-stations.npz');d={k:z[k]for k in z.files};d['profiles']=d['physical_vertices'][d['rings']];d['source']=Mesh.load(bindings[name]['mesh']);d['partId']=name;data[key]=d
    reader=AtlasReader();atlas=json.loads(CAT.read_bytes());parts={p['id']:p for p in atlas['parts']};brain={}
    for pid in ['ZA:Medulla oblongata.r','ZA:Pons.r','ZA:Occipital bone:Bone-6']:brain[pid]=reader.read_part(atlas,parts[pid],frame=data['l']['source'].frame,catalogue_path=CAT,catalogue_sha=EXPECTED)
    upper,solve=upper_coupled(data,clivus,brain);atomic_json(OUT/'coupled-solve-before-construction.json',json.loads(json.dumps(solve,default=clean)));results=[];parent={};bone=Mesh.load(bindings['ZA:Occipital bone:Bone-6']['mesh'])
    for key,d in data.items():
        c=d['centers'];target=upper[key];lowerProof=None
        if key!='b':lower,lowerProof=lower_course(d['physical_vertices'],d['rings'],guides['transverseForamina'],key,guides['foramenMagnumGuide'],bone);target[:171]=lower[:171]
        profiles=d['profiles'];rot=redirect_local_sections(profiles[None],target,dtype=float)[1];rot[171:]=solve['rotations'][key][171:]if key!='b'else rot[171:];changed=np.any(target!=c,axis=1);unchanged=~changed
        # Exact source neighborhoods have exact identity orientation, while
        # the shared constant-translation union keeps its original casts.
        if key!='b':rot[:9]=np.eye(3);rot[156:170]=np.eye(3);rot[203:]=solve['rootRotation']
        else:rot=solve['rotations'][key];rot[:2]=solve['rootRotation'];rot[20:]=np.eye(3)
        changed|=np.any(rot!=np.eye(3),axis=(1,2));mesh,fields,metrics=apply(d['source'],d['physical_vertices'],d['source_to_physical'],[(d['rings'],target,rot,changed)]);pid=d['partId'];path=OUT/(key+'-artery.npz');mesh.save(path);arrays=dict(physical_vertices=d['physical_vertices'],source_to_physical=d['source_to_physical'],source_faces=d['source'].faces)
        for fidx,field in enumerate(fields):arrays.update({str(fidx)+'_'+k:v for k,v in field.items()})
        fp=OUT/(key+'-field.npz');atomic_npz(fp,**arrays);selfgate,_=self_contact_summary(mesh,max_pairs=250000,max_tests=2000000);assert selfgate['complete_within_pair_budget'];results.append(dict(partId=pid,source=bindings[pid],candidate=str(path),candidateSha256=digest_file(path),field=str(fp),fieldSha256=digest_file(fp),components=1,metrics=metrics,wholeSelfContact=selfgate,topology=topology(mesh),lowerCourseProof=lowerProof));parent[key]=(c,target,rot)
    # Only native AICA and the two measured lateral-pontine roots fall in the
    # changed basilar support. Their distal components stay literal.
    pivot0=data['b']['centers'];targetB=upper['b'];rotationB=parent['b'][2]
    for pid in ['ZA:'+name+'.'+side for name in ['Anterior inferior cerebellar artery','Lateral pontine branches of basilar artery']for side in ['l','r']]:
        src=Mesh.load(bindings[pid]['mesh']);v,inv,courses=station_data(src,pivot0[8]);updates=[];dependency=[]
        for course in courses:
            rr=course['rings'];cc=course['centers'];j=int(np.argmin(np.linalg.norm(pivot0-cc[0],axis=1)));assert j in [4,12,17];pivot=pivot0[j];R=rotationB[j];D=targetB[j]-pivot;arc=arclength(cc);hold=min(.002,arc[-1]/8);end=min(.035,arc[-1]*.8);weight=1-smooth((arc-hold)/(end-hold));common=(cc-pivot)@R.T+pivot+D;tt=cc+weight[:,None]*(common-cc);rot=redirect_local_sections(v[rr][None],tt,dtype=float)[1];relative=Rotation.from_matrix(R@rot.transpose(0,2,1)).as_rotvec();rot=Rotation.from_rotvec(weight[:,None]*relative).as_matrix()@rot;same=weight==0;rot[same]=np.eye(3);tt[same]=cc[same];updates.append((rr,tt,rot,~same));dependency.append(dict(componentSourceFaceIds=course['sourceFaceIds'],parentStation=j,rootCommonProperPose=True,literalDistalStartsAtArc_m=end,originalCastNotAnnulus=True))
        mesh,fields,metrics=apply(src,v,inv,updates);path=OUT/(identity(pid)[:16]+'.npz');mesh.save(path);arrays=dict(physical_vertices=v,source_to_physical=inv,source_faces=src.faces)
        for fidx,field in enumerate(fields):arrays.update({str(fidx)+'_'+k:z for k,z in field.items()})
        fp=path.with_name(path.stem+'-field.npz');atomic_npz(fp,**arrays);sg,_=self_contact_summary(mesh,max_pairs=250000,max_tests=2000000);assert sg['complete_within_pair_budget'];results.append(dict(partId=pid,source=bindings[pid],candidate=str(path),candidateSha256=digest_file(path),field=str(fp),fieldSha256=digest_file(fp),components=len(fields),metrics=metrics,dependencies=dependency,wholeSelfContact=sg,topology=topology(mesh)))
    supports=[]
    for row in results:
        src=Mesh.load(row['source']['mesh']);dst=Mesh.load(row['candidate']);changed=np.any(src.vertices!=dst.vertices,axis=1);ff=np.flatnonzero(changed[src.faces].any(1));points=np.r_[src.vertices[src.faces[ff]].reshape(-1,3),dst.vertices[dst.faces[ff]].reshape(-1,3)];supports.append(dict(partId=row['partId'],sourceFaceIds=ff,sourceAndCandidateBounds_m=np.array([points.min(0),points.max(0)]),maximumDisplacement_m=float(np.linalg.norm(dst.vertices-src.vertices,axis=1).max())))
    result=dict(schema='one-coupled-vertebrobasilar-proper-ring-C2-diagnostic/v6',sourceCatalogueSha256=EXPECTED,completeSourceInventoriesBeforeConstruction=[digest_file(BASE/'authority-inventory.json'),digest_file(BASE/'necessary-basilar-source-closure-v5/authority-inventory.json')],sourceGuideSha256=digest_file(BASE/'source-derived-throat-guides-independent-v5.json'),sourceClivalPlanesSha256=digest_file(BASE/'posterior-clivus-patch-qualified-rays-independent-v1.json'),movingOwners=results,changedSupport=supports,skeletonHeartAndAllOtherAnatomyFixed=True,sourceCalibreAndFacesAndUVPreserved=True,regionalPartitionsNotYetConstructed=True,finiteWallsAndBloodLumensNotEstablished=True,sourceInterfacesRemainUnmatchedOuterCasts=True,actualCompleteContextGatePending=True,acceptedForInstallation=False,elapsedSeconds=time.monotonic()-began)
    atomic_json(OUT/'evidence.json',json.loads(json.dumps(result,default=clean)));assert digest_file(CAT)==EXPECTED
    print(json.dumps(dict(evidence=str(OUT/'evidence.json'),elapsedSeconds=result['elapsedSeconds'],owners=len(results),selfCrossings=[(q['partId'],q['wholeSelfContact']['crossing_face_pairs'])for q in results],maximumDisplacement_m=max(q['maximumDisplacement_m']for q in supports))))

if __name__=='__main__':
    try:main()
    except Exception as e:
        if OUT.exists():atomic_json(OUT/'failed-stage.json',dict(error=type(e).__name__,message=str(e),notInstallable=True))
        raise
