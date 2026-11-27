"""Bounded detours of source sections around a closed forbidden envelope.

Every section moves by one translation; its polygonal cross-section is not
scaled. A quadratic course energy spreads the detour between neighboring
stations. This is an engineering proposal, never anatomical acceptance: the
caller must check the actual triangles, protected ports and branch junctions.
"""
import numpy as np


def closest_triangles(points,triangles,*,batch_size=48):
    """Exact exhaustive point/triangle distances, in bounded-memory batches."""
    p=np.asarray(points,float);t=np.asarray(triangles,float)
    if len(t)==0:raise ValueError('A real triangle envelope is required')
    a,b,c=t[:,0],t[:,1],t[:,2];ab=b-a;ac=c-a
    normal=np.cross(ab,ac)
    # The Gram determinant aa*cc-bb² cancels for valid sliver faces in exact
    # Boolean boundaries. Cross-product area retains their literal geometry.
    den=np.einsum('ij,ij->i',normal,normal)
    if np.any(den<=0):raise ValueError('Degenerate obstacle faces must be repaired first')
    output=[];indices=[]
    for start in range(0,len(p),batch_size):
        point=p[start:start+batch_size];ap=point[:,None]-a
        u=np.einsum('mfi,fi->mf',np.cross(ap,ac),normal)/den
        v=np.einsum('mfi,fi->mf',np.cross(ab,ap),normal)/den
        projected=a+u[:,:,None]*ab+v[:,:,None]*ac
        inside=(u>=0)&(v>=0)&(u+v<=1)
        best=np.where(inside[:,:,None],projected,np.inf)
        distance=np.where(inside,np.sum((projected-point[:,None])**2,axis=2),np.inf)
        for x,y in [(a,b),(b,c),(c,a)]:
            edge=y-x;l2=np.einsum('ij,ij->i',edge,edge)
            w=np.clip(np.einsum('mfi,fi->mf',point[:,None]-x,edge)/l2,0.,1.)
            candidate=x+w[:,:,None]*edge;d=np.sum((candidate-point[:,None])**2,axis=2)
            take=d<distance;best[take]=candidate[take];distance[take]=d[take]
        j=np.argmin(distance,axis=1);output.append(best[np.arange(len(point)),j]);indices.append(j)
    return np.concatenate(output),np.concatenate(indices)


def smooth_course_metric(stations,*,smoothing=64.,fixed=()):
    """Inverse quadratic metric, with explicit zero-motion source stations."""
    if stations<2:raise ValueError('At least two source sections are required')
    D=np.diff(np.eye(stations),n=2,axis=0);H=np.eye(stations)+float(smoothing)*(D.T@D)
    free=np.ones(stations,bool);free[list(fixed)]=False
    K=np.zeros((stations,stations));K[np.ix_(free,free)]=np.linalg.inv(H[np.ix_(free,free)])
    return K


def project_course_halfspaces(metric,stations,normals,offsets,*,initial=None,iterations=200,residual_m=1e-10):
    """Hildreth projection in a positive quadratic metric.

    A row means normal dot displacement[station] >= offset. Fixed stations
    have a zero metric row and cannot be silently released to make a fit pass.
    """
    K=np.asarray(metric,float);s=np.asarray(stations,int);n=np.asarray(normals,float);b=np.asarray(offsets,float)
    d=np.zeros((len(K),3))if initial is None else np.asarray(initial,float).copy()
    if len(s)==0:return d,dict(iterations=0,maximumResidualMeters=0.)
    if not np.allclose(np.linalg.norm(n,axis=1),1.,rtol=0.,atol=1e-12):raise ValueError('Unit constraint normals required')
    dual=np.zeros(len(s))
    for loop in range(iterations):
        for i,j in enumerate(s):
            gap=float(b[i]-n[i]@d[j]);diagonal=K[j,j]
            if diagonal<=0:
                if gap>residual_m:raise ValueError('A protected source station conflicts with the requested obstacle constraint')
                continue
            new=max(0.,dual[i]+gap/diagonal);delta=new-dual[i];dual[i]=new
            if delta:d+=delta*K[:,j,None]*n[i]
        residual=float(np.max(b-np.einsum('ij,ij->i',n,d[s])))
        if residual<=residual_m:return d,dict(iterations=loop+1,maximumResidualMeters=max(0.,residual))
    raise ValueError('Course constraints did not converge within the caller\'s iteration bound')


def propose_section_detour(positions,faces,rings,obstacle_triangles,*,contains,crossings,
                           origin,fixed_stations=(),smoothing=64.,clearance_m=.00001,
                           maximum_motion_m=.008,course_rounds=5,face_rounds=4,
                           center_constraint_stations=None):
    """A bounded proposal; callbacks provide independent inside/crossing tests.

    The complete source polygon is translated at each station. Explicit fixed
    stations stay literal, including storage. All detected face conflicts are
    added as actual obstacle-plane constraints after centerline initialization.
    Unresolved conflicts remain a rejection, never an increased tolerance.
    An optional center-constraint set limits this conservative initialization
    to a moving support whose fixed surroundings the caller separately checks.
    Actual whole-face refinement and all protected-coordinate checks remain.
    """
    v=np.asarray(positions,float);f=np.asarray(faces,int);obstacle=np.asarray(obstacle_triangles,float)
    origin=np.asarray(origin,float);fixed=set(map(int,fixed_stations));centers=np.array([v[q].mean(0)for q in rings]);offsets=np.array([v[q]-c for q,c in zip(rings,centers)])
    center_ids=None if center_constraint_stations is None else set(map(int,center_constraint_stations))
    if center_ids is not None and any(j<0 or j>=len(rings)for j in center_ids):raise ValueError('Center initialization station outside the source tube')
    station_of=np.full(len(v),-1,int)
    for j,ring in enumerate(rings):station_of[ring]=j
    if np.any(station_of[np.unique(f)]<0):raise ValueError('Every source face needs its actual source station')
    K=smooth_course_metric(len(rings),smoothing=smoothing,fixed=fixed);d=np.zeros_like(centers);stations=[];normals=[];bounds=[];rounds=[];face_evidence=[];rejected=None
    for attempt in range(course_rounds):
        current=centers+d;closest,_=closest_triangles(current,obstacle,batch_size=16)
        inside,ambiguous=contains(current)
        if np.any(ambiguous):raise ValueError('Ambiguous obstacle membership needs an explicit side')
        normal=current-closest;distance=np.linalg.norm(normal,axis=1)
        if np.any(distance<1e-12):
            rejected='A course station is on the obstacle boundary; an explicit side is required';break
        normal/=distance[:,None];normal[inside]*=-1
        support=-np.einsum('nki,ni->nk',offsets,normal).min(1);signed=np.where(inside,-distance,distance)
        active=np.flatnonzero(signed<support+clearance_m)
        if center_ids is not None:active=np.array([j for j in active if int(j)in center_ids],int)
        rounds.append(dict(attempt=attempt,activeStations=len(active),insideCenters=int(inside.sum())))
        if not len(active):break
        for j in active:
            stations.append(j);normals.append(normal[j]);bounds.append(float(support[j]+clearance_m-normal[j]@(centers[j]-closest[j])))
        try:d,proof=project_course_halfspaces(K,stations,normals,bounds)
        except ValueError as error:rejected=str(error);break
        rounds[-1]['projection']=proof
        if np.linalg.norm(d,axis=1).max()>maximum_motion_m:
            rejected='Course exceeds the declared motion domain';break
    def construct():
        candidate=v.copy();profile=[]
        for j,ring in enumerate(rings):
            if j in fixed:continue
            raw=v[ring]+d[j];stored=(raw-origin).astype('<f4').astype(float)+origin;candidate[ring]=stored
            before=v[ring]-v[ring].mean(0);after=stored-stored.mean(0);delta=float(np.linalg.norm(stored-raw,axis=1).max())
            if np.linalg.norm(after-before,axis=1).max()>2*delta+32*np.finfo(float).eps:raise ValueError('Stored section profile changed beyond actual serialization arithmetic')
            area=lambda p:float(np.linalg.norm(np.cross(p,np.roll(p,-1,axis=0)).sum(0))/2)
            area_before,area_after=area(before),area(after);perimeter=float(np.linalg.norm(before-np.roll(before,-1,axis=0),axis=1).sum());area_bound=2*delta*perimeter+len(ring)*delta**2+256*np.finfo(float).eps*area_before
            if abs(area_after-area_before)>area_bound:raise ValueError('Section area exceeds its measured storage bound')
            profile.append(dict(station=j,areaBeforeMetersSquared=area_before,areaAfterMetersSquared=area_after,storageProfileBoundMeters=2*delta,storageAreaBoundMetersSquared=area_bound))
        return candidate,profile
    candidate,profile=construct()
    if rejected is None:
        for attempt in range(face_rounds):
            hit=crossings(obstacle,candidate[f]);face_evidence.append(dict(attempt=attempt,properCrossings=hit['properOrCoplanarIntersections']))
            if not hit['examples']:break
            for obstacle_face,source_face in hit['examples']:
                tri=obstacle[obstacle_face];normal=np.cross(tri[1]-tri[0],tri[2]-tri[0]);normal/=np.linalg.norm(normal)
                for vertex in f[source_face]:
                    stations.append(int(station_of[vertex]));normals.append(normal);bounds.append(float(clearance_m-normal@(v[vertex]-tri[0])))
            try:d,proof=project_course_halfspaces(K,stations,normals,bounds)
            except ValueError as error:rejected=str(error);break
            face_evidence[-1]['projection']=proof
            if np.linalg.norm(d,axis=1).max()>maximum_motion_m:
                rejected='Actual face constraints exceed the declared motion domain';break
            candidate,profile=construct()
    protected=np.concatenate([rings[j]for j in sorted(fixed)])if fixed else np.zeros(0,int)
    if not np.array_equal(candidate[protected],v[protected]):raise ValueError('A protected port moved')
    return candidate,d,dict(initialization=rounds,actualFaceConstraints=face_evidence,rejectedReason=rejected,
        sectionProfiles=profile,protectedStations=sorted(fixed),sourceIndicesUnchanged=True,
        centerInitializationStations=None if center_ids is None else sorted(center_ids),
        maximumMotionMeters=float(np.linalg.norm(candidate-v,axis=1).max()),maximumAllowedMotionMeters=maximum_motion_m,
        prescribedClearanceMeters=clearance_m,anatomicallyAccepted=False)
