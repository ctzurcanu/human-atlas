"""Bound one station's attainable motion and retain the limiting dual weights.

This diagnoses a caller's current convex fitting model. It is not a global
anatomical impossibility proof or an automatically installable mesh proposal.
"""
import time
import clarabel
import numpy as np
from scipy import sparse


def diagnose_direction(hard,rhs,station_count,station,direction,*,maximum_motion_m,
                       distortion=None,distortion_bound=.075,time_limit_seconds=5.):
    started=time.monotonic();A=sparse.csc_matrix(hard);b=np.asarray(rhs,float)
    n=3*station_count;g=np.asarray(direction,float)
    if A.shape!=(len(b),n) or g.shape!=(3,) or not 0<=station<station_count:
        raise ValueError('Explicit source station, direction and corresponding hard system required')
    if not np.isfinite(np.r_[A.data,b,g]).all() or np.linalg.norm(g)<=0:
        raise ValueError('Finite source constraints and nonzero direction required')
    if not 0<maximum_motion_m<=.01 or not 0<time_limit_seconds<=10:
        raise ValueError('Finite motion and native runtime limits required')
    g=g/np.linalg.norm(g);q=np.zeros(n);q[3*station:3*station+3]=-g
    variable=np.arange(n);rows=np.repeat(np.arange(station_count)*4,3)+np.tile([1,2,3],station_count)
    balls=sparse.coo_matrix((-np.ones(n),(rows,variable)),shape=(4*station_count,n)).tocsc()
    constraints=[-A,balls];bounds=np.r_[-b/.001,np.tile([maximum_motion_m/.001,0,0,0],station_count)]
    cones=[clarabel.NonnegativeConeT(len(b))]+[clarabel.SecondOrderConeT(4) for _ in range(station_count)]
    C=None;face_count=0
    if distortion is not None:
        C=sparse.coo_matrix(distortion)
        if C.shape[1]!=n or C.shape[0]%6 or not 0<distortion_bound<1:raise ValueError('Corresponding six-row facet chart constraints required')
        face_count=C.shape[0]//6;mapped=C.row//6*7+1+C.row%6
        D=sparse.coo_matrix((-.001*C.data,(mapped,C.col)),shape=(7*face_count,n)).tocsc()
        constraints.append(D);bounds=np.r_[bounds,np.tile([distortion_bound,0,0,0,0,0,0],face_count)]
        cones.extend(clarabel.SecondOrderConeT(7) for _ in range(face_count))
    M=sparse.vstack(constraints,format='csc');P=sparse.eye(n,format='csc')*1e-9
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.max_iter=150;settings.time_limit=time_limit_seconds
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-9
    result=clarabel.DefaultSolver(P,q,M,bounds,cones,settings).solve()
    report=dict(complete=True,status=str(result.status),nativeSeconds=float(result.solve_time),
                timeLimitSeconds=time_limit_seconds,feasible=False,optimalityCertified=False,
                anatomicalAcceptance=False,globalRoutingImpossibilityEstablished=False)
    if str(result.status)!='Solved':return None,report
    x=np.asarray(result.x);d=x.reshape(station_count,3)*.001;dual=np.asarray(result.z)
    residual=float(np.maximum(b-A@d.ravel(),0).max(initial=0));motion=float(np.linalg.norm(d,axis=1).max())
    chart=float(np.linalg.norm((C@d.ravel()).reshape(-1,6),axis=1).max()) if C is not None else 0.
    assert residual<=1e-10 and motion<=maximum_motion_m+1e-10 and chart<=distortion_bound+1e-9
    stationarity=float(abs(P@x+q+M.T@dual).max())
    assert stationarity<1e-6
    start=len(b)+4*station_count
    report.update(feasible=True,optimalityCertified=True,station=int(station),direction=g.tolist(),
        attainedDirectionalMotion_m=float(g@d[station]),maximumHardResidual_m=residual,
        maximumMotion_m=motion,maximumFacetTangentDistortion=chart,
        dualStationarityResidual=stationarity,primalObjective=float(result.obj_val),dualObjective=float(result.obj_val_dual),
        hardRowDualWeights=dual[:len(b)].tolist(),
        facetConeDualWeights=dual[start::7].tolist() if face_count else [],
        elapsedSeconds=time.monotonic()-started)
    return d,report
