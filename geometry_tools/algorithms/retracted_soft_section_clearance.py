"""Bounded section clearance with hard context preservation and soft repair goals."""
import time
import clarabel
import numpy as np
from scipy import sparse


def prune_motion_ball_rows(matrix,rhs,maximum_motion_m):
    """Drop only halfspaces already guaranteed by the station motion balls."""
    A=sparse.csr_matrix(matrix);A.sum_duplicates();b=np.asarray(rhs,float)
    assert A.shape[1]%3==0 and A.shape[0]==len(b)and maximum_motion_m>0
    row=np.repeat(np.arange(len(b)),np.diff(A.indptr));count=A.shape[1]//3
    key=row*count+A.indices//3
    ids,inv=np.unique(key,return_inverse=True)
    lengths=np.sqrt(np.bincount(inv,weights=A.data**2,minlength=len(ids)))
    bound=np.bincount(ids//count,weights=lengths,minlength=len(b))*maximum_motion_m
    margin=64*np.finfo(float).eps*(1+abs(b)+bound)
    keep=b>-bound-margin
    return A[keep],b[keep],dict(originalRows=len(b),retainedRows=int(keep.sum()),droppedGuaranteedRows=int((~keep).sum()),
                               rule='Cauchy-Schwarz per independent station ball, with conservative Float64 margin',maximumMotion_m=maximum_motion_m)


def solve_soft_sections(stiffness,hard,hard_rhs,soft,soft_rhs,*,maximum_motion_m,soft_weight=10000.,time_limit_seconds=5.,distortion=None,distortion_bound=.075):
    tick=time.monotonic();S=sparse.csc_matrix(stiffness);A=sparse.csc_matrix(hard);B=sparse.csc_matrix(soft)
    b=np.asarray(hard_rhs,float);c=np.asarray(soft_rhs,float);count=S.shape[1];k=len(c);n=3*count
    assert count>0 and k>0 and A.shape==(len(b),n)and B.shape==(k,n)
    assert np.isfinite(S.data).all()and np.isfinite(A.data).all()and np.isfinite(B.data).all()and np.isfinite(b).all()and np.isfinite(c).all()
    assert 0<maximum_motion_m<=.01 and 1<=soft_weight<=1e6 and 0<time_limit_seconds<=10
    H=sparse.eye(count,format='csc')+S.T@S
    P=sparse.triu(sparse.block_diag([sparse.kron(H,sparse.eye(3)),sparse.eye(k)*soft_weight]),format='csc')
    hard_rows=sparse.hstack([-A,sparse.csc_matrix((len(b),k))])
    soft_rows=sparse.hstack([-B,-sparse.eye(k)])
    slack_rows=sparse.hstack([sparse.csc_matrix((k,n)),-sparse.eye(k)])
    variable=np.arange(n);rows=np.repeat(np.arange(count)*4,3)+np.tile([1,2,3],count)
    balls=sparse.coo_matrix((-np.ones(n),(rows,variable)),shape=(4*count,n+k)).tocsc()
    constraints=[hard_rows,soft_rows,slack_rows,balls]
    rhs=np.r_[-b/.001,-c/.001,np.zeros(k),np.tile([maximum_motion_m/.001,0,0,0],count)]
    cones=[clarabel.NonnegativeConeT(len(b)+2*k)]+[clarabel.SecondOrderConeT(4)for _ in range(count)]
    distortion_count=0
    if distortion is not None:
        C=sparse.coo_matrix(distortion);assert C.shape[1]==n and C.shape[0]%6==0 and 0<distortion_bound<1
        distortion_count=C.shape[0]//6
        mapped=(C.row//6)*7+1+C.row%6
        D=sparse.coo_matrix((-.001*C.data,(mapped,C.col)),shape=(distortion_count*7,n+k)).tocsc()
        constraints.append(D);rhs=np.r_[rhs,np.tile([distortion_bound,0,0,0,0,0,0],distortion_count)]
        cones += [clarabel.SecondOrderConeT(7)for _ in range(distortion_count)]
    M=sparse.vstack(constraints,format='csc')
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=time_limit_seconds;settings.max_iter=150
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-9
    result=clarabel.DefaultSolver(P,np.zeros(n+k),M,rhs,cones,settings).solve()
    proof=dict(method='Sparse rigid-section ridge/Gram QP with nonnegative soft repair slack and hard preservation halfspaces',
               status=str(result.status),feasible=False,stationCount=count,hardRows=len(b),softRows=k,softWeight=soft_weight,
               nativeSeconds=float(result.solve_time),timeLimitSeconds=time_limit_seconds,iterations=int(result.iterations),
               maximumMotionBound_m=maximum_motion_m,sourceFacetCones=distortion_count,elapsedSeconds=time.monotonic()-tick)
    # A reduced-accuracy optimum is not an accepted geometric result by its
    # status. Inspect its literal point against every original hard halfspace,
    # motion ball and facet cone below. Infeasible/time-limited native states
    # never emit a point; all actual geometric residual limits stay identical.
    if str(result.status)not in ['Solved','AlmostSolved']:return None,proof
    proof['optimalityCertified']=str(result.status)=='Solved'
    proof['independentPrimalFeasibilityRequired']=True
    d=np.asarray(result.x[:n]).reshape(count,3)*.001;s=np.asarray(result.x[n:])*.001
    from .origin_primal_retraction import retract_origin_feasible
    try:
        d,retraction=retract_origin_feasible(d,A,b,maximum_motion_m=maximum_motion_m,distortion=distortion,distortion_bound=distortion_bound)
    except ValueError as exc:
        proof['retractionRejected']=str(exc);return None,proof
    proof['originRetraction']=retraction
    if retraction['applied']:
        proof['optimalityCertified']=False
        s=np.maximum(c-B@d.ravel(),0.)
    residual=float(np.maximum(b-A@d.ravel(),0).max(initial=0))
    soft_residual=float(np.maximum(c-B@d.ravel()-s,0).max(initial=0))
    motion=float(np.linalg.norm(d,axis=1).max())
    passed=np.isfinite(d).all()and np.isfinite(s).all()and residual<=1e-10 and soft_residual<=1e-10 and s.min(initial=0)>=-1e-10 and motion<=maximum_motion_m+1e-10
    if distortion is not None:
        deformation=np.asarray(distortion@d.ravel()).reshape(-1,6)
        maximum=float(np.linalg.norm(deformation,axis=1).max())
        passed=passed and maximum<=distortion_bound+1e-9
        proof.update(maximumFacetTangentDistortion=maximum,facetTangentDistortionBound=distortion_bound)
    proof.update(feasible=bool(passed),maximumHardResidual_m=residual,maximumSoftResidual_m=soft_residual,
                 maximumMotion_m=motion,maximumRemainingSlack_m=float(s.max()),remainingSlackNorm_m=float(np.linalg.norm(s)),
                 softGoalsAreNotClearanceAcceptance=True,elapsedSeconds=time.monotonic()-tick)
    return d if passed else None,proof
