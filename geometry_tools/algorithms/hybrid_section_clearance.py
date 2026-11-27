"""Bounded section clearance with hard context preservation and soft repair goals."""
import time
import clarabel
import numpy as np
from scipy import sparse
from .spectral_tangent_cones import mixed_tangent_cones,mixed_tangent_residuals


def solve_hybrid_sections(stiffness,hard,hard_rhs,soft,soft_rhs,*,maximum_motion_m,soft_weight=10000.,time_limit_seconds=5.,distortion=None,distortion_bound=.075,numerical_cone_reserve=0.,distortion_offset=None):
    assert 0<=numerical_cone_reserve<=1e-5
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
    assert numerical_cone_reserve < distortion_bound
    if distortion is not None:
        C=sparse.coo_matrix(distortion);assert C.shape[1]==n and C.shape[0]%6==0 and 0<distortion_bound<1
        distortion_count=C.shape[0]//6
        offset=np.zeros((distortion_count,6)) if distortion_offset is None else np.asarray(distortion_offset,float).reshape(distortion_count,6)
        D,cone_rhs,extra_cones,spectral_mask,hybrid_proof=mixed_tangent_cones(C,offset,bound=distortion_bound,reserve=numerical_cone_reserve,extra_columns=k)
        constraints.append(D);rhs=np.r_[rhs,cone_rhs];cones+=extra_cones
    M=sparse.vstack(constraints,format='csc')
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=time_limit_seconds;settings.max_iter=150
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-9
    # Clarabel documented termination callback; no point from a time-limited
    # native run is eligible. Stop deliberately only at a tight primal point,
    # then independently replay literal halfspaces, balls and source cones.
    # https://clarabel.org/stable/examples/py/example_callback/
    history=[];callback_stopped=False
    def terminate(info):
        nonlocal callback_stopped
        row=dict(iteration=int(info.iterations),primalResidual=float(info.res_primal),dualResidual=float(info.res_dual),absoluteGap=float(info.gap_abs),relativeGap=float(info.gap_rel))
        history.append(row)
        callback_stopped=bool(info.iterations>=5 and info.res_primal<=1e-11 and info.res_dual<=1e-7 and (info.gap_abs<=1e-4 or info.gap_rel<=1e-6))
        return callback_stopped
    native=clarabel.DefaultSolver(P,np.zeros(n+k),M,rhs,cones,settings)
    native.set_termination_callback(terminate);result=native.solve()

    proof=dict(method='Sparse rigid-section ridge/Gram QP with nonnegative soft repair slack and hard preservation halfspaces',
               status=str(result.status),feasible=False,stationCount=count,hardRows=len(b),softRows=k,softWeight=soft_weight,
               nativeSeconds=float(result.solve_time),timeLimitSeconds=time_limit_seconds,iterations=int(result.iterations),
               maximumMotionBound_m=maximum_motion_m,sourceFacetCones=distortion_count,terminationHistory=history,intentionalPrimalFeasibilityTermination=callback_stopped,numericalConeReserve=numerical_cone_reserve,elapsedSeconds=time.monotonic()-tick)
    # A reduced-accuracy optimum is not an accepted geometric result by its
    # status. Inspect its literal point against every original hard halfspace,
    # motion ball and facet cone below. Infeasible/time-limited native states
    # never emit a point; all actual geometric residual limits stay identical.
    if str(result.status)not in ['Solved','AlmostSolved'] and not (str(result.status)=='CallbackTerminated' and callback_stopped):return None,proof
    proof['optimalityCertified']=str(result.status)=='Solved'
    proof['independentPrimalFeasibilityRequired']=True
    d=np.asarray(result.x[:n]).reshape(count,3)*.001;s=np.asarray(result.x[n:])*.001
    residual=float(np.maximum(b-A@d.ravel(),0).max(initial=0))
    soft_residual=float(np.maximum(c-B@d.ravel()-s,0).max(initial=0))
    motion=float(np.linalg.norm(d,axis=1).max())
    passed=np.isfinite(d).all()and np.isfinite(s).all()and residual<=1e-10 and soft_residual<=1e-10 and s.min(initial=0)>=-1e-10 and motion<=maximum_motion_m+1e-10
    if distortion is not None:
        deformation=np.asarray(distortion@d.ravel()).reshape(-1,6)+offset
        maximum=float(mixed_tangent_residuals(deformation,spectral_mask).max())
        proof.update(mixedTangentCones=hybrid_proof,maximumRawFrobeniusResidual=float(np.linalg.norm(deformation,axis=1).max()))
        passed=passed and maximum<=distortion_bound+1e-9
        proof.update(maximumOriginalCorotatedResidual=maximum,originalMetricNotReset=True,facetTangentDistortionBound=distortion_bound)
    proof.update(feasible=bool(passed),maximumHardResidual_m=residual,maximumSoftResidual_m=soft_residual,
                 maximumMotion_m=motion,maximumRemainingSlack_m=float(s.max()),remainingSlackNorm_m=float(np.linalg.norm(s)),
                 softGoalsAreNotClearanceAcceptance=True,elapsedSeconds=time.monotonic()-tick)
    return d if passed else None,proof
