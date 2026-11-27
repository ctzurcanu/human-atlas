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


def solve_soft_sections(stiffness,hard,hard_rhs,soft,soft_rhs,*,maximum_motion_m,soft_weight=10000.,time_limit_seconds=5.,distortion=None,distortion_bound=.075,numerical_cone_reserve=0.,distortion_offset=None,condition_hard_rows=False,regularizer_weight=1.,linear_section_objective=None,objective_gap_tolerance=None):
    assert 0<=numerical_cone_reserve<=1e-5
    tick=time.monotonic();S=sparse.csc_matrix(stiffness);A=sparse.csc_matrix(hard);B=sparse.csc_matrix(soft)
    b=np.asarray(hard_rhs,float);c=np.asarray(soft_rhs,float);count=S.shape[1];k=len(c);n=3*count
    assert count>0 and k>0 and A.shape==(len(b),n)and B.shape==(k,n)
    assert np.isfinite(S.data).all()and np.isfinite(A.data).all()and np.isfinite(B.data).all()and np.isfinite(b).all()and np.isfinite(c).all()
    assert 0<maximum_motion_m<=.01 and 1<=soft_weight<=1e6 and 0<time_limit_seconds<=10
    H=sparse.eye(count,format='csc')+S.T@S
    assert 0<=regularizer_weight<=1
    P=sparse.triu(sparse.block_diag([sparse.kron(H,sparse.eye(3))*regularizer_weight,sparse.eye(k)*soft_weight]),format='csc');P.eliminate_zeros()
    objective=np.zeros(n+k)
    if linear_section_objective is not None:
        v=np.asarray(linear_section_objective,float);assert v.shape==(count,3) and np.isfinite(v).all()
        objective[:n]=v.ravel()
    assert regularizer_weight>0 or linear_section_objective is not None
    # Positive row rescaling changes neither the halfspace nor the replay gate.
    # Near-contact rows can otherwise span many coefficient orders of magnitude.
    row_norm=np.sqrt(np.asarray(A.multiply(A).sum(axis=1)).ravel())
    row_scale=1/np.maximum(row_norm,1e-12) if condition_hard_rows else np.ones(len(b))
    native_A=sparse.diags(row_scale)@A
    hard_rows=sparse.hstack([-native_A,sparse.csc_matrix((len(b),k))])
    soft_rows=sparse.hstack([-B,-sparse.eye(k)])
    slack_rows=sparse.hstack([sparse.csc_matrix((k,n)),-sparse.eye(k)])
    variable=np.arange(n);rows=np.repeat(np.arange(count)*4,3)+np.tile([1,2,3],count)
    balls=sparse.coo_matrix((-np.ones(n),(rows,variable)),shape=(4*count,n+k)).tocsc()
    constraints=[hard_rows,soft_rows,slack_rows,balls]
    rhs=np.r_[-b*row_scale/.001,-c/.001,np.zeros(k),np.tile([maximum_motion_m/.001,0,0,0],count)]
    cones=[clarabel.NonnegativeConeT(len(b)+2*k)]+[clarabel.SecondOrderConeT(4)for _ in range(count)]
    distortion_count=0
    assert numerical_cone_reserve < distortion_bound
    if distortion is not None:
        C=sparse.coo_matrix(distortion);assert C.shape[1]==n and C.shape[0]%6==0 and 0<distortion_bound<1
        distortion_count=C.shape[0]//6
        mapped=(C.row//6)*7+1+C.row%6
        D=sparse.coo_matrix((-.001*C.data,(mapped,C.col)),shape=(distortion_count*7,n+k)).tocsc()
        offset=np.zeros((distortion_count,6)) if distortion_offset is None else np.asarray(distortion_offset,float).reshape(distortion_count,6)
        assert np.isfinite(offset).all()
        cone_rhs=np.c_[np.full(distortion_count,distortion_bound-numerical_cone_reserve),offset].ravel()
        constraints.append(D);rhs=np.r_[rhs,cone_rhs]
        cones += [clarabel.SecondOrderConeT(7)for _ in range(distortion_count)]
    M=sparse.vstack(constraints,format='csc')
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=time_limit_seconds;settings.max_iter=150
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-12 if condition_hard_rows else 1e-9
    if objective_gap_tolerance is not None:
        assert 0<objective_gap_tolerance<=1e-8
        settings.tol_gap_abs=settings.tol_gap_rel=objective_gap_tolerance
    if condition_hard_rows:
        settings.static_regularization_constant=1e-12
        settings.iterative_refinement_abstol=1e-14
        settings.iterative_refinement_reltol=1e-14
    result=clarabel.DefaultSolver(P,objective,M,rhs,cones,settings).solve()
    proof=dict(method='Sparse rigid-section ridge/Gram QP with nonnegative soft repair slack and hard preservation halfspaces',
               status=str(result.status),feasible=False,stationCount=count,hardRows=len(b),softRows=k,softWeight=soft_weight,
               nativeSeconds=float(result.solve_time),timeLimitSeconds=time_limit_seconds,iterations=int(result.iterations),
               maximumMotionBound_m=maximum_motion_m,sourceFacetCones=distortion_count,numericalConeReserve=numerical_cone_reserve,elapsedSeconds=time.monotonic()-tick)
    proof.update(positiveHardRowNormalization=condition_hard_rows,nativeFeasibilityTolerance=settings.tol_feas,
                 nativeObjectiveGapTolerance=settings.tol_gap_abs,
                 originalUnscaledResidualsRequired=True,regularizerWeight=regularizer_weight,
                 nonzeroLinearObjectiveEntries=int(np.count_nonzero(objective)))
    for field in ['obj_val','obj_val_dual','r_prim','r_dual']:
        if hasattr(result,field):proof['native_'+field]=float(getattr(result,field))
    # A reduced-accuracy optimum is not an accepted geometric result by its
    # status. Inspect its literal point against every original hard halfspace,
    # motion ball and facet cone below. Infeasible/time-limited native states
    # never emit a point; all actual geometric residual limits stay identical.
    if str(result.status)not in ['Solved','AlmostSolved']:return None,proof
    proof['optimalityCertified']=str(result.status)=='Solved'
    proof['independentPrimalFeasibilityRequired']=True
    d=np.asarray(result.x[:n]).reshape(count,3)*.001;s=np.asarray(result.x[n:])*.001
    residual=float(np.maximum(b-A@d.ravel(),0).max(initial=0))
    soft_residual=float(np.maximum(c-B@d.ravel()-s,0).max(initial=0))
    motion=float(np.linalg.norm(d,axis=1).max())
    passed=np.isfinite(d).all()and np.isfinite(s).all()and residual<=1e-10 and soft_residual<=1e-10 and s.min(initial=0)>=-1e-10 and motion<=maximum_motion_m+1e-10
    if distortion is not None:
        deformation=np.asarray(distortion@d.ravel()).reshape(-1,6)+offset
        maximum=float(np.linalg.norm(deformation,axis=1).max())
        passed=passed and maximum<=distortion_bound+1e-9
        proof.update(maximumOriginalCorotatedResidual=maximum,originalMetricNotReset=True,facetTangentDistortionBound=distortion_bound)
    proof.update(feasible=bool(passed),maximumHardResidual_m=residual,maximumSoftResidual_m=soft_residual,
                 maximumMotion_m=motion,maximumRemainingSlack_m=float(s.max()),remainingSlackNorm_m=float(np.linalg.norm(s)),
                 softGoalsAreNotClearanceAcceptance=True,elapsedSeconds=time.monotonic()-tick)
    if len(b):
        worst=int(np.argmax(b-A@d.ravel()))
        proof.update(worstOriginalHardRow=worst,worstOriginalHardRowNorm=float(row_norm[worst]),
                     maximumNormalizedHardResidual_m=float(np.maximum((b-A@d.ravel())*row_scale,0).max(initial=0)))
    return d if passed else None,proof
