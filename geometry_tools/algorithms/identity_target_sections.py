"""Exact quadratic compilation of paired identity soft targets.

For [I;-I] targets [t;-t], the two nonnegative hinge slacks have squared sum
||x-t||^2. Eliminating them preserves the entire optimization problem and
removes six slack variables/section. Other target matrices are rejected.
"""
import time
import clarabel
import numpy as np
from scipy import sparse


def solve_identity_target_sections(stiffness,hard,hard_rhs,soft,soft_rhs,*,
        maximum_motion_m,soft_weight=10000.,time_limit_seconds=2.,**kwargs):
    if any(v is not None for k,v in kwargs.items() if k in ['distortion','distortion_offset']):
        raise ValueError('Source-shape limits must use the independently replayed outer stages')
    tick=time.monotonic();S=sparse.csc_matrix(stiffness);count=S.shape[1];n=3*count
    A=sparse.csc_matrix(hard);b=np.asarray(hard_rhs,float);B=sparse.csr_matrix(soft);c=np.asarray(soft_rhs,float)
    identity=sparse.eye(n,format='csr')
    if B.shape!=(2*n,n) or c.shape!=(2*n,) or (B[:n]-identity).nnz or (B[n:]+identity).nnz or not np.array_equal(c[n:],-c[:n]):
        raise ValueError('Literal paired identity targets required for exact slack elimination')
    if A.shape!=(len(b),n) or not np.isfinite(np.r_[A.data,b,c,S.data]).all():
        raise ValueError('Finite matching original hard constraints required')
    if not 0<maximum_motion_m<=.01 or not 0<time_limit_seconds<=10 or not 1<=soft_weight<=1e6:
        raise ValueError('Finite explicit pose, target weight and native wall limits required')
    H=sparse.eye(count,format='csc')+S.T@S
    P=sparse.triu(sparse.kron(H,sparse.eye(3))+soft_weight*sparse.eye(n),format='csc')
    q=-soft_weight*c[:n]/.001
    variable=np.arange(n);rows=np.repeat(np.arange(count)*4,3)+np.tile([1,2,3],count)
    balls=sparse.coo_matrix((-np.ones(n),(rows,variable)),shape=(4*count,n)).tocsc()
    M=sparse.vstack([-A,balls],format='csc')
    rhs=np.r_[-b/.001,np.tile([maximum_motion_m/.001,0,0,0],count)]
    cones=[clarabel.NonnegativeConeT(len(b))]+[clarabel.SecondOrderConeT(4) for _ in range(count)]
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=time_limit_seconds
    settings.max_iter=150;settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-9
    result=clarabel.DefaultSolver(P,q,M,rhs,cones,settings).solve()
    proof=dict(complete=True,method='Exact paired-identity slack elimination',status=str(result.status),feasible=False,
        stationCount=count,hardRows=len(b),softRows=2*n,eliminatedSlackVariables=2*n,
        eliminatedNonnegativeConstraints=4*n,objectiveAndPhysicalBoundsUnchanged=True,
        nativeSeconds=float(result.solve_time),timeLimitSeconds=time_limit_seconds,iterations=int(result.iterations),
        anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-tick)
    if str(result.status) not in ['Solved','AlmostSolved']:
        return None,proof
    d=np.asarray(result.x).reshape(count,3)*.001
    residual=float(np.maximum(b-A@d.ravel(),0).max(initial=0))
    motion=float(np.linalg.norm(d,axis=1).max())
    feasible=np.isfinite(d).all() and residual<=1e-10 and motion<=maximum_motion_m+1e-10
    proof.update(feasible=bool(feasible),maximumHardResidual_m=residual,maximumMotion_m=motion,
        squaredTargetDeviation_m2=float(np.linalg.norm(d.ravel()-c[:n])**2),optimalityCertified=str(result.status)=='Solved',
        independentPrimalFeasibilityRequired=True,elapsedSeconds=time.monotonic()-tick)
    return d if feasible else None,proof
