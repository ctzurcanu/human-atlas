"""Sparse source-section clearance with finite native solve/motion limits."""
import time
import clarabel
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components
from .section_course_socp import _evidence_number


def solve_sparse_section_socp(stiffness_rows,matrix,offsets,*,maximum_motion_m,time_limit_seconds=5.,maximum_iterations=100,residual_m=1e-10):
    """Minimize ||d||²+||S d||² under A d>=b and station motion balls.

    The explicit ridge plus Gram construction is positive definite without a
    dense inverse/eigensolve. Distances are conditioned in millimetres. Stored
    sparse rows stay sparse; no dense face-by-station matrix is constructed.
    """
    tick=time.monotonic();S=sparse.csc_matrix(stiffness_rows,dtype=float);A=sparse.csc_matrix(matrix,dtype=float);b=np.asarray(offsets,float)
    count=S.shape[1]
    if (not count or A.shape!=(len(b),3*count)or b.ndim!=1 or not np.isfinite(S.data).all()or not np.isfinite(A.data).all()or not np.isfinite(b).all()
        or not np.isfinite(maximum_motion_m)or maximum_motion_m<=0 or not 0<time_limit_seconds<=30
        or not isinstance(maximum_iterations,int)or maximum_iterations<1 or not np.isfinite(residual_m)or residual_m<=0):
        raise ValueError('Finite sparse stiffness/clearance rows and bounded solve/motion required')
    H=sparse.eye(count,format='csc')+S.T@S
    P=sparse.triu(sparse.kron(H,sparse.eye(3)),format='csc')
    variable=np.arange(3*count);rows=(np.repeat(np.arange(count)*4,3)+np.tile([1,2,3],count))
    balls=sparse.coo_matrix((-np.ones(len(variable)),(rows,variable)),shape=(4*count,3*count)).tocsc()
    constraint=sparse.vstack([-A,balls],format='csc');rhs=np.r_[-b/.001,np.tile([maximum_motion_m/.001,0.,0.,0.],count)]
    cones=[clarabel.NonnegativeConeT(len(b))]+[clarabel.SecondOrderConeT(4)for _ in range(count)]
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.time_limit=time_limit_seconds;settings.max_iter=maximum_iterations
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-9
    result=clarabel.DefaultSolver(P,np.zeros(3*count),constraint,rhs,cones,settings).solve()
    proof=dict(method='Sparse ridge-plus-Gram source section SOCP',status=str(result.status),feasible=False,
        stationCount=count,halfspaces=len(b),inputNonzeroCoefficients=A.nnz,stiffnessNonzeroCoefficients=S.nnz,
        sparseConstraintBytes=constraint.data.nbytes+constraint.indices.nbytes+constraint.indptr.nbytes,
        denseFaceByStationMatrixAllocated=False,precisionEigenvalueLowerBound=1.,
        maximumMotionBound_m=maximum_motion_m,timeLimitSeconds=time_limit_seconds,iterations=int(result.iterations),
        solverTimeSeconds=_evidence_number(result.solve_time),solverPrimalObjective=_evidence_number(result.obj_val),
        solverDualObjective=_evidence_number(result.obj_val_dual),anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-tick)
    dual=np.asarray(result.z)[:len(b)];order=np.argsort(-np.abs(dual))[:min(20,len(b))]
    proof['leadingHalfspaceDualRows']=[dict(row=int(i),dual=_evidence_number(dual[i]))for i in order]
    if str(result.status)!='Solved':return None,proof
    d=np.asarray(result.x).reshape(count,3)*.001
    if not np.isfinite(d).all():proof['status']='NonfiniteSolution';return None,proof
    error=float(np.max(b-A@d.ravel(),initial=0.));motion=float(np.linalg.norm(d,axis=1).max())
    passed=error<=residual_m and motion<=maximum_motion_m+residual_m
    proof.update(feasible=bool(passed),maximumResidual_m=max(error,0.),maximumActualMotion_m=motion)
    return d if passed else None,proof


def solve_independent_sparse_sections(stiffness_rows,matrix,offsets,*,maximum_motion_m,time_limit_seconds=5.,wall_seconds=15.,residual_m=1e-10):
    """Split only algebraically independent objective/constraint components.

    This preserves the original feasible set and objective, unlike dividing an
    anatomical assembly by name. Every component has its own native deadline
    inside one finite overall deadline. Partial solves emit no whole candidate.
    """
    tick=time.monotonic();S=sparse.csr_matrix(stiffness_rows,dtype=float);A=sparse.csr_matrix(matrix,dtype=float);A.eliminate_zeros();b=np.asarray(offsets,float);count=S.shape[1]
    if not 0<wall_seconds<=30 or not count or A.shape!=(len(b),3*count)or not np.isfinite(A.data).all()or not np.isfinite(S.data).all()or not np.isfinite(b).all():
        raise ValueError('Finite sparse source rows and a finite overall deadline required')
    rows=np.repeat(np.arange(len(b)),np.diff(A.indptr));incidence=sparse.coo_matrix((np.ones(A.nnz),(rows,A.indices//3)),shape=(len(b),count)).tocsr();incidence.data[:]=1
    stiffness=S.copy();stiffness.eliminate_zeros();stiffness.data[:]=1
    graph=(incidence.T@incidence+stiffness.T@stiffness).tocsr();components,labels=connected_components(graph,directed=False)
    constant=np.diff(incidence.indptr)==0
    proof=dict(method='Exact sparse objective/halfspace component decomposition',componentCount=int(components),components=[],
        feasible=False,timeLimitPerComponentSeconds=time_limit_seconds,overallWallSeconds=wall_seconds,anatomicalAcceptance=False)
    if np.any(constant&(b>residual_m)):
        proof.update(status='ProtectedConstraintInfeasible',constantConflictRows=np.flatnonzero(constant&(b>residual_m)).tolist(),elapsedSeconds=time.monotonic()-tick)
        return None,proof
    result=np.zeros((count,3));witnesses=[]
    for label in range(components):
        remaining=wall_seconds-(time.monotonic()-tick)
        if remaining<.05:
            proof.update(status='ComponentWallDeadline',elapsedSeconds=time.monotonic()-tick);return None,proof
        ids=np.flatnonzero(labels==label);selected=np.flatnonzero(np.asarray(incidence[:,ids].getnnz(1)).ravel()>0);srows=np.flatnonzero(np.asarray(stiffness[:,ids].getnnz(1)).ravel()>0)
        columns=(ids[:,None]*3+np.arange(3)).ravel()
        assert np.all(np.isin(A[selected].indices//3,ids))and np.all(np.isin(S[srows].indices,ids))
        d,record=solve_sparse_section_socp(S[srows][:,ids],A[selected][:,columns],b[selected],maximum_motion_m=maximum_motion_m,
            time_limit_seconds=min(time_limit_seconds,remaining),residual_m=residual_m)
        proof['components'].append(dict(stations=ids.tolist(),halfspaces=len(selected),solver=record))
        witnesses.extend(dict(row=int(selected[w['row']]),dual=w['dual'],component=label)for w in record.get('leadingHalfspaceDualRows',[]))
        if d is None:
            proof.update(status=record['status'],leadingHalfspaceDualRows=witnesses,elapsedSeconds=time.monotonic()-tick);return None,proof
        result[ids]=d
    residual=float(np.max(b-A@result.ravel(),initial=0.));motion=float(np.linalg.norm(result,axis=1).max())
    passed=residual<=residual_m and motion<=maximum_motion_m+residual_m
    proof.update(status='Solved',feasible=bool(passed),maximumResidual_m=max(0.,residual),maximumActualMotion_m=motion,
        leadingHalfspaceDualRows=witnesses,elapsedSeconds=time.monotonic()-tick)
    return result if passed else None,proof
