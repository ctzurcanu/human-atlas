"""Bounded sparse quadratic fit under actual contact/face halfspaces.

All planes are feasible at zero motion. The source-cell graph regularizes
neighbouring translations, avoiding the isolated kinks of sparse L1 changes.
Explicit separating planes and linearized face guards remain hard. Exact
emitted triangles remain a gate. OSQP uses an explicit native time limit.
"""
import time
import numpy as np
from scipy import sparse
import osqp


def _compact_linear_rows(matrix,lower,upper):
    """Remove only exact duplicate inequalities and vacuous zero rows."""
    A=matrix.tocsr();A.sum_duplicates();A.eliminate_zeros();A.sort_indices()
    keep=[];seen=set();zero=0;duplicates=0
    for i in range(A.shape[0]):
        start,stop=A.indptr[i:i+2]
        if start==stop:
            if not lower[i]<=0<=upper[i]:raise ValueError('A fixed source constraint is infeasible')
            zero+=1;continue
        key=(float(lower[i]),float(upper[i]),A.indices[start:stop].tobytes(),A.data[start:stop].tobytes())
        if key in seen:duplicates+=1;continue
        seen.add(key);keep.append(i)
    ids=np.asarray(keep,int)
    return A[ids].tocsc(),lower[ids],upper[ids],dict(inputRows=A.shape[0],retainedRows=len(ids),
        removedVacuousFixedRows=zero,removedExactDuplicateRows=duplicates,tolerancesChanged=False)


def project_cell_direction(vertices,faces,owner,fixed,direction,contact_owner,normal,rhs,*,step_m=.002,wall_seconds=2.,objective_matrix=None,objective_rhs=None,solver_backend='osqp'):
    start=time.monotonic();v=np.asarray(vertices,float);f=np.asarray(faces,int);o=np.asarray(owner,int)
    fixed=np.asarray(fixed,bool);d=np.asarray(direction,float);live=np.flatnonzero(~fixed)
    compact=np.full(len(fixed),-1,int);compact[live]=np.arange(len(live));size=3*len(live)
    if not 0<wall_seconds<=3 or np.any(np.asarray(rhs)>1e-9):raise ValueError('Finite projection time and zero-feasible contact planes required')
    rows=[];cols=[];values=[]
    for end,sign in [(0,1.),(1,-1.)]:
        group=contact_owner[:,end];valid=~fixed[group]
        for axis in range(3):
            rows.extend(np.flatnonzero(valid).tolist());cols.extend((compact[group[valid]]*3+axis).tolist());values.extend((sign*normal[valid,axis]).tolist())
    contacts=sparse.coo_matrix((values,(rows,cols)),shape=(len(rhs),size)).tocsr()
    tri=v[f];cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1)
    if np.any(area<1e-20):raise ValueError('No orientation plane for a degenerate source face')
    n=cross/area[:,None];g1=np.cross(tri[:,2]-tri[:,0],n);g2=np.cross(n,tri[:,1]-tri[:,0]);g0=-g1-g2
    gradients=np.stack([g0,g1,g2],axis=1);unit=np.linalg.norm(gradients,axis=2).max(1)
    rows=[];cols=[];values=[]
    for corner in range(3):
        group=o[f[:,corner]];valid=~fixed[group]
        for axis in range(3):
            rows.extend(np.flatnonzero(valid).tolist());cols.extend((compact[group[valid]]*3+axis).tolist());values.extend((-gradients[valid,corner,axis]/unit[valid]).tolist())
    orientation=sparse.coo_matrix((values,(rows,cols)),shape=(len(f),size)).tocsr()
    # A completely fixed face, or a face wholly inside one translating rigid
    # cell, has exactly invariant area/orientation. Its derivative row is
    # vacuous; avoid thousands of numerically redundant solver constraints.
    active_face=np.any(~fixed[o[f]],axis=1)&np.any(o[f]!=o[f[:,0,None]],axis=1)
    orientation=orientation[active_face]
    contact_bounds=-np.asarray(rhs);orientation_bounds=(.75*area/unit)[active_face]
    ident=sparse.eye(size,format='csr')
    A=sparse.vstack([contacts,orientation,ident],format='csc')
    upper=np.r_[contact_bounds,orientation_bounds,np.full(size,step_m)]/.001
    lower=np.r_[np.full(len(rhs)+int(active_face.sum()),-np.inf),np.full(size,-step_m)]/.001
    A,lower,upper,row_proof=_compact_linear_rows(A,lower,upper)
    if objective_matrix is None:
        objective_matrix=ident;objective_rhs=d[live].ravel()
    objective=sparse.csr_matrix(objective_matrix);target=np.asarray(objective_rhs,float)
    if objective.shape[1]!=size or target.shape!=(objective.shape[0],):raise ValueError('Source-bound sparse least-squares objective required')
    P=sparse.triu(objective.T@objective,format='csc');linear=-(objective.T@target)/.001
    if solver_backend=='osqp':
        solver=osqp.OSQP();solver.setup(P=P,q=linear,A=A,l=lower,u=upper,verbose=False,
            time_limit=wall_seconds,max_iter=6000,eps_abs=1e-7,eps_rel=1e-7,polishing=True,
            adaptive_rho=True,check_termination=10)
        result=solver.solve(raise_error=False)
        status=result.info.status;iterations=int(result.info.iter);values=result.x
        solver_name,version,iteration_limit='OSQP',osqp.__version__,6000
        solved=status=='solved'
    elif solver_backend=='clarabel':
        import clarabel
        finite=np.isfinite(lower)
        conic_A=sparse.vstack([A,-A[finite]],format='csc')
        conic_b=np.r_[upper,-lower[finite]]
        settings=clarabel.DefaultSettings();settings.verbose=False
        settings.time_limit=wall_seconds;settings.max_iter=100
        settings.tol_gap_abs=settings.tol_gap_rel=settings.tol_feas=1e-9
        settings.tol_infeas_abs=settings.tol_infeas_rel=1e-9
        result=clarabel.DefaultSolver(P,linear,conic_A,conic_b,
            [clarabel.NonnegativeConeT(len(conic_b))],settings).solve()
        status=str(result.status);iterations=int(result.iterations);values=result.x
        solver_name,version,iteration_limit='Clarabel',clarabel.__version__,100
        solved=status=='Solved'
    else:raise ValueError('Explicit supported QP backend required: osqp or clarabel')
    if not solved:raise ValueError('Hard contact quadratic fit did not finish: '+status)
    q=np.zeros_like(d);q[live]=np.asarray(values).reshape(-1,3)*.001
    maximum=np.linalg.norm(q,axis=1).max(initial=0.)
    if maximum>step_m:q*=step_m/maximum
    residual=contacts@q[live].ravel()-contact_bounds
    if residual.max(initial=0)>2e-9:raise ValueError('Projected direction violates its actual contact plane')
    return q,dict(status=status,solver=solver_name,solverVersion=version,iterations=iterations,nativeIterationLimit=iteration_limit,hardContactRows=len(rhs),
        maximumContactPlaneViolation_m=float(residual.max(initial=0)),fixedCellsLiteral=bool(np.all(q[fixed]==0)),
        sourceFaceOrientationRowsExcluded=int((~active_face).sum()),linearConstraints=row_proof,
        elapsedSeconds=time.monotonic()-start,nativeWallSeconds=wall_seconds,anatomicalAcceptance=False)
