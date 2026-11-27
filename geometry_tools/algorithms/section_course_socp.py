"""Weighted source-section halfspaces with Euclidean motion bounds.

The bound is inside the conic optimization, rather than checked after an
unbounded smooth-course projection. This is a local linearized constraint
solve; it never certifies the resulting anatomical geometry.
"""
import time

import clarabel
import numpy as np
from scipy import sparse


def _evidence_number(value):
    """Unavailable solver values stay explicit and serializable on failure."""
    value=float(value)
    return value if np.isfinite(value) else None


def solve_weighted_section_socp(metric, weights, normals, offsets, *,
                                maximum_motion_m=None, residual_m=1e-10,
                                time_limit_seconds=2., maximum_iterations=100,
                                precision_matrix=None, vector_rows=None, vector_offsets=None,
                                target_displacement=None,deduplicate_halfspaces=True,
                                elastic_row_groups=None,elastic_penalty=1000.,
                                motion_weights=None):
    """Return displacement and evidence, or ``None`` on infeasibility/failure.

    Rows mean ``normal[i] @ (weights[i] @ displacement) >= offsets[i]``.
    Positive diagonal metric rows are movable; zero rows remain exactly zero.
    With a finite ``maximum_motion_m``, minimize the same quadratic energy as
    the inverse supplied metric, subject to a ball at every movable station.
    With ``None``, instead minimize a common radius, to diagnose the smallest
    possible maximum motion under these particular linearized halfspaces.
    All internal distances are millimetres for numerical conditioning.
    ``motion_weights`` maps parameter vectors to the actual source stations;
    their norms, rather than arbitrary basis-coordinate norms, are bounded.
    """
    start = time.monotonic()
    K = np.asarray(metric, float)
    w = np.asarray(weights, float)
    n = np.asarray(normals, float)
    b = np.asarray(offsets, float)
    if K.ndim != 2 or K.shape[0] != K.shape[1] or not len(K):
        raise ValueError('A nonempty square inverse course metric is required')
    if w.shape != (len(b), len(K)) or n.shape != (len(b), 3):
        raise ValueError('One section-weight row and unit normal per bound required')
    if not all(np.isfinite(x).all() for x in [K, w, n, b]):
        raise ValueError('All metric and halfspace values must be finite')
    if not np.allclose(K, K.T, rtol=0., atol=1e-12):
        raise ValueError('The inverse course metric must be symmetric')
    if not np.allclose(np.linalg.norm(n, axis=1), 1., rtol=0., atol=1e-12):
        raise ValueError('Unit constraint normals are required')
    if not np.isfinite(residual_m) or residual_m <= 0.:
        raise ValueError('A positive finite residual bound is required')
    if not np.isfinite(time_limit_seconds) or not 0. < time_limit_seconds <= 30.:
        raise ValueError('A finite solver deadline between zero and 30 seconds is required')
    if not isinstance(maximum_iterations, int) or maximum_iterations < 1:
        raise ValueError('A positive finite iteration limit is required')
    if maximum_motion_m is not None and (
            not np.isfinite(maximum_motion_m) or maximum_motion_m <= 0.):
        raise ValueError('A positive finite motion radius is required')
    free = np.flatnonzero(np.diag(K) > 0.)
    fixed = np.flatnonzero(np.diag(K) <= 0.)
    if np.any(K[fixed] != 0.) or np.any(K[:, fixed] != 0.):
        raise ValueError('Protected metric rows and columns must be literally zero')
    inverse = K[np.ix_(free, free)]
    explicit_precision = None
    if precision_matrix is not None:
        H = np.asarray(precision_matrix, float)
        if len(fixed) or H.shape != K.shape or not np.isfinite(H).all():
            raise ValueError('Explicit precision requires a finite all-free matching matrix')
        if not np.allclose(H, H.T, rtol=0., atol=1e-12) or not np.allclose(H@K, np.eye(len(K)), rtol=0., atol=1e-8):
            raise ValueError('Explicit precision does not invert the declared course metric')
        explicit_precision = H
    if len(free):
        np.linalg.cholesky(inverse)
    motion = np.eye(len(K))[free] if motion_weights is None else np.asarray(motion_weights,float)
    if motion.ndim!=2 or motion.shape[1]!=len(K) or not np.isfinite(motion).all():
        raise ValueError('Finite motion weights need one coefficient per parameter station')
    if len(free)and np.linalg.matrix_rank(motion[:,free])!=len(free):
        raise ValueError('Motion weights must bound every movable parameter direction')
    scale = .001
    count = 3 * len(free)
    minimize_radius = maximum_motion_m is None
    groups=np.full(len(b),-1,int) if elastic_row_groups is None else np.asarray(elastic_row_groups)
    if groups.shape!=(len(b),) or groups.dtype.kind not in 'iu' or np.any(groups < -1):
        raise ValueError('One integer elastic group or -1 hard marker per halfspace required')
    group_labels=np.unique(groups[groups>=0])
    if len(group_labels) and (minimize_radius or not np.isfinite(elastic_penalty) or elastic_penalty<=0):
        raise ValueError('Elastic route proposals require finite motion and positive penalty')
    group_rank={int(g):i for i,g in enumerate(group_labels)}
    elastic_count=len(group_labels)
    if target_displacement is not None and minimize_radius:
        raise ValueError('A shape target requires an explicit finite motion bound')
    variables = count + int(minimize_radius) + elastic_count
    section_A = (w[:, free, None] * n[:, None]).reshape(len(b), count)
    row_support = np.sum(np.abs(w[:, free]), axis=1)
    if vector_rows is not None:
        extra = np.asarray(vector_rows, float)
        eb = np.asarray(vector_offsets, float)
        if extra.shape != (len(eb), len(K), 3) or not np.isfinite(extra).all() or not np.isfinite(eb).all():
            raise ValueError('Finite vector constraint rows need one 3D coefficient per control')
        section_A = np.r_[section_A, extra[:, free].reshape(len(eb), count)]
        row_support = np.r_[row_support, np.linalg.norm(extra[:, free], axis=2).sum(1)]
        b = np.r_[b, eb]
        groups=np.r_[groups,np.full(len(eb),-1,int)]
    if motion_weights is not None and len(free):
        source_coefficients=np.einsum('smj,mn->snj',section_A.reshape(len(b),len(free),3),np.linalg.pinv(motion[:,free]))
        row_support=np.linalg.norm(source_coefficients,axis=2).sum(1)
    lower = np.divide(np.maximum(b, 0.), row_support,
                      out=np.full(len(b), np.inf), where=row_support > 0.)
    lower[(row_support == 0.) & (b <= 0.)] = 0.
    protected_conflicts = np.flatnonzero((row_support == 0.) & (b > residual_m) & (groups<0))
    proof = dict(method='Clarabel weighted halfspaces with per-station second-order cones',
                 minimizeCommonRadius=minimize_radius, fixedStations=fixed.tolist(),
                 freeStations=free.tolist(), requestedMotionBoundMeters=maximum_motion_m,
                 constraintRowNecessaryRadiusLowerBoundMeters=[_evidence_number(x)for x in lower],
                 necessaryRadiusLowerBoundMeters=_evidence_number(lower.max(initial=0.)),
                 infiniteNecessaryRadius=bool(np.isinf(lower).any()),
                 timeLimitSeconds=time_limit_seconds, residualBoundMeters=residual_m,
                 geometryModified=False, anatomicalAcceptance=False)
    proof.update(elasticProposal=bool(elastic_count),elasticGroups=group_labels.tolist())
    proof.update(mappedSourceMotion=motion_weights is not None,controlledMotionRows=len(motion))
    if len(protected_conflicts):
        proof.update(status='ProtectedConstraintInfeasible', feasible=False,
                     protectedConflictRows=protected_conflicts.tolist(),
                     elapsedSeconds=time.monotonic() - start)
        return None, proof
    if not variables or (not len(b) and target_displacement is None):
        proof.update(status='Solved', feasible=True, iterations=0,
                     minimumRadiusMeters=0., maximumActualMotionMeters=0.,
                     maximumResidualMeters=0., elapsedSeconds=time.monotonic() - start)
        return np.zeros((len(K), 3)), proof
    # Source aliases and the two incident triangles of an edge often submit
    # the identical finite constraint repeatedly. Remove literal duplicates
    # only; retain every original row for the independent residual replay and
    # map solver witnesses back to their actual input-row identity.
    active_rows=np.arange(len(b))
    if deduplicate_halfspaces and len(b):
        _,first=np.unique(np.c_[section_A,b,groups],axis=0,return_index=True)
        active_rows=np.sort(first)
    solve_A=section_A[active_rows];solve_b=b[active_rows]
    proof.update(inputHalfspaceRows=len(b),uniqueHalfspaceRows=len(active_rows),
                 exactDuplicateRowsRemoved=len(b)-len(active_rows))
    blocks = []
    rhs = []
    cones = []
    first = sparse.csc_matrix(-solve_A)
    if minimize_radius:
        first = sparse.hstack([first, sparse.csc_matrix((len(solve_b), 1))], format='csc')
    elif elastic_count:
        chosen=groups[active_rows];soft=np.flatnonzero(chosen>=0)
        extra=sparse.csc_matrix((-np.ones(len(soft)),(soft,[group_rank[int(chosen[i])] for i in soft])),shape=(len(solve_b),elastic_count))
        first=sparse.hstack([first,extra],format='csc')
    blocks.append(first)
    rhs.append(-solve_b / scale)
    cones.append(clarabel.NonnegativeConeT(len(solve_b)))
    if elastic_count:
        blocks.append(sparse.hstack([sparse.csc_matrix((elastic_count,count)),-sparse.eye(elastic_count)],format='csc'))
        rhs.append(np.zeros(elastic_count));cones.append(clarabel.NonnegativeConeT(elastic_count))
    if minimize_radius:
        nonnegative = sparse.csc_matrix(([-1.], ([0], [count])), shape=(1, variables))
        blocks.append(nonnegative)
        rhs.append(np.zeros(1))
        cones.append(clarabel.NonnegativeConeT(1))
    for coefficients in motion[:,free]:
        columns=[];rows=[];values=[]
        for j in np.flatnonzero(coefficients):
            columns.extend([3*int(j),3*int(j)+1,3*int(j)+2]);rows.extend([1,2,3]);values.extend([-float(coefficients[j])]*3)
        if minimize_radius:
            columns.append(count)
            rows.append(0)
            values.append(-1.)
        blocks.append(sparse.csc_matrix((values, (rows, columns)), shape=(4, variables)))
        rhs.append(np.r_[0. if minimize_radius else maximum_motion_m / scale, np.zeros(3)])
        cones.append(clarabel.SecondOrderConeT(4))
    if minimize_radius:
        P = sparse.csc_matrix((variables, variables))
        q = np.r_[np.zeros(count), 1.]
    else:
        H = np.linalg.inv(inverse) if explicit_precision is None else explicit_precision
        H = (H + H.T) / 2.
        P = sparse.triu(sparse.kron(H, sparse.eye(3)), format='csc')
        if elastic_count:P=sparse.block_diag([P,elastic_penalty*sparse.eye(elastic_count)],format='csc')
        q = np.zeros(variables)
        if target_displacement is not None:
            target = np.asarray(target_displacement, float)
            if target.shape != (len(K), 3) or not np.isfinite(target).all():
                raise ValueError('Finite target displacement required for every source control')
            q[:count] = -(H@target[free]/scale).ravel()
    settings = clarabel.DefaultSettings()
    settings.verbose = False
    settings.time_limit = time_limit_seconds
    settings.max_iter = maximum_iterations
    settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-9
    settings.tol_infeas_abs = settings.tol_infeas_rel = 1e-9
    solution = clarabel.DefaultSolver(P, q, sparse.vstack(blocks, format='csc'),
                                     np.concatenate(rhs), cones, settings).solve()
    status = str(solution.status)
    dual = np.asarray(solution.z)[:len(solve_b)]
    order = np.argsort(-np.abs(dual))[:min(12, len(solve_b))]
    proof.update(status=status, solverIterations=int(solution.iterations),
                 solverTimeSeconds=float(solution.solve_time),
                 solverPrimalObjective=_evidence_number(solution.obj_val),
                 solverDualObjective=_evidence_number(solution.obj_val_dual),
                 nonfiniteSolverEvidence=bool(not np.isfinite([solution.obj_val,solution.obj_val_dual]).all()or not np.isfinite(dual).all()),
                 leadingHalfspaceDualRows=[dict(row=int(active_rows[j]), dual=_evidence_number(dual[j]),
                     necessaryRadiusLowerBoundMeters=_evidence_number(lower[active_rows[j]])) for j in order],
                 elapsedSeconds=time.monotonic() - start)
    if status != 'Solved':
        proof['feasible'] = False
        return None, proof
    x = np.asarray(solution.x)
    if not np.isfinite(x).all():
        proof.update(status='NonfiniteSolution',solverStatus=status,feasible=False)
        return None,proof
    d = np.zeros((len(K), 3))
    d[free] = x[:count].reshape(-1, 3) * scale
    radius = float(np.linalg.norm(motion@d, axis=1).max(initial=0.))
    raw_residual=b-section_A@d[free].ravel()
    adjusted_residual=raw_residual.copy()
    if elastic_count:
        allowances=x[count:]*scale
        soft=np.flatnonzero(groups>=0)
        adjusted_residual[soft]-=np.array([allowances[group_rank[int(groups[i])]]for i in soft])
        proof.update(elasticAllowancesMeters=allowances.tolist(),
                     maximumUnrelaxedResidualMeters=float(raw_residual.max(initial=0.)),
                     originalHalfspacesSatisfied=bool(raw_residual.max(initial=0.)<=residual_m),
                     proposalOnly=True)
    residual = float(adjusted_residual.max(initial=0.))
    requested = float(x[-1] * scale) if minimize_radius else maximum_motion_m
    bound_error = max(0., radius - requested)
    feasible = bool(np.isfinite(d).all() and residual <= residual_m and
                    bound_error <= residual_m)
    proof.update(feasible=feasible, maximumResidualMeters=max(0., residual),
                 maximumMotionBoundErrorMeters=bound_error,
                 maximumActualMotionMeters=radius,maximumParameterNormMeters=float(np.linalg.norm(d,axis=1).max(initial=0.)),protectedStationsLiteral=True)
    if minimize_radius:
        proof.update(minimumRadiusMeters=requested,
                     minimumRadiusDualLowerBoundMeters=float(solution.obj_val_dual * scale))
    return d if feasible else None, proof
