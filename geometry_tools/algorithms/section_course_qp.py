"""Independent primal solve for fixed-profile section translation constraints.

Separate feasibility from slow convergence of cyclic dual projections. Each
section's local halfspaces are checked by a three-variable linear program;
the full, coupled quadratic metric then chooses a smooth feasible course.
This solver cannot make incompatible obstacle-side choices feasible.
"""
import numpy as np
from scipy.optimize import linprog, minimize, LinearConstraint


def project_section_halfspaces(metric, stations, normals, offsets, *,
                               residual_m=1e-10, maximum_iterations=300):
    K = np.asarray(metric, float)
    s = np.asarray(stations, int)
    n = np.asarray(normals, float)
    b = np.asarray(offsets, float)
    if not np.allclose(np.linalg.norm(n, axis=1), 1., atol=1e-12, rtol=0.):
        raise ValueError('Unit constraint normals required')
    d = np.zeros((len(K), 3))
    if not len(s):
        return d, dict(iterations=0, maximumResidualMeters=0.)
    free = np.flatnonzero(np.diag(K) > 0.)
    scale = .001  # Metres to millimetres; no geometry tolerance is changed.
    feasible = np.zeros_like(d)
    for j in np.unique(s):
        take = s == j
        if K[j, j] <= 0.:
            if np.max(b[take]) > residual_m:
                raise ValueError('A protected source station conflicts with the obstacle constraint')
            continue
        lp = linprog(np.zeros(3), A_ub=-n[take], b_ub=-b[take] / scale,
            bounds=[(None, None)] * 3, method='highs',
            options={'primal_feasibility_tolerance': 1e-10,
                     'dual_feasibility_tolerance': 1e-10})
        if not lp.success:
            raise ValueError(f'Incompatible local obstacle halfspaces at source station {int(j)} (LP status {lp.status})')
        feasible[j] = lp.x
    # Ordering is section-major XYZ. Fixed rows never enter the optimization.
    mapping = {int(j): i for i, j in enumerate(free)}
    A = np.zeros((len(s), 3 * len(free)))
    for i, j in enumerate(s):
        if int(j) in mapping:
            k = mapping[int(j)] * 3
            A[i, k:k + 3] = n[i]
    # Whiten the positive metric. Direct SLSQP in section coordinates spends
    # its iteration budget learning the fourth-difference energy curvature.
    # Here displacement/scale = L @ x and the objective is exactly |x|²/2.
    L = np.linalg.cholesky(K[np.ix_(free, free)])
    B = A @ np.kron(L, np.eye(3))
    def objective(x):
        return float(.5 * (x @ x))
    def gradient(x):
        return x
    initial = np.linalg.solve(L, feasible[free]).ravel()
    result = minimize(objective, initial, jac=gradient,
        constraints=[LinearConstraint(B, b / scale, np.inf)], method='SLSQP',
        options={'maxiter': maximum_iterations, 'ftol': 1e-13})
    d[free] = (L @ result.x.reshape(-1, 3)) * scale
    residual = float(np.max(b - np.einsum('ij,ij->i', n, d[s])))
    if not np.isfinite(d).all() or residual > residual_m:
        raise ValueError(f'Primal course solve did not pass: {result.message}; residual {residual:g} m')
    return d, dict(method='independent local LP feasibility, whitened coupled primal SLSQP',
        iterations=int(result.nit), maximumResidualMeters=max(0., residual),
        minimumConverged=bool(result.success), solverMessage=str(result.message),
        objective=float(result.fun), protectedStationsLiteral=True,
        residualBoundMeters=residual_m)


def project_weighted_section_halfspaces(metric, weights, normals, offsets, *,
                                        residual_m=1e-10, maximum_iterations=300):
    """Project barycentric contact constraints shared by several sections.

    Unlike assigning a whole facet plane to every incident vertex, a row
    constrains only the actual material contact point. Negative weights also
    express separation between two moving source features.
    """
    K = np.asarray(metric, float); w = np.asarray(weights, float)
    n = np.asarray(normals, float); b = np.asarray(offsets, float)
    if w.shape != (len(b), len(K)) or n.shape != (len(b), 3):
        raise ValueError('One section-weight row and unit normal per constraint required')
    if not np.allclose(np.linalg.norm(n, axis=1), 1., atol=1e-12, rtol=0.):
        raise ValueError('Unit constraint normals required')
    d = np.zeros((len(K), 3)); free = np.flatnonzero(np.diag(K) > 0.)
    scale = .001
    A = (w[:, free, None] * n[:, None]).reshape(len(b), -1)
    protected_rows = np.linalg.norm(A, axis=1) == 0.
    if np.any(b[protected_rows] > residual_m):
        raise ValueError('A protected source contact conflicts with the obstacle constraint')
    if not len(b):
        return d, dict(iterations=0, maximumResidualMeters=0.)
    keep = ~protected_rows; A = A[keep]; rhs = b[keep] / scale
    lp = linprog(np.zeros(A.shape[1]), A_ub=-A, b_ub=-rhs,
        bounds=[(None, None)] * A.shape[1], method='highs',
        options={'primal_feasibility_tolerance': 1e-10,
                 'dual_feasibility_tolerance': 1e-10})
    if not lp.success:
        raise ValueError(f'Incompatible weighted obstacle halfspaces (LP status {lp.status})')
    L = np.linalg.cholesky(K[np.ix_(free, free)])
    B = A @ np.kron(L, np.eye(3))
    initial = np.linalg.solve(L, lp.x.reshape(-1, 3)).ravel()
    result = minimize(lambda x: .5 * (x @ x), initial, jac=lambda x: x,
        constraints=[LinearConstraint(B, rhs, np.inf)], method='SLSQP',
        options={'maxiter': maximum_iterations, 'ftol': 1e-13})
    d[free] = (L @ result.x.reshape(-1, 3)) * scale
    residual = float(np.max(b - np.einsum('ij,ij->i', n, w @ d)))
    if not np.isfinite(d).all() or residual > residual_m:
        raise ValueError(f'Weighted course residual exceeds the declared bound: {residual:g} m')
    return d, dict(method='barycentric contact, whitened primal QP',
        iterations=int(result.nit), minimumConverged=bool(result.success),
        solverMessage=str(result.message), maximumResidualMeters=max(0., residual),
        residualBoundMeters=residual_m, protectedStationsLiteral=True)
