"""Bounded, coupled translation proposals for rigid vessel sections.

Surface witnesses are linearized local constraints, never containment proof.
Every unsatisfied inequality is returned explicitly. No profile is resized.
"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix, hstack, vstack, eye


def clearance_translation_step(profile_count, edges, constraint_profiles,
                               normals, deficits_m, *, fixed_profiles=(),
                               maximum_step_m=.003, wall_seconds=10.,
                               clearance_weight=1000., smoothness_weight=2.):
    """Minimize L1 movement/differential movement and explicit violation slacks.

    Inequalities have the convention ``normal dot delta >= deficit_m``.
    Variables and reported slacks use metres; the numerical LP uses mm.
    Fixed sections remain exactly fixed. Positive slacks forbid acceptance.
    """
    n = int(profile_count)
    edges = np.asarray(edges, dtype=int).reshape(-1, 2)
    owners = np.asarray(constraint_profiles, dtype=int)
    normals = np.asarray(normals, dtype=float)
    deficits = np.asarray(deficits_m, dtype=float)
    fixed = np.unique(np.asarray(fixed_profiles, dtype=int))
    if n < 1 or normals.shape != (len(owners), 3) or deficits.shape != owners.shape:
        raise ValueError('One finite surface inequality per owning profile required')
    if (not np.isfinite(normals).all() or not np.isfinite(deficits).all()
            or np.any(owners < 0) or np.any(owners >= n)
            or np.any(edges < 0) or np.any(edges >= n)
            or np.any(fixed < 0) or np.any(fixed >= n)):
        raise ValueError('Invalid profile ownership or surface witnesses')
    if (not 0 < wall_seconds <= 30 or not 0 < maximum_step_m <= .01
            or not np.isfinite([wall_seconds, maximum_step_m,
                                clearance_weight, smoothness_weight]).all()
            or clearance_weight <= 0 or smoothness_weight < 0):
        raise ValueError('A finite wall deadline and explicit local trust region required')
    lengths = np.linalg.norm(normals, axis=1)
    if not np.allclose(lengths, 1., atol=1e-10, rtol=0):
        raise ValueError('Surface normals must be unit length')
    d = 3*n
    m = len(owners)
    e = 3*len(edges)
    # x, |x|, |edge difference|, surface violation slack
    cols = np.repeat(owners*3, 3) + np.tile(np.arange(3), m)
    H = coo_matrix((normals.ravel(), (np.repeat(np.arange(m), 3), cols)),
                   shape=(m, d)).tocsr()
    row = np.repeat(np.arange(e), 2)
    edge_cols = (edges[:, None, :]*3 + np.arange(3)[None, :, None]).reshape(-1)
    D = coo_matrix((np.tile([1., -1.], e), (row, edge_cols)),
                   shape=(e, d)).tocsr()
    Z = lambda r, c: coo_matrix((r, c)).tocsr()
    I = eye(d, format='csr')
    E = eye(e, format='csr')
    M = eye(m, format='csr')
    A = vstack([
        hstack([I, -I, Z(d, e), Z(d, m)]),
        hstack([-I, -I, Z(d, e), Z(d, m)]),
        hstack([D, Z(e, d), -E, Z(e, m)]),
        hstack([-D, Z(e, d), -E, Z(e, m)]),
        hstack([-H, Z(m, d), Z(m, e), -M]),
    ], format='csr')
    b = np.r_[np.zeros(2*d+2*e), -deficits*1000.]
    bounds = [(-maximum_step_m*1000., maximum_step_m*1000.)] * d
    for i in fixed:
        bounds[3*i:3*i+3] = [(0., 0.)]*3
    bounds += [(0., None)]*(d+e+m)
    objective = np.r_[np.zeros(d), np.ones(d),
                      np.full(e, smoothness_weight), np.full(m, clearance_weight)]
    result = linprog(objective, A_ub=A, b_ub=b, bounds=bounds,
                     method='highs', options={'time_limit': float(wall_seconds)})
    if not result.success:
        # A timed-out solver iterate is never emitted as candidate geometry.
        raise RuntimeError(f'Clearance proposal failed ({result.status}): {result.message}')
    delta = result.x[:d].reshape(n, 3)/1000.
    residual = np.maximum(0., deficits-np.einsum('ij,ij->i', normals, delta[owners]))
    if len(fixed) and np.any(delta[fixed] != 0):
        raise RuntimeError('Solver changed a protected profile')
    return dict(translations_m=delta, inequality_residual_m=residual,
                explicit_solver_slack_m=result.x[-m:]/1000. if m else np.empty(0),
                iterations=int(result.nit), maximum_axis_step_m=maximum_step_m,
                maximum_residual_m=float(residual.max(initial=0.)),
                surface_containment_proved=False, profile_radius_changed=False)
