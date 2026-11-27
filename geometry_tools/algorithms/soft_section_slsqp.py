"""Bounded small-compartment convex solve with independently replayed limits."""
import time
import numpy as np
from scipy import sparse
from scipy.optimize import minimize


def solve_small_compartment(stiffness, A, b, B, target, distortion, *, maximum_motion_m,
                            distortion_bound, soft_weight=10000., wall_seconds=5.):
    started = time.monotonic()
    if not 0 < wall_seconds <= 10:
        raise ValueError('Small-compartment solve needs a finite deadline up to ten seconds')
    A, B, C = [sparse.csr_matrix(m) for m in [A, B, distortion]]
    n = A.shape[1]
    if n % 3 or n > 600 or B.shape[1] != n or C.shape[1] != n or C.shape[0] % 6:
        raise ValueError('Small explicit three-coordinate stations and six-row facet fields required')
    dense_A, dense_B = A.toarray(), B.toarray()
    bmm, tmm = np.asarray(b)*1000., np.asarray(target)*1000.
    blocks = (C.toarray()*.001).reshape(-1, 6, n)
    gram = sparse.kron(sparse.csr_matrix(stiffness).T@sparse.csr_matrix(stiffness), sparse.eye(3)).toarray()
    gram += np.eye(n)*1e-6
    radius2 = (maximum_motion_m*1000.)**2

    def guard():
        if time.monotonic()-started > wall_seconds:
            raise TimeoutError('Finite small-compartment SLSQP deadline')

    def objective(x):
        guard()
        slack = np.maximum(tmm-dense_B@x, 0.)
        return .5*x@gram@x+.5*soft_weight*slack@slack

    def gradient(x):
        guard()
        return gram@x-soft_weight*dense_B.T@np.maximum(tmm-dense_B@x, 0.)

    def radius(x):
        guard()
        return radius2-np.sum(x.reshape(-1, 3)**2, axis=1)

    def radius_jac(x):
        result = np.zeros((n//3, n))
        for i in range(n//3):
            result[i, i*3:i*3+3] = -2*x[i*3:i*3+3]
        return result

    def facets(x):
        guard()
        return distortion_bound**2-np.sum((blocks@x)**2, axis=1)

    def facets_jac(x):
        return -2*np.einsum('fi,fij->fj', blocks@x, blocks)

    try:
        answer = minimize(objective, np.zeros(n), jac=gradient, method='SLSQP',
                          constraints=[dict(type='ineq', fun=lambda x:dense_A@x-bmm, jac=lambda x:dense_A),
                                       dict(type='ineq', fun=radius, jac=radius_jac),
                                       dict(type='ineq', fun=facets, jac=facets_jac)],
                          options=dict(maxiter=100, ftol=1e-12))
    except TimeoutError:
        return None, dict(complete=True, feasible=False, timedOut=True,
                          elapsedSeconds=time.monotonic()-started, wallSeconds=wall_seconds)
    d = answer.x.reshape(-1, 3)/1000.
    hard = float(np.maximum(np.asarray(b)-A@d.ravel(), 0).max(initial=0.))
    tangent = float(np.linalg.norm((C@d.ravel()).reshape(-1, 6), axis=1).max(initial=0.))
    motion = float(np.linalg.norm(d, axis=1).max(initial=0.))
    feasible = bool(answer.success and hard <= 1e-10 and tangent <= distortion_bound+1e-9 and motion <= maximum_motion_m+1e-10)
    proof = dict(complete=True, feasible=feasible, timedOut=False, method='Small source-compartment SLSQP; identical hard intervals and facet/motion cones',
                 nativeSuccess=bool(answer.success), nativeStatus=int(answer.status), nativeMessage=str(answer.message), iterations=int(answer.nit),
                 maximumHardResidual_m=hard, maximumFacetDistortion=tangent, maximumMotion_m=motion,
                 maximumFacetDistortionBound=distortion_bound, maximumMotionBound_m=maximum_motion_m,
                 elapsedSeconds=time.monotonic()-started, wallSeconds=wall_seconds,
                 maximumRemainingSoftSlack_m=float(np.maximum(np.asarray(target)-B@d.ravel(), 0).max(initial=0.)),
                 anatomicalAcceptance=False)
    return (d if feasible else None), proof
