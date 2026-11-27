"""Small-step common landmark flow, with explicit derivative/time bounds.

The same smooth map acts on every member of an assembly. Positive sampled
Jacobians are a numerical guard, not proof of global injectivity or anatomy.
"""
import time

import numpy as np
from scipy.linalg import solve


def landmark_flow(points, landmarks, targets, *, steps=32, wall_seconds=10.,
                  fixed_above=1.08, full_below=.96):
    started = time.monotonic()
    q, c, target = [np.array(x, dtype=float, copy=True) for x in (points, landmarks, targets)]
    if q.ndim != 2 or q.shape[1] != 3 or c.ndim != 2 or c.shape[1] != 3 or target.shape != c.shape:
        raise ValueError('Finite corresponding 3-D points and landmarks required')
    if not np.isfinite(np.r_[q.ravel(), c.ravel(), target.ravel()]).all():
        raise ValueError('Finite geometry required')
    if not 1 <= steps <= 64 or not 0 < wall_seconds <= 15 or not full_below < fixed_above:
        raise ValueError('Explicit finite flow limits required')
    if len(np.unique(c, axis=0)) != len(c):
        raise ValueError('Repeated landmarks are not independent constraints')
    J = np.repeat(np.eye(3)[None], len(q), axis=0)
    minimum = 1.; records = []

    def fade(x):
        t = np.clip((fixed_above-x[:, 1])/(fixed_above-full_below), 0, 1)
        w = t**3*(10-15*t+6*t*t)
        dy = -30*t*t*(1-t)**2/(fixed_above-full_below)
        return w, dy

    for step in range(steps):
        if time.monotonic()-started >= wall_seconds:
            raise TimeoutError('Native landmark flow deadline; no candidate accepted')
        centre = c.mean(0); scale = max(np.ptp(c, axis=0).max(), .01)
        x = (c-centre)/scale
        d = x[:, None]-x[None]
        radius = np.linalg.norm(d, axis=2)
        K = radius**2*np.log(np.maximum(radius, 1e-30))
        P = np.c_[np.ones(len(x)), x]
        A = np.block([[K+np.eye(len(x))*1e-10, P], [P.T, np.zeros((4, 4))]])
        velocity = (target-c)/(steps-step)
        coef = solve(A, np.r_[velocity, np.zeros((4, 3))], assume_a='sym')

        def evaluate(p):
            u = (p-centre)/scale; delta = u[:, None]-x[None]
            r = np.linalg.norm(delta, axis=2)
            kernel = r*r*np.log(np.maximum(r, 1e-30))
            v = kernel@coef[:len(x)]+np.c_[np.ones(len(u)), u]@coef[len(x):]
            gradient = delta*(2*np.log(np.maximum(r, 1e-30))+1)[:, :, None]/scale
            D = np.einsum('nki,kj->nji', gradient, coef[:len(x)])+coef[len(x)+1:].T
            w, dy = fade(p)
            D *= w[:, None, None]; D[:, :, 1] += v*dy[:, None]
            return v*w[:, None], D

        v, D = evaluate(q)
        # Bound every sampled incremental deformation; do not silently change
        # the requested number of steps or run an unbounded retry sequence.
        norm = np.linalg.svd(D, compute_uv=False)[:, 0].max(initial=0)
        if norm > .30:
            raise ValueError('Landmark increment exceeds the declared derivative bound')
        inc = D+np.eye(3)
        determinants = np.linalg.det(inc)
        if determinants.min(initial=1) <= 0:
            raise ValueError('Landmark increment inverts a sampled control')
        J = inc@J; q += v
        cv, _ = evaluate(c); c += cv
        minimum = min(minimum, float(determinants.min(initial=1)))
        records.append(dict(step=step, derivativeNorm=float(norm),
                            remainingLandmarkError_m=float(np.linalg.norm(c-target, axis=1).max())))
    if time.monotonic()-started >= wall_seconds:
        raise TimeoutError('Native landmark flow deadline; no candidate accepted')
    return q, J, dict(steps=records, minimumSampledIncrementJacobian=minimum,
                      landmarkError_m=float(np.linalg.norm(c-target, axis=1).max()),
                      nativeWallSeconds=wall_seconds, elapsedSeconds=time.monotonic()-started,
                      globalInjectivityNotEstablished=True)
