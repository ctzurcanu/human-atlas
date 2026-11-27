"""Bounded local/global ARAP for a declared free surface core.

Uniform positive edge weights avoid negative cotangent weights on source slivers.
Attachments are exact Dirichlet controls. This minimizes edge distortion but does
not guarantee injectivity, thickness, calibre, collision freedom or anatomy.
"""
import time
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import factorized


def arap_core(source, edges, boundary, positions, initial, *, maximum_steps=32, wall_seconds=2.,cache_fixed_rotations=True):
    started = time.monotonic()
    p, q = np.asarray(source, float), np.array(initial, float, copy=True)
    e, fixed = np.asarray(edges, int), np.asarray(boundary, int)
    target = np.asarray(positions, float)
    if p.ndim != 2 or p.shape[1] != 3 or q.shape != p.shape or target.shape != (len(fixed), 3):
        raise ValueError('Corresponding source, initial and attachment positions required')
    if not np.isfinite(np.r_[p.ravel(), q.ravel(), target.ravel()]).all() or maximum_steps < 1 or not 0 < wall_seconds <= 10:
        raise ValueError('Finite geometry and bounded positive iteration/time budgets required')
    if e.ndim != 2 or e.shape[1] != 2 or np.any(e < 0) or np.any(e >= len(p)):
        raise ValueError('Explicit source edges required')
    free = np.setdiff1d(np.arange(len(p)), fixed)
    q[fixed] = target
    if not len(free):
        return q, dict(iterations=0, converged=True, elapsed_seconds=time.monotonic()-started,
                       fixed_controls_literal=True, energy=[])
    i, j = e.T
    adj = coo_matrix((np.ones(2*len(e)), (np.r_[i,j],np.r_[j,i])), shape=(len(p),len(p))).tocsr()
    L = diags(np.asarray(adj.sum(1)).ravel()) - adj
    solve = factorized(L[free][:,free].tocsc())
    fixed_rhs = L[free][:,fixed] @ target
    dp = p[i] - p[j]
    # A vertex with no free neighbour has a constant covariance: both ends
    # of every incident edge stay fixed. Its optimal rotation can be reused.
    # Boundary vertices adjacent to free controls remain active, using their
    # complete edge stars, including fixed/fixed incident edges.
    free_mask=np.zeros(len(p),bool);free_mask[free]=True
    touching=free_mask[i]|free_mask[j]
    active=np.unique(e[touching])if cache_fixed_rotations else np.arange(len(p))
    active_lookup=np.full(len(p),-1,int);active_lookup[active]=np.arange(len(active))
    active_i=active_lookup[i]>=0;active_j=active_lookup[j]>=0
    R=None
    energy, converged = [], False
    for step in range(maximum_steps):
        if time.monotonic()-started >= wall_seconds:
            raise TimeoutError('Native ARAP core deadline; incomplete candidate rejected')
        dq = q[i] - q[j]
        outer = np.einsum('ni,nj->nij', dp, dq)
        if R is None:
            covariance = np.zeros((len(p),3,3))
            np.add.at(covariance, i, outer)
            np.add.at(covariance, j, outer)
        else:
            covariance=np.zeros((len(active),3,3))
            np.add.at(covariance,active_lookup[i[active_i]],outer[active_i])
            np.add.at(covariance,active_lookup[j[active_j]],outer[active_j])
        U, _, Vt = np.linalg.svd(covariance)
        V = Vt.transpose(0,2,1)
        sign = np.linalg.det(V @ U.transpose(0,2,1))
        V[:,:,2] *= sign[:,None]
        if R is None:R=V @ U.transpose(0,2,1)
        else:R[active]=V @ U.transpose(0,2,1)
        desired = np.einsum('nij,nj->ni', (R[i]+R[j])*.5, dp)
        rhs = np.zeros_like(p)
        np.add.at(rhs, i, desired)
        np.add.at(rhs, j, -desired)
        previous = q.copy()
        q[free] = solve(rhs[free]-fixed_rhs)
        q[fixed] = target
        if not np.isfinite(q).all():
            raise ValueError('Non-finite ARAP core solve')
        dq = q[i]-q[j]
        score = .5 * float(np.sum((dq-np.einsum('nij,nj->ni', R[i], dp))**2)
                             + np.sum((dq-np.einsum('nij,nj->ni', R[j], dp))**2))
        energy.append(score)
        if np.linalg.norm(q-previous,axis=1).max(initial=0) < 1e-10:
            converged = True
            break
    if time.monotonic()-started >= wall_seconds:
        raise TimeoutError('Native ARAP core deadline; incomplete candidate rejected')
    return q, dict(iterations=step+1, converged=converged, energy=energy,
                   maximum_steps=maximum_steps, native_wall_seconds=wall_seconds,
                   fixed_controls_literal=bool(np.array_equal(q[fixed],target)),
                   variable_rotation_controls=len(active),cached_rotation_controls=len(p)-len(active),
                   elapsed_seconds=time.monotonic()-started,
                   injectivity_not_guaranteed=True)
