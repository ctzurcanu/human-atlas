"""Bounded local/global rigid-profile graph fitting with explicit fixed nodes.

The caller supplies actual source adjacency and any reviewed interface groups.
No adjacency is inferred from spatial proximity. This numerical graph energy
does not prove emitted tissue validity or anatomical placement.

Local/global rotation fitting follows Sorkine and Alexa, SGP 2007:
https://igl.ethz.ch/projects/ARAP/
"""
import time
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import splu


def fit_profile_network(source, targets, edges, *, fixed=(), edge_weights=None,
                        target_weights=None, reference_rotations=None,
                        rotation_regularization=.05, iterations=20,
                        wall_seconds=10., fixed_at_targets=False):
    p = np.asarray(source, float)
    t = np.asarray(targets, float)
    e = np.asarray(edges, int).reshape(-1, 2)
    n = len(p)
    f = np.unique(np.asarray(fixed, int))
    if p.shape != (n, 3) or t.shape != p.shape or not n:
        raise ValueError('Corresponding source and target profile graph nodes required')
    if not np.isfinite(np.r_[p.ravel(),t.ravel()]).all():
        raise ValueError('Non-finite profile graph')
    if np.any(e < 0) or np.any(e >= n) or np.any(e[:,0] == e[:,1]) or np.any(f < 0) or np.any(f >= n):
        raise ValueError('Invalid source graph or fixed node')
    if not 1 <= iterations <= 100 or not 0 < wall_seconds <= 30:
        raise ValueError('Finite native iteration and wall limits required')
    if not np.isfinite(rotation_regularization) or rotation_regularization < 0:
        raise ValueError('Nonnegative finite rotation regularization required')
    if not isinstance(fixed_at_targets,bool):
        raise ValueError('Explicit boolean fixed-at-target policy required')
    w = np.ones(len(e)) if edge_weights is None else np.asarray(edge_weights,float)
    a = np.ones(n) if target_weights is None else np.asarray(target_weights,float)
    R0 = np.repeat(np.eye(3)[None],n,axis=0) if reference_rotations is None else np.asarray(reference_rotations,float)
    if w.shape != (len(e),) or a.shape != (n,) or np.any(w <= 0) or np.any(a <= 0) or not np.isfinite(np.r_[w,a]).all():
        raise ValueError('Positive finite graph and target weights required')
    if R0.shape != (n,3,3) or not np.isfinite(R0).all() or not np.allclose(R0.transpose(0,2,1)@R0,np.eye(3),atol=1e-11,rtol=0) or np.any(np.linalg.det(R0)<1-1e-11):
        raise ValueError('Reference profile rotations must be proper')
    start = time.monotonic()
    i,j = e.T
    degree = np.bincount(np.r_[i,j],weights=np.r_[w,w],minlength=n)
    adjacency = coo_matrix((np.r_[w,w],(np.r_[i,j],np.r_[j,i])),shape=(n,n)).tocsr()
    H = diags(degree+a)-adjacency
    free = np.setdiff1d(np.arange(n),f)
    factor = splu(H[free][:,free].tocsc()) if len(free) else None
    dx = p[i]-p[j]
    squared = np.einsum('ij,ij->i',dx,dx)*w
    rw = rotation_regularization*np.bincount(np.r_[i,j],weights=np.r_[squared,squared],minlength=n)
    q = t.copy()
    pinned=t[f].copy() if fixed_at_targets else p[f].copy()
    pinned_rotations=R0[f].copy() if fixed_at_targets else np.repeat(np.eye(3)[None],len(f),axis=0)
    q[f] = pinned
    R = R0.copy()
    R[f] = pinned_rotations
    energies = []
    for _ in range(iterations):
        if time.monotonic()-start > wall_seconds:
            raise TimeoutError('Native profile network wall deadline reached')
        dy = q[i]-q[j]
        covariance = rw[:,None,None]*R0
        outer = np.einsum('ni,nj,n->nij',dy,dx,w)
        np.add.at(covariance,i,outer)
        np.add.at(covariance,j,outer)
        u,_,vt = np.linalg.svd(covariance)
        D = np.repeat(np.eye(3)[None],n,axis=0)
        D[:,-1,-1] = np.linalg.det(u@vt)
        R = u@D@vt
        R[f] = pinned_rotations
        rotated = np.einsum('nij,nj->ni',.5*(R[i]+R[j]),dx)*w[:,None]
        rhs = a[:,None]*t
        np.add.at(rhs,i,rotated)
        np.add.at(rhs,j,-rotated)
        if len(free):
            q[free] = factor.solve(rhs[free]-H[free][:,f]@pinned)
        q[f] = pinned
        edge_error = (q[i]-q[j])-np.einsum('nij,nj->ni',.5*(R[i]+R[j]),dx)
        energies.append(float(np.sum(w[:,None]*edge_error**2)+np.sum(a[:,None]*(q-t)**2)))
    return dict(centres=q,rotations=R,energies=energies,
                target_residual_m=np.linalg.norm(q-t,axis=1),
                fixed_nodes_literal=bool(np.array_equal(q[f],pinned)),
                fixed_node_reference='explicit-target' if fixed_at_targets else 'source',
                fixed_rotations_literal=bool(np.array_equal(R[f],pinned_rotations)),
                adjacency_inferred=False,anatomical_acceptance=False,
                elapsed_seconds=time.monotonic()-start)
