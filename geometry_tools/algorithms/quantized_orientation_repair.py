"""Finite integer-word orientation repair; source boundaries stay fixed.

This fixes encoding, not anatomy. Full emitted topology and intersection
audits are separate, mandatory obligations before using the result.
"""
import time
import numpy as np
from scipy.optimize import milp, Bounds, LinearConstraint
from scipy.sparse import lil_matrix, hstack, eye, vstack


def repair_orientation_words(positions, faces, protected, *, maximum_m=2e-8,
                             seconds=25., iterations=5, corner_normals=None):
    start = time.monotonic()
    p, f = np.asarray(positions, float), np.asarray(faces, int)
    q = p.astype('<f4').astype(float)
    anchor = q.copy()
    fixed = np.asarray(protected, bool)
    tri = p[f]
    reference = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    reference /= np.linalg.norm(reference, axis=1)[:, None]
    reference = reference[:, None]
    if corner_normals is not None:
        field = np.asarray(corner_normals, float)
        if field.shape != (len(f), 3, 3) or np.any(np.linalg.norm(field, axis=2) == 0):
            raise ValueError('Three nonzero source corner normals are required per face')
        field = field/np.linalg.norm(field, axis=2)[:, :, None]
        reference = np.concatenate([reference, field], axis=1)
    def area(at):
        t = at[f]
        return np.einsum('fki,fi->fk', reference,
                        np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])).min(1)
    records = []
    for iteration in range(iterations):
        bad = np.flatnonzero(area(q) <= 0)
        if not len(bad):
            break
        if time.monotonic()-start >= seconds:
            raise TimeoutError('Quantized orientation repair deadline')
        seeds = np.unique(f[bad])
        near = np.flatnonzero(np.any(np.isin(f, seeds), axis=1))
        moving = np.unique(f[near])
        moving = moving[~fixed[moving]]
        if not len(moving) or len(moving) > 256:
            raise ValueError('No bounded free internal patch for integer-word repair')
        incident = np.flatnonzero(np.any(np.isin(f, moving), axis=1))
        t = q[f[incident]]
        n = reference[incident]
        length = np.linalg.norm(t-np.roll(t, 1, axis=1), axis=2).max(1)
        current = np.einsum('fki,fi->fk', n,
                           np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]))/length[:, None]/1e-9
        derivative = np.stack([np.cross((t[:, 1]-t[:, 2])[:, None], n),
                               np.cross((t[:, 2]-t[:, 0])[:, None], n),
                               np.cross(n, (t[:, 1]-t[:, 0])[:, None])], axis=2)/length[:, None, None, None]
        spacing = abs(np.spacing(q[moving].astype('<f4')).astype(float))
        # Avoid an unbounded word count at near-zero coordinates.
        if spacing.min() < 1e-14:
            raise ValueError('Near-zero coordinate needs an explicit alternate frame')
        columns = {int(v): i for i, v in enumerate(moving)}
        width = len(moving)*3
        conditions = n.shape[1]
        a = lil_matrix((len(incident)*conditions, width))
        for row, face in enumerate(f[incident]):
            for corner, vertex in enumerate(face):
                if int(vertex) in columns:
                    column = columns[int(vertex)]
                    a[row*conditions:(row+1)*conditions, column*3:column*3+3] = derivative[row, :, corner]*spacing[column]/1e-9
        a = a.tocsr()
        limits = np.floor((maximum_m/np.sqrt(3))/spacing).ravel()
        remaining_lo = np.ceil((anchor[moving]-maximum_m/np.sqrt(3)-q[moving])/spacing).ravel()
        remaining_hi = np.floor((anchor[moving]+maximum_m/np.sqrt(3)-q[moving])/spacing).ravel()
        lo, hi = np.maximum(-limits, remaining_lo), np.minimum(limits, remaining_hi)
        identity = eye(width, format='csr')
        matrix = vstack([hstack([a, lil_matrix((len(incident)*conditions, width))]),
                         hstack([identity, -identity]), hstack([-identity, -identity])]).tocsr()
        lower = np.r_[.000001-current.ravel(), np.full(width*2, -np.inf)]
        upper = np.r_[np.full(len(incident)*conditions, np.inf), np.zeros(width*2)]
        result = milp(np.r_[np.zeros(width), spacing.ravel()/1e-9],
                      integrality=np.r_[np.ones(width), np.zeros(width)],
                      bounds=Bounds(np.r_[lo, np.zeros(width)], np.r_[hi, limits]),
                      constraints=LinearConstraint(matrix, lower, upper),
                      options={'time_limit': min(5., seconds-(time.monotonic()-start)),
                               'node_limit': 1000, 'mip_rel_gap': .01})
        if result.x is None:
            raise ValueError('Integer-word orientation constraints are infeasible in the bounded patch')
        words = np.rint(result.x[:width]).reshape(-1, 3)
        trial = q.copy()
        trial[moving] = (q[moving]+spacing*words).astype('<f4').astype(float)
        if np.linalg.norm(trial-anchor, axis=1).max() > maximum_m:
            raise ValueError('Integer-word correction exceeds its immutable bound')
        records.append(dict(iteration=iteration, invalid_before=len(bad),
                            invalid_after=int(np.sum(area(trial) <= 0)),
                            moving_controls=moving.tolist(), solver_status=int(result.status)))
        q = trial
    if np.any(area(q) <= 0):
        raise ValueError('Finite integer-word iterations leave invalid facets')
    if not np.array_equal(q[fixed], anchor[fixed]):
        raise ValueError('A protected boundary control moved')
    return q, dict(passed=True, method='bounded integer Float32 words and actual nonlinear orientation replay',
        iterations=records, protected_controls_literal=True,
        maximum_word_correction_m=float(np.linalg.norm(q-anchor, axis=1).max()),
        maximum_total_encoding_error_m=float(np.linalg.norm(q-p, axis=1).max()),
        full_emitted_audits_required=True, elapsed_seconds=time.monotonic()-start)
