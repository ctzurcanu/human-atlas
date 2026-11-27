"""Bounded cutting-plane exchange with full original-row acceptance.

Motion support is only a scheduling key. No inequality is approximated,
discarded, or inferred from another row; every original row is replayed after
every subproblem, and each added exact row remains active thereafter.
"""
import time
import numpy as np
from scipy import sparse
from .hybrid_section_clearance import solve_hybrid_sections


def violated_support_representatives(matrix, residual, active, *, tolerance=1e-10):
    matrix = sparse.csr_matrix(matrix)
    chosen = {}
    for i in np.flatnonzero((residual > tolerance) & ~active):
        columns = matrix.indices[matrix.indptr[i]:matrix.indptr[i + 1]]
        key = tuple(np.unique(columns // 3))
        previous = chosen.get(key)
        if previous is None or residual[i] > residual[previous]:
            chosen[key] = int(i)
    return np.asarray(sorted(chosen.values()), dtype=np.int64)


def solve_sections_generated(stiffness, hard, hard_rhs, soft, soft_rhs, *,
        maximum_rounds=6, stage_seconds=3., wall_seconds=24., checkpoint=None,
        initial_active_ids=None, retain_active_ids=False, subproblem_solver=None, **kwargs):
    if not 1 <= maximum_rounds <= 32 or not 0 < stage_seconds <= 8 or not 0 < wall_seconds <= 30:
        raise ValueError('Explicit finite stage and total budgets required')
    A = sparse.csr_matrix(hard)
    A.eliminate_zeros()
    b = np.asarray(hard_rhs, float)
    tick = time.monotonic()
    active = b >= -1e-10
    if initial_active_ids is not None:
        seed = np.asarray(initial_active_ids)
        if seed.ndim != 1 or seed.dtype.kind not in 'iu' or np.any(seed < 0) or np.any(seed >= len(b)):
            raise ValueError('Exact original-row indices required for saved active constraints')
        active[seed] = True
    records = []
    for number in range(maximum_rounds):
        if time.monotonic() - tick + stage_seconds > wall_seconds:
            reason = 'Total native stage budget reached'
            break
        ids = np.flatnonzero(active)
        solver = solve_hybrid_sections if subproblem_solver is None else subproblem_solver
        d, proof = solver(stiffness, A[ids], b[ids], soft, soft_rhs,
                                        time_limit_seconds=stage_seconds, **kwargs)
        row = dict(round=number, activeRows=len(ids), solver=proof)
        if retain_active_ids:
            row['activeOriginalRowIds'] = ids.tolist()
        if d is None:
            records.append(row)
            if checkpoint is not None:
                checkpoint(number, None, row)
            reason = 'Finite native subproblem did not return a feasible point'
            break
        residual = b - A @ d.ravel()
        bad = residual > 1e-10
        added = violated_support_representatives(A, residual, active)
        row.update(fullOriginalRows=len(b), violatedOriginalRows=int(bad.sum()),
                   fullMaximumResidual_m=float(np.maximum(residual, 0).max(initial=0)),
                   exactRowsScheduled=len(added))
        if retain_active_ids:
            row['nextActiveOriginalRowIds'] = np.union1d(ids, added).tolist()
        records.append(row)
        if checkpoint is not None:
            checkpoint(number, d, row)
        if not bad.any():
            return d, dict(complete=True, feasible=True, allOriginalHalfspacesChecked=True,
                           fullOriginalRows=len(b), maximumHardResidual_m=row['fullMaximumResidual_m'],
                           rounds=records, elapsedSeconds=time.monotonic() - tick,
                           nativeSeconds=sum(r['solver']['nativeSeconds'] for r in records),
                           activeOriginalRowIds=ids.tolist(),
                           anatomicalAcceptance=False)
        if not len(added):
            reason = 'Active rows failed independent full replay'
            break
        active[added] = True
    else:
        reason = 'Declared constraint-generation rounds exhausted'
    return None, dict(complete=True, feasible=False, reason=reason, rounds=records,
                      elapsedSeconds=time.monotonic() - tick)
