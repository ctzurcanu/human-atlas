"""Generate source-shape cones only where full exact replay requires them.

The caller's check must replay exact section poses against every original
source facet, not the linearized source operator. No source mesh or physical
budget is reset. A point is returned only after the full original halfspace
system AND the complete independent source-shape check pass.
"""
import time
import numpy as np
from scipy import sparse
from .support_exchange_sections import solve_sections_generated


def solve_lazy_source_shape(stiffness, hard, hard_rhs, soft, soft_rhs, *,
        all_distortion, all_offsets, exact_shape_violations, maximum_stages=4,
        wall_seconds=28., checkpoint=None, **kwargs):
    if not 1 <= maximum_stages <= 8 or not 0 < wall_seconds <= 30:
        raise ValueError('Explicit finite source-shape stage and wall limits required')
    C = sparse.csr_matrix(all_distortion)
    offset = np.asarray(all_offsets).reshape(-1, 6)
    if C.shape[0] != 6 * len(offset):
        raise ValueError('One six-row source operator per original moving facet required')
    selected = np.zeros(len(offset), bool)
    seed = None
    records = []
    tick = time.monotonic()
    for stage in range(maximum_stages):
        remaining = wall_seconds - (time.monotonic() - tick)
        if remaining < 2.:
            reason = 'Total declared shape-generation time budget reached'
            break
        ids = np.flatnonzero(selected)
        rows = (ids[:, None] * 6 + np.arange(6)).ravel()
        d, proof = solve_sections_generated(stiffness, hard, hard_rhs, soft, soft_rhs,
            maximum_rounds=24, stage_seconds=2., wall_seconds=remaining,
            initial_active_ids=seed,
            distortion=C[rows] if len(ids) else None,
            distortion_offset=offset[ids] if len(ids) else None,
            checkpoint=(lambda number, value, report: checkpoint(stage, number, value, report))
                       if checkpoint is not None else None, **kwargs)
        record = dict(stage=stage, activeOriginalFacetCones=len(ids), solver=proof)
        if d is None:
            records.append(record)
            reason = 'Finite constrained stage did not return a fully halfspace-feasible point'
            break
        bad, replay = exact_shape_violations(d)
        bad = np.asarray(bad, dtype=np.int64)
        if bad.ndim != 1 or np.any(bad < 0) or np.any(bad >= len(offset)):
            raise ValueError('Exact source replay must return original facet indices')
        record.update(exactOriginalShapeReplay=replay, violatingOriginalFacets=len(bad))
        records.append(record)
        if not len(bad):
            return d, dict(complete=True, feasible=True, fullOriginalHalfspacesChecked=True,
                           fullOriginalSourceShapeChecked=True, stages=records,
                           elapsedSeconds=time.monotonic() - tick, anatomicalAcceptance=False)
        if np.all(selected[bad]):
            reason = 'Exact pose still violates a constrained source facet; relinearization required'
            break
        selected[bad] = True
        seed = np.asarray(proof['activeOriginalRowIds'], dtype=np.int64)
    else:
        reason = 'Declared shape-generation stages exhausted'
    return None, dict(complete=True, feasible=False, reason=reason, stages=records,
                      elapsedSeconds=time.monotonic() - tick)
