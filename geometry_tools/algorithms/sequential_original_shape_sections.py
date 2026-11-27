"""Finite local source-shape proposals, accepted only by complete exact replay."""
import time
import numpy as np
from scipy import sparse
from .support_exchange_sections import solve_sections_generated


def solve_sequential_original_shape(stiffness, hard, hard_rhs, soft, soft_rhs, *,
        exact_source_shape_rows, maximum_stages=8, wall_seconds=28., checkpoint=None,
        state_checkpoint=None, initial_local_rows=None, initial_active_ids=None, additional_exact_rows=None, **kwargs):
    if not 1 <= maximum_stages <= 12 or not 0 < wall_seconds <= 30:
        raise ValueError('Explicit finite sequential-shape stage and wall limits required')
    A, b = sparse.csr_matrix(hard), np.asarray(hard_rhs, float)
    local_A, local_b = sparse.csr_matrix((0, A.shape[1])), np.array([])
    if initial_local_rows is not None:
        local_A, local_b = initial_local_rows
        local_A, local_b = sparse.csr_matrix(local_A), np.asarray(local_b, float)
        if local_A.shape != (len(local_b), A.shape[1]):
            raise ValueError('Saved local rows must match the unchanged original pose variables')
        A, b = sparse.vstack([A, local_A], format='csr'), np.r_[b, local_b]
    seed = initial_active_ids
    records = []
    tick = time.monotonic()
    for stage in range(maximum_stages):
        remaining = wall_seconds-(time.monotonic()-tick)
        if remaining < 2.:
            reason = 'Declared sequential source-shape time budget reached'
            break
        d, proof = solve_sections_generated(stiffness, A, b, soft, soft_rhs,
            maximum_rounds=24, stage_seconds=2., wall_seconds=remaining,
            initial_active_ids=seed,
            checkpoint=(lambda number, value, row: checkpoint(stage, number, value, row))
                       if checkpoint is not None else None, **kwargs)
        record = dict(stage=stage, originalAndLocalRows=len(b), solver=proof)
        if d is None:
            records.append(record)
            reason = 'Finite local stage did not return a fully halfspace-feasible point'
            break
        next_A, next_b, replay = exact_source_shape_rows(d)
        next_A = sparse.csr_matrix(next_A)
        next_b = np.asarray(next_b, float)
        if next_A.shape != (len(next_b), A.shape[1]) or not replay.get('allOriginalFacetsChecked'):
            raise ValueError('Complete original-facet replay and corresponding local rows required')
        extra_bad=0
        if additional_exact_rows is not None:
            extra_A,extra_b,extra_replay=additional_exact_rows(d)
            extra_A=sparse.csr_matrix(extra_A);extra_b=np.asarray(extra_b,float)
            if extra_A.shape!=(len(extra_b),A.shape[1]) or not extra_replay.get('complete'):
                raise ValueError('Complete additional exact replay and matching local rows required')
            extra_bad=int(extra_replay['violatingExactConstraints'])
            next_A=sparse.vstack([next_A,extra_A],format='csr');next_b=np.r_[next_b,extra_b]
            record['additionalExactConstraintReplay']=extra_replay
        record.update(exactSourceShapeReplay=replay, localSourceShapeRowsAdded=len(next_b))
        records.append(record)
        if replay['violatingOriginalFacets'] == 0 and extra_bad==0:
            return d, dict(complete=True, feasible=True, fullOriginalHalfspacesChecked=True,
                           fullOriginalSourceShapeChecked=True, stages=records,
                           elapsedSeconds=time.monotonic()-tick, anatomicalAcceptance=False)
        if not len(next_b):
            reason = 'Source-shape violation has no resolved local pose correction'
            break
        seed = np.r_[np.asarray(proof['activeOriginalRowIds'], dtype=np.int64),
                     np.arange(len(b), len(b)+len(next_b))]
        A, b = sparse.vstack([A, next_A], format='csr'), np.r_[b, next_b]
        local_A, local_b = sparse.vstack([local_A, next_A], format='csr'), np.r_[local_b, next_b]
        if state_checkpoint is not None:
            state_checkpoint(stage, d, local_A, local_b, seed, record)
    else:
        reason = 'Declared sequential source-shape stages exhausted'
    return None, dict(complete=True, feasible=False, reason=reason, stages=records,
                      elapsedSeconds=time.monotonic()-tick)
