"""Finite constraint generation; every original halfspace is replayed.

Subproblems contain a growing exact subset. A proposal is returned only when
it satisfies the complete original matrix as well as the subproblem's motion
and immutable-source metric cones. This reduces factorization work, never the
final acceptance scope or tolerances.
"""
import time
import numpy as np
from scipy import sparse
from .hybrid_section_clearance import solve_hybrid_sections as solve_soft_sections


def solve_sections_generated(stiffness,hard,hard_rhs,soft,soft_rhs,*,
        maximum_rounds=6,stage_seconds=2.,wall_seconds=18.,checkpoint=None,**kwargs):
    if not 1<=maximum_rounds<=8 or not 0<stage_seconds<=8 or not 0<wall_seconds<=30:
        raise ValueError('Explicit finite stage and total budgets required')
    A=sparse.csr_matrix(hard);b=np.asarray(hard_rhs,float);tick=time.monotonic()
    active=b>=-1e-10;records=[]
    for round_number in range(maximum_rounds):
        if time.monotonic()-tick+stage_seconds>wall_seconds:
            return None,dict(complete=True,feasible=False,reason='Total native stage budget reached',rounds=records)
        ids=np.flatnonzero(active)
        d,proof=solve_soft_sections(stiffness,A[ids],b[ids],soft,soft_rhs,time_limit_seconds=stage_seconds,**kwargs)
        row=dict(round=round_number,activeRows=len(ids),solver=proof)
        if d is None:
            records.append(row)
            if checkpoint is not None:checkpoint(round_number,None,row)
            return None,dict(complete=True,feasible=False,reason='Finite native subproblem did not return a feasible point',rounds=records,elapsedSeconds=time.monotonic()-tick)
        residual=b-A@d.ravel();bad=residual>1e-10
        row.update(fullOriginalRows=len(b),fullMaximumResidual_m=float(np.maximum(residual,0).max(initial=0)),violatedOriginalRows=int(bad.sum()))
        records.append(row)
        if checkpoint is not None:checkpoint(round_number,d,row)
        if not bad.any():
            return d,dict(complete=True,feasible=True,allOriginalHalfspacesChecked=True,fullOriginalRows=len(b),
                maximumHardResidual_m=row['fullMaximumResidual_m'],rounds=records,
                nativeSeconds=sum(x['solver']['nativeSeconds'] for x in records),elapsedSeconds=time.monotonic()-tick,
                optimalityCertified=proof.get('optimalityCertified',False),anatomicalAcceptance=False)
        active|=bad
    return None,dict(complete=True,feasible=False,reason='Declared constraint-generation rounds exhausted',rounds=records,elapsedSeconds=time.monotonic()-tick)
