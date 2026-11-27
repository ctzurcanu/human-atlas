"""Eliminate declared fixed three-vector variables without penalty or soft locks."""
import numpy as np
from scipy import sparse

def solve_with_fixed_sections(solver, stiffness, hard, hard_rhs, soft, soft_rhs, *, fixed, **kwargs):
    fixed=np.asarray(fixed,bool)
    S=sparse.csc_matrix(stiffness)
    if fixed.shape!=(S.shape[1],) or np.all(fixed):
        raise ValueError('A declared mask and at least one moving section are required')
    free=np.flatnonzero(~fixed)
    columns=(free[:,None]*3+np.arange(3)).ravel()
    options=dict(kwargs)
    if options.get('distortion') is not None:
        options['distortion']=sparse.csc_matrix(options['distortion'])[:,columns]
    value,proof=solver(S[:,free],sparse.csc_matrix(hard)[:,columns],hard_rhs,
        sparse.csc_matrix(soft)[:,columns],soft_rhs,**options)
    proof=dict(proof,declaredFixedSectionVectors=int(fixed.sum()),
        fixedSectionMethod='Exact parameter elimination; no penalty or approximate lock')
    if value is None:return None,proof
    expanded=np.zeros((len(fixed),3))
    expanded[free]=value
    assert not expanded[fixed].any()
    return expanded,proof
