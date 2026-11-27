"""Remove only byte-identical sparse left-hand sides, keeping strongest RHS.

No rounding, direction clustering or weakened geometric tolerance is used.
The returned source-row map permits independent replay of every original row.
"""
import numpy as np
from scipy import sparse


def exact_halfspace_dominance(matrix,rhs):
    A=sparse.csr_matrix(matrix,copy=True);A.sum_duplicates();A.eliminate_zeros();A.sort_indices()
    b=np.asarray(rhs,float)
    if b.shape!=(A.shape[0],) or not np.isfinite(b).all() or not np.isfinite(A.data).all():
        raise ValueError('Finite sparse halfspaces required')
    # Normalize the sign bit of numerical zero only; zero coefficients were
    # eliminated and every retained coefficient stays literally unchanged.
    representatives=[];selected=[];mapping=np.empty(len(b),int);keys={}
    for i in range(len(b)):
        lo,hi=A.indptr[i:i+2];key=(A.indices[lo:hi].tobytes(),A.data[lo:hi].tobytes())
        j=keys.get(key)
        if j is None:
            j=len(selected);keys[key]=j;selected.append(i);representatives.append(i)
        elif b[i]>b[selected[j]]:selected[j]=i
        mapping[i]=j
    selected=np.asarray(selected,int);reduced=A[selected];bounds=b[selected]
    assert np.all(bounds[mapping]>=b)
    return reduced,bounds,mapping,dict(complete=True,originalRows=len(b),retainedRows=len(bounds),
        removedIdenticalLhsRows=len(b)-len(bounds),coefficientRoundingUsed=False,
        rule='Identical sparse coefficients; retain the largest original right-hand side',
        selectedSourceRows=selected.tolist())
