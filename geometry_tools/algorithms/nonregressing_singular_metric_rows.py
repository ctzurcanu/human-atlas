"""Conservative incremental cones bounded by inherited original singular strain."""
import numpy as np
from scipy import sparse
from .incremental_source_metric_rows import incremental_source_metric_rows

def nonregressing_singular_metric_rows(reference,current,stations,station_count,*,nominal_budget=.075):
    C,offset,proof=incremental_source_metric_rows(reference,current,stations,station_count,nominal_budget=nominal_budget)
    E=(reference.triangles[:,1:]-reference.triangles[:,:1]).transpose(0,2,1)
    Q=(current.triangles[:,1:]-current.triangles[:,:1]).transpose(0,2,1)
    _,s,v=np.linalg.svd(E,full_matrices=False);R=np.swapaxes(v,1,2)/s[:,None,:]
    strain=np.abs(np.linalg.svd(Q@R,compute_uv=False)-1).max(1)
    budget=np.maximum(nominal_budget,strain)
    factor=np.asarray(proof['sourceFacetBudgets'])/budget
    C=sparse.diags(np.repeat(factor,6))@C
    offset=offset*factor[:,None]
    return C,offset,dict(proof,sourceFacetBudgets=budget.tolist(),
        normalizationFactors=(nominal_budget/budget).tolist(),inheritedOriginalSingularStrain=strain.tolist(),
        rule='A proper-original-frame Frobenius cone bounded by max(nominal, inherited maximum singular strain) guarantees non-worsening original singular strain',
        sourceReferenceNeverReset=True)
