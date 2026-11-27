"""Incremental bend-frame rows with immutable-source, per-facet budgets."""
import numpy as np
from scipy import sparse
from .corotated_source_tangent import corotated_tangent_rows

def incremental_source_metric_rows(reference,current,stations,station_count,*,nominal_budget=.075):
    C,offset,proof=corotated_tangent_rows(reference,current,stations,station_count)
    t,q=reference.triangles,current.triangles
    E=(t[:,1:]-t[:,:1]).transpose(0,2,1);Q=(q[:,1:]-q[:,:1]).transpose(0,2,1)
    _,s,v=np.linalg.svd(E,full_matrices=False);R=np.swapaxes(v,1,2)/s[:,None,:]
    offset=offset+((Q-E)@R).transpose(0,2,1).reshape(-1,6)
    inherited=np.linalg.norm(offset,axis=1);budgets=np.maximum(nominal_budget,inherited)
    assert np.isfinite(budgets).all() and np.all(budgets<1)
    factors=nominal_budget/budgets
    scaled=sparse.diags(np.repeat(factors,6))@C
    return scaled,offset*factors[:,None],dict(proof,
        incrementalMotionAboutCurrentPose=True,referenceGeometryUnchanged=True,
        nominalBudget=nominal_budget,sourceFacetBudgets=budgets.tolist(),
        inheritedCorotatedResiduals=inherited.tolist(),normalizationFactors=factors.tolist(),
        inheritedAboveNominalFacets=int((inherited>nominal_budget).sum()),
        rule='Each original-source residual budget is max(nominal, inherited); it is never added to inherited distortion')
