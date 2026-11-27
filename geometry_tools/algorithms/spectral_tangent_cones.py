"""Exact operator-norm cones for inherited two-axis source distortion.

Most facets keep the cheaper sufficient Frobenius cone. If that cone rejects
the current, otherwise permitted inherited source strain, use the exact LMI
[[e I3, D], [D.T, e I2]] >= 0 instead. Its eigenvalues establish ||D||_2 <= e.
Whitening, proper bend frames and original facet budgets are caller inputs;
this operation never resets references or increases a physical shape budget.
"""
import clarabel
import numpy as np
from scipy import sparse


def mixed_tangent_cones(matrix,offset,*,bound,reserve,extra_columns=0):
    C=sparse.csr_matrix(matrix);o=np.asarray(offset,float).reshape(-1,6);count=len(o)
    if C.shape[0]!=6*count or not np.isfinite(C.data).all() or not np.isfinite(o).all() or not 0<=reserve<bound<1:
        raise ValueError('Finite original-frame six-row facet operators and explicit bounds required')
    mask=np.linalg.norm(o,axis=1)>bound-reserve+1e-10
    regular=np.flatnonzero(~mask);special=np.flatnonzero(mask);blocks=[];bounds=[];cones=[]
    if len(regular):
        a=C[(regular[:,None]*6+np.arange(6)).ravel()].tocoo()
        rows=a.row//6*7+1+a.row%6
        D=sparse.coo_matrix((-.001*a.data,(rows,a.col)),shape=(len(regular)*7,C.shape[1]+extra_columns)).tocsc()
        blocks.append(D);bounds.append(np.c_[np.full(len(regular),bound-reserve),o[regular]].ravel());cones.extend(clarabel.SecondOrderConeT(7) for _ in regular)
    if len(special):
        a=C[(special[:,None]*6+np.arange(6)).ravel()].tocoo();positions=np.array([6,7,8,10,11,12])
        rows=a.row//6*15+positions[a.row%6]
        D=sparse.coo_matrix((-.001*np.sqrt(2)*a.data,(rows,a.col)),shape=(len(special)*15,C.shape[1]+extra_columns)).tocsc()
        # Reserve unused budget, not an improvement that the inherited state
        # never promised. In particular, a facet already exactly at its
        # original-source budget must retain the feasible no-change state.
        # Never enlarge the caller's physical budget, including for an input
        # that already exceeds it.
        inherited=np.linalg.svd(o[special].reshape(-1,2,3).transpose(0,2,1),compute_uv=False)[:,0]
        limits=np.maximum(bound-reserve,np.minimum(bound,inherited))
        r=np.zeros((len(special),15));r[:,[0,2,5,9,14]]=limits[:,None];r[:,positions]=np.sqrt(2)*o[special]
        blocks.append(D);bounds.append(r.ravel());cones.extend(clarabel.PSDTriangleConeT(5) for _ in special)
    return sparse.vstack(blocks,format='csc'),np.concatenate(bounds),cones,mask,dict(complete=True,
        frobeniusFacets=len(regular),exactSpectralFacets=len(special),originalReferencesAndBudgetsUnchanged=True,
        numericalReserveConsumesOnlyUnusedInheritedBudget=True,
        packedUpperTriangularColumnOrder=True,offDiagonalScale=float(np.sqrt(2)),
        documentation='https://clarabel.org/stable/api_cone_types/')


def mixed_tangent_residuals(deformation,spectral_mask):
    d=np.asarray(deformation,float).reshape(-1,6);mask=np.asarray(spectral_mask,bool)
    if mask.shape!=(len(d),):raise ValueError('One preserved cone classification per original facet required')
    norms=np.linalg.norm(d,axis=1)
    if mask.any():norms[mask]=np.linalg.svd(d[mask].reshape(-1,2,3).transpose(0,2,1),compute_uv=False)[:,0]
    return norms
