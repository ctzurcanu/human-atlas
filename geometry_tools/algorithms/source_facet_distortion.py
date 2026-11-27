"""Sparse affine tangent changes for explicitly assigned translation stations."""
import numpy as np
from scipy import sparse


def facet_tangent_rows(mesh,stations,station_count):
    """Return six rows per triangle: vec((E_new-E_old) R).

    E_old R is an orthonormal source tangent frame, obtained by a two-column
    SVD. Bounding each six-row Frobenius norm by epsilon<1 bounds all tangent
    singular values in [1-epsilon,1+epsilon]. This is a source-facet shape
    condition, not a volumetric Jacobian or anatomical tissue acceptance.
    """
    stations=np.asarray(stations,int)
    assert stations.shape==(len(mesh.vertices),)and np.all(stations>=-1)and np.all(stations<station_count)
    t=mesh.triangles;edges=np.stack([t[:,1]-t[:,0],t[:,2]-t[:,0]],axis=2)
    u,s,v=np.linalg.svd(edges,full_matrices=False)
    assert np.isfinite(s).all()and np.all(s[:,1]>0),'Source-degenerate facets need separate repair'
    R=np.swapaxes(v,1,2)/s[:,None,:]
    coefficients=np.stack([-R.sum(1),R[:,0],R[:,1]],axis=1)
    rows=[];cols=[];data=[]
    for corner in range(3):
        group=stations[mesh.faces[:,corner]];live=group>=0
        for column in range(2):
            for component in range(3):
                rows.append(np.flatnonzero(live)*6+column*3+component)
                cols.append(group[live]*3+component)
                data.append(coefficients[live,corner,column])
    C=sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(t)*6,station_count*3)).tocsr()
    C.sum_duplicates();C.eliminate_zeros()
    return C,dict(sourceFaces=len(t),rows=C.shape[0],stationCount=station_count,
                  minimumSourceTangentSingularValue_m=float(s[:,1].min()),
                  rule='Frobenius norm of affine source tangent perturbation',anatomicalAcceptance=False)
