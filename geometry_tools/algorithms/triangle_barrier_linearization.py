"""Analytic oriented triangle-area derivatives for explicit translation fields.

This is a local engineering barrier. It does not certify self-intersection,
finite volume, anatomical interfaces, or continuous deformation injectivity.
"""
import numpy as np
from scipy.sparse import coo_matrix,csr_matrix,hstack


def oriented_triangle_barrier(vertices,faces,reference_normals,area_scale,
                              vertex_influence):
    v=np.asarray(vertices,float);f=np.asarray(faces,int)
    n=np.asarray(reference_normals,float);s=np.asarray(area_scale,float)
    W=csr_matrix(vertex_influence)
    if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or n.shape!=(len(f),3) or s.shape!=(len(f),) or W.shape[0]!=len(v):
        raise ValueError('Literal triangle arrays, reference normals and explicit influence required')
    if not np.isfinite(np.r_[v.ravel(),n.ravel(),s,W.data]).all() or np.any(s<=0) or not np.allclose(np.linalg.norm(n,axis=1),1.,atol=1e-10,rtol=0):
        raise ValueError('Finite nondegenerate source areas and unit reference normals required')
    t=v[f];e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0]
    value=np.einsum('ij,ij->i',n,np.cross(e1,e2))/s
    g1=np.cross(e2,n)/s[:,None];g2=np.cross(n,e1)/s[:,None]
    grad=np.stack([-g1-g2,g1,g2],axis=1)
    axes=[]
    for k in range(3):
        G=coo_matrix((grad[:,:,k].ravel(),(np.repeat(np.arange(len(f)),3),f.ravel())),shape=(len(f),len(v))).tocsr()
        axes.append(G@W)
    # Stack axis-major sparse columns, then permute to profile-major xyz.
    J=hstack(axes,format='csr')
    m=W.shape[1]
    J=J[:,np.arange(3*m).reshape(3,m).T.ravel()]
    return value,J
