"""Project an explicit source triangle onto its original hinge circle.

This preserves its intrinsic edge lengths while allowing a local fold angle
to change. It neither identifies an anatomical joint nor accepts a mesh.
"""
import numpy as np


def project_hinge_isometry(source, proposed, vertex, edge):
    p=np.asarray(source,float);q=np.asarray(proposed,float)
    edge=np.asarray(edge,int);vertex=int(vertex)
    if p.shape!=q.shape or p.ndim!=2 or p.shape[1]!=3 or edge.shape!=(2,):
        raise ValueError('Corresponding3D controls and one explicit hinge edge required')
    if len(set([vertex,*edge]))!=3 or min(vertex,*edge)<0 or max(vertex,*edge)>=len(p):
        raise ValueError('Three distinct existing controls required')
    if not np.isfinite(np.r_[p.ravel(),q.ravel()]).all() or not np.array_equal(p[edge],q[edge]):
        raise ValueError('Finite controls and literal fixed hinge endpoints required')
    a,b=p[edge];axis=b-a;length=np.linalg.norm(axis)
    if length<=1e-14:raise ValueError('Nonzero source hinge edge required')
    axis/=length;centre=a+axis*np.dot(p[vertex]-a,axis)
    radius=np.linalg.norm(p[vertex]-centre)
    radial=q[vertex]-centre;radial-=axis*np.dot(radial,axis);r=np.linalg.norm(radial)
    if min(radius,r)<=1e-14:raise ValueError('Nondegenerate source and proposed hinge radii required')
    result=centre+radius*radial/r
    np.testing.assert_allclose(np.linalg.norm(result-p[edge],axis=1),np.linalg.norm(p[vertex]-p[edge],axis=1),rtol=1e-10,atol=1e-14)
    return result


def intrinsic_facet_strains(source, proposed, faces):
    """Two singular strains per corresponding nondegenerate triangle."""
    x=np.asarray(source,float)[faces];y=np.asarray(proposed,float)[faces]
    E=np.swapaxes(x[:,1:]-x[:,:1],1,2);F=np.swapaxes(y[:,1:]-y[:,:1],1,2)
    if np.any(np.linalg.norm(np.cross(E[:,:,0],E[:,:,1]),axis=1)<=1e-24):
        raise ValueError('Nondegenerate source facets required')
    # Restrict the ambient affine map to the source tangent chart.
    u,_,_=np.linalg.svd(E,full_matrices=False)
    mapped=F@np.linalg.pinv(E)@u
    return np.linalg.svd(mapped,compute_uv=False)-1
