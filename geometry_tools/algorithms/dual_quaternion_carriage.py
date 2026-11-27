"""Blend explicit proper bone poses without linear rotation shrinkage.

Kavan et al., Geometric Skinning with Approximate Dual Quaternion Blending
(2008), https://users.cs.utah.edu/~ladislav/kavan08geometric/
This is a placement tool, not a collision or anatomical acceptance certificate.
Quaternions use scipy's xyzw convention. Inputs contain no inferred anatomy.
"""
import numpy as np
from scipy.spatial.transform import Rotation


def _product(a,b):
    return np.concatenate((a[...,-1,None]*b[...,:3]+b[...,-1,None]*a[...,:3]
                           +np.cross(a[...,:3],b[...,:3]),
                           (a[...,-1]*b[...,-1]-np.sum(a[...,:3]*b[...,:3],axis=-1))[...,None]),axis=-1)


def blend_proper_poses(points,weights,rotations,translations):
    """Return posed points and the actual blended proper rotations.

    Nonnegative weights must sum to one. Each query independently chooses
    its largest-weight quaternion for antipodal alignment; an unrelated
    distant bone cannot force the sign at a joint.
    """
    p,w,R,t=map(lambda x:np.asarray(x,float),(points,weights,rotations,translations))
    if p.ndim!=2 or p.shape[1:]!=(3,) or R.ndim!=3 or R.shape[1:]!=(3,3):
        raise ValueError('Explicit query points and 3D bone rotations required')
    if t.shape!=(len(R),3) or w.shape!=(len(p),len(R)):
        raise ValueError('One translation per bone and one weight per query/bone required')
    if not all(np.isfinite(x).all() for x in (p,w,R,t)) or np.any(w<0) or not np.allclose(w.sum(1),1,atol=1e-12,rtol=0):
        raise ValueError('Finite poses and normalized nonnegative weights required')
    if not np.allclose(R.transpose(0,2,1)@R,np.eye(3),atol=1e-10,rtol=0) or np.any(np.linalg.det(R)<1-1e-10):
        raise ValueError('Reflection, shear and scaling are not proper rigid poses')
    real=Rotation.from_matrix(R).as_quat()
    dual=.5*_product(np.c_[t,np.zeros(len(t))],real)
    reference=real[np.argmax(w,axis=1)]
    sign=np.where(reference@real.T<0,-1.,1.)
    a=(w*sign)@real;b=(w*sign)@dual
    norm=np.linalg.norm(a,axis=1)
    if np.any(norm<1e-10):raise ValueError('Quaternion blend is singular')
    a/=norm[:,None];b/=norm[:,None]
    b-=a*np.sum(a*b,axis=1)[:,None]
    conjugate=a*np.array([-1.,-1.,-1.,1.])
    displacement=2*_product(b,conjugate)[:,:3]
    rotation=Rotation.from_quat(a)
    return rotation.apply(p)+displacement,rotation.as_matrix()
