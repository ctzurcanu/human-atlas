"""Proper rigid source-group rotations about a fixed sampled anatomical axis."""
import numpy as np
from scipy.spatial.transform import Rotation


def phase_rigid_groups(vertices,labels,axis,tangents,angle_degrees,protected_raw,fade_distance_m=.015):
    v=np.asarray(vertices,float);labels=np.asarray(labels,int);axis=np.asarray(axis,float);t=np.asarray(tangents,float)
    if labels.shape!=(len(v),)or axis.ndim!=2 or axis.shape[1]!=3 or t.shape!=axis.shape:
        raise ValueError('Explicit original group ownership and matching fixed axis/tangents required')
    if not np.isfinite(angle_degrees)or abs(angle_degrees)>30 or not 0<fade_distance_m<=.05:
        raise ValueError('Explicit phase within 30 degrees and bounded fade distance required')
    if not np.isfinite(np.r_[v.ravel(),axis.ravel(),t.ravel()]).all()or not np.allclose(np.linalg.norm(t,axis=1),1.,atol=1e-12,rtol=0):
        raise ValueError('Finite source coordinates and unit axis tangents required')
    _,owner=np.unique(labels,return_inverse=True);count=int(owner.max())+1
    centres=np.array([v[owner==k].mean(0)for k in range(count)])
    protected_raw=np.asarray(protected_raw,int)
    if not len(protected_raw)or protected_raw.min()<0 or protected_raw.max()>=len(v):
        raise ValueError('Actual protected source mouth controls required')
    protected=np.unique(owner[protected_raw]);arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(axis,axis=0),axis=1))]
    if np.any(np.diff(arc)<=0):raise ValueError('Distinct ordered anatomical axis samples required')
    nearest=np.argmin(np.linalg.norm(centres[:,None]-axis[None],axis=2),axis=1)
    start=float(arc[nearest[protected]].mean());u=np.clip((arc[nearest]-start)/fade_distance_m,0,1)
    fade=u*u*u*(10-15*u+6*u*u);fade[protected]=0.
    R=Rotation.from_rotvec(t[nearest]*np.radians(angle_degrees)*fade[:,None]).as_matrix()
    q=axis[nearest][owner]+np.einsum('nij,nj->ni',R[owner],v-axis[nearest][owner])
    fixed=np.isin(owner,protected);q[fixed]=v[fixed]
    if angle_degrees==0:q=v.copy()
    return q,R[owner],dict(angleDegrees=float(angle_degrees),rigidSourceGroups=count,protectedSourceGroups=len(protected),
                           fadeDistance_m=fade_distance_m,sourceCoordinatesWelded=False,axisMoved=False,
                           sampledAxisOnly=True,anatomicalAcceptance=False)
