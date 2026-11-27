"""Proper source-profile rotation about a literal protected contact edge."""
import numpy as np
from scipy.spatial.transform import Rotation


def contact_edge_pivot(points, fixed_indices):
    p=np.asarray(points,float);fixed=np.asarray(fixed_indices,int)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()or fixed.ndim!=1 or np.any(fixed<0)or np.any(fixed>=len(p)):
        raise ValueError('Finite source profile and valid fixed source indices required')
    q=np.unique(p[fixed],axis=0)
    if len(q)<2:raise ValueError('A literal contact edge requires two distinct fixed points')
    origin=q[0];direction=q[-1]-origin;direction/=np.linalg.norm(direction)
    error=np.linalg.norm(np.cross(q-origin,direction),axis=1).max(initial=0)
    if error>1e-12:raise ValueError('Noncollinear fixed points lock the complete rigid profile')
    return origin,direction


def pose_about_contact_edge(points, fixed_indices, angle_radians, *, maximum_angle_radians):
    p=np.asarray(points,float);fixed=np.asarray(fixed_indices,int)
    if not np.isfinite(angle_radians)or not np.isfinite(maximum_angle_radians)or maximum_angle_radians<=0 or abs(angle_radians)>maximum_angle_radians:
        raise ValueError('A finite angle inside the predeclared bound required')
    center,axis=contact_edge_pivot(p,fixed);rotation=Rotation.from_rotvec(axis*angle_radians).as_matrix();raw=(p-center)@rotation.T+center
    error=float(np.linalg.norm(raw[fixed]-p[fixed],axis=1).max(initial=0))
    if error>3e-16:raise ValueError('Pivot roundoff exceeds the rigid source replay tolerance')
    result=raw.copy();result[fixed]=p[fixed]
    return result,rotation,center,dict(complete=True,properRotation=True,angleRadians=float(angle_radians),
                                      maximumLiteralRestoreRoundoff_m=error,calibreChanged=False,
                                      anatomicalAcceptance=False)
