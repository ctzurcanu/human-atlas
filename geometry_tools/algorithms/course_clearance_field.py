"""Small cubic displacement fields with literal shared-source dependencies.

Only numerical controls are introduced. Source faces and vertices are never
reconstructed. All members of an explicit section/branch closure share one
translation. The convex basis bounds the displacement everywhere it is used.
"""
import numpy as np
from scipy.interpolate import BSpline
from .tube_sections import arclength


def course_stations(points, axis):
    points,axis=np.asarray(points,float),np.asarray(axis,float)
    arc=arclength(axis);d=np.diff(axis,axis=0)
    t=np.clip(np.einsum('nki,ki->nk',points[:,None]-axis[:-1],d)/np.sum(d*d,axis=1),0,1)
    q=axis[:-1]+t[:,:,None]*d
    k=np.argmin(np.sum((points[:,None]-q)**2,axis=2),axis=1)
    return (arc[k]+t[np.arange(len(points)),k]*np.diff(arc)[k])/arc[-1]


def cubic_carrier(points, axis, groups, *, anchors=None, controls=12):
    """Return raw-point weights; overlap and exact aliases form one closure."""
    p=np.asarray(points,float)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all() or not 6<=controls<=32:
        raise ValueError('Finite source points and 6–32 explicit course controls required')
    _,inverse=np.unique(p,axis=0,return_inverse=True)
    parent=np.arange(int(inverse.max())+1)
    def root(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for group in groups:
        raw=np.asarray(group,int)
        if not len(raw) or np.any(raw<0) or np.any(raw>=len(p)):
            raise ValueError('Explicit nonempty original-source closure required')
        keys=np.unique(inverse[raw]);a=root(int(keys[0]))
        for k in keys[1:]:parent[root(int(k))]=a
    owner=np.array([root(int(i))for i in inverse])
    station=course_stations(p,axis)
    anchor_groups={} if anchors is None else anchors
    anchor_by_owner={}
    for gi,anchor in anchor_groups.items():
        key=int(owner[np.asarray(groups[gi],int)[0]])
        anchor_by_owner.setdefault(key,[]).append(np.asarray(anchor,float))
    for key in np.unique(owner):
        ids=np.flatnonzero(owner==key)
        if len(ids)==1 and key not in anchor_by_owner:continue
        anchor=np.mean(anchor_by_owner[key],axis=0)if key in anchor_by_owner else p[ids].mean(0)
        station[ids]=course_stations(anchor[None],axis)[0]
    degree=3
    knots=np.r_[np.zeros(degree+1),np.linspace(0,1,controls-degree+1)[1:-1],np.ones(degree+1)]
    weights=BSpline.design_matrix(station,knots,degree).toarray()
    if np.min(weights)<-1e-14 or not np.allclose(weights.sum(1),1,rtol=0,atol=1e-13):
        raise ValueError('The carrier must be a nonnegative partition of unity')
    D=np.diff(np.eye(controls),n=2,axis=0)
    precision=np.eye(controls)*.2+10*(D.T@D)
    return weights,station,np.linalg.inv(precision),dict(
        sourceControls=len(p),numericalCourseControls=controls,dependencyClosures=len(np.unique(owner)),
        sourceFacesAdded=0,sourceFacesRemoved=0,maximumPartitionError=float(np.max(abs(weights.sum(1)-1))),
        convexMotionBound=True,exactSourceAliasesCoMove=True,anatomicalAcceptance=False)
