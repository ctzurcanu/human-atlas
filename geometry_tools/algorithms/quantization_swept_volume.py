"""Conservative volume bound for moving triangular surface controls.

Every intermediate triangle of linear vertex interpolation lies inside the
convex hull of its six endpoint controls. An oriented box containing those six
points therefore bounds its entire sweep. Summing boxes bounds the volume in
which an oriented closed surface can change its winding. This is a numerical
storage diagnostic, never an anatomical clearance allowance.
"""
import numpy as np

def swept_volume_bound(before,after,faces):
    before=np.asarray(before,float);after=np.asarray(after,float);faces=np.asarray(faces)
    if before.shape!=after.shape or before.ndim!=2 or before.shape[1]!=3 or not np.isfinite(before).all()or not np.isfinite(after).all():raise ValueError('Finite corresponding positions required')
    if faces.ndim!=2 or faces.shape[1]!=3 or not np.issubdtype(faces.dtype,np.integer)or not len(faces)or faces.min()<0 or faces.max()>=len(before):raise ValueError('Valid corresponding triangle indices required')
    a=before[faces];b=after[faces];u=a[:,1]-a[:,0];n=np.cross(u,a[:,2]-a[:,0]);ul=np.linalg.norm(u,axis=1);nl=np.linalg.norm(n,axis=1)
    if np.any(ul==0)or np.any(nl==0):raise ValueError('Degenerate source triangles do not define a frame')
    u=u/ul[:,None];n=n/nl[:,None];v=np.cross(n,u);axes=np.stack([u,v,n],axis=2);six=np.concatenate([a,b],axis=1)-a[:,0,None];q=np.einsum('tki,tij->tkj',six,axes);width=np.ptp(q,axis=1)
    # Compensate double-precision projection arithmetic from local coordinates;
    # this deliberately dominates the few dot/add roundings in the boxes.
    guard=128*np.finfo(float).eps*np.max(np.linalg.norm(six,axis=2),axis=1);width+=2*guard[:,None]
    return dict(maximumControlMotionMeters=float(np.linalg.norm(after-before,axis=1).max()),sweptVolumeUpperBoundMetersCubed=float(np.prod(width,axis=1).sum()),triangles=len(faces),method='Six-control swept triangle boxes, in local orthonormal frames with explicit arithmetic guard')
