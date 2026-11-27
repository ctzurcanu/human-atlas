"""Pose an existing parent annulus opposite a declared child's inlet direction."""
import numpy as np


def pose_parent_collar(source_mouth, source_ring, target_mouth, child_next, *, gauge_tolerance_m=1e-8):
    arrays=[np.asarray(x,float)for x in (source_mouth,source_ring,target_mouth,child_next)]
    a,ring,b,child=arrays
    if any(x.ndim!=2 or x.shape[1]!=3 or len(x)<3 or not np.isfinite(x).all()for x in arrays):
        raise ValueError('Four finite, explicitly matched annular point sets are required')
    if a.shape!=b.shape or not 0<gauge_tolerance_m<=1e-6:
        raise ValueError('Matched source/target mouths and explicit export tolerance required')
    distance=lambda x:np.linalg.norm(x[:,None]-x[None],axis=2)
    gauge=float(np.max(abs(distance(a)-distance(b))))
    if gauge>gauge_tolerance_m:
        raise ValueError('Restore the original mouth shape before constructing its collar')
    ac,bc=a.mean(0),b.mean(0)
    u,s,vt=np.linalg.svd((a-ac).T@(b-bc))
    if s[1]<1e-20:raise ValueError('A collinear mouth cannot determine an annular pose')
    vt[-1]*=np.linalg.det(vt.T@u.T);R=vt.T@u.T
    delta=child.mean(0)-bc;length=float(np.linalg.norm(delta))
    source_length=float(np.linalg.norm(ring.mean(0)-ac))
    if min(length,source_length)<=1e-8:
        raise ValueError('A collar requires nonzero actual parent and child intervals')
    direction=delta/length
    centre=bc-direction*source_length
    result=(ring-ring.mean(0))@R.T+centre
    error=float(np.max(abs(distance(result)-distance(ring))))
    if error>1e-12:raise ValueError('Parent collar lost source annular calibre')
    return result,dict(sourceInterval_m=source_length,childInterval_m=length,
        originalMouthGaugeError_m=gauge,collarGaugeError_m=error,
        parentChildDirectionDot=float((centre-bc)@direction/source_length),
        syntheticFaces=0,anatomicalAcceptance=False)
