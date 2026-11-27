"""Explicit cross-course separation axes within a supplied tissue layer.

This geometry helper does not establish an anatomical course, blood junction,
or clearance. The caller must supply reviewed layer normals, source-course
tangents and a consistent side reference; unresolved directions are rejected.
"""
import numpy as np


def layer_cross_course_axes(layer_normals,first_tangents,second_tangents,*,side_reference):
    n,a,b=[np.asarray(v,float).copy() for v in [layer_normals,first_tangents,second_tangents]]
    if n.ndim!=2 or n.shape[1]!=3 or a.shape!=n.shape or b.shape!=n.shape:
        raise ValueError('Matching explicit layer normals and source-course tangents required')
    side=np.broadcast_to(np.asarray(side_reference,float),n.shape)
    if not np.isfinite(np.r_[n.ravel(),a.ravel(),b.ravel(),side.ravel()]).all():
        raise ValueError('Finite direction fields required')
    for v in [n,a,b]:
        lengths=np.linalg.norm(v,axis=1)
        if np.any(lengths<1e-12):raise ValueError('Unresolved zero direction')
        v/=lengths[:,None]
    b*=np.where(np.einsum('ij,ij->i',a,b)<0,-1.,1.)[:,None]
    course=a+b;course/=np.linalg.norm(course,axis=1)[:,None]
    tangent=course-np.einsum('ij,ij->i',course,n)[:,None]*n
    lengths=np.linalg.norm(tangent,axis=1)
    if np.any(lengths<1e-3):
        raise ValueError('Course enters the layer nearly normally; no reviewed in-layer lateral axis')
    tangent/=lengths[:,None]
    axis=np.cross(n,tangent);axis/=np.linalg.norm(axis,axis=1)[:,None]
    orientation=np.einsum('ij,ij->i',axis,side)
    if np.any(abs(orientation)<1e-8):
        raise ValueError('Supplied side reference does not resolve lateral orientation')
    axis*=np.where(orientation<0,-1.,1.)[:,None]
    return axis,dict(complete=True,axes=len(axis),sourceGeometryModified=False,
        maximumLayerNormalComponent=float(abs(np.einsum('ij,ij->i',axis,n)).max(initial=0)),
        maximumCourseComponent=float(abs(np.einsum('ij,ij->i',axis,course)).max(initial=0)),
        actualClearanceEstablished=False,anatomicalCourseEstablished=False,
        sideChosenOnlyFromExplicitReference=True)
