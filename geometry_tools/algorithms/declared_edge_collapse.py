"""Collapse one explicitly owned interior edge, retaining a PL distance bound.

This is a remeshing operation, not a tissue deformation or anatomy classifier.
The caller must still review sections, fields, contacts and the complete mesh.
"""
from collections import Counter
import numpy as np


def collapse_edge(vertices,faces,edge,*,protected=(),owned_faces,maximum_error_m):
    p=np.asarray(vertices,float);f=np.asarray(faces,int);a,b=map(int,edge)
    if a==b or min(a,b)<0 or max(a,b)>=len(p) or not np.isfinite(p).all():
        raise ValueError('Two distinct finite source controls required')
    if maximum_error_m<=0 or {a,b}&set(map(int,protected)):
        raise ValueError('Positive error bound and unprotected interior edge required')
    for v in [a,b]:
        if np.sum(np.all(p==p[v],axis=1))!=1:
            raise ValueError('Explicit raw-alias ownership is required before collapse')
    touched=np.flatnonzero(np.any(np.isin(f,[a,b]),axis=1))
    if not set(touched)<=set(map(int,owned_faces)):
        raise ValueError('A complete source-owned vertex star is required')
    shared=np.flatnonzero(np.any(f==a,axis=1)&np.any(f==b,axis=1))
    if len(shared)!=2:raise ValueError('Interior edge must have exactly two incident triangles')
    neighbors=lambda v:set(map(int,f[np.any(f==v,axis=1)].ravel()))-{a,b}
    link=neighbors(a)&neighbors(b)
    opposite=set(map(int,f[shared].ravel()))-{a,b}
    if link!=opposite or len(opposite)!=2:
        raise ValueError('Source edge fails the manifold link condition')
    target=(p[a]+p[b])*.5;error=float(np.linalg.norm(p[[a,b]]-target,axis=1).max())
    if error>maximum_error_m:raise ValueError('Midpoint collapse exceeds the declared PL error bound')
    q=p.copy();q[[a,b]]=target
    keep=np.setdiff1d(np.arange(len(f)),shared);ff=f[keep]
    tri=q[ff];area=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)
    if np.any(area<=1e-24):raise ValueError('A retained triangle became degenerate')
    # Each removed triangle maps onto an edge still covered by the target mesh.
    physical,inv=np.unique(q,axis=0,return_inverse=True)
    edges=Counter(tuple(sorted(map(int,[x,y])))for t in inv[ff]for x,y in zip(t,np.roll(t,-1)))
    for face in f[shared]:
        image=np.unique(inv[face])
        if len(image)!=2 or tuple(image)not in edges:
            raise ValueError('Collapsed triangle image is not covered by a retained edge')
    proof=dict(complete=True,edge=[a,b],removedSourceFaceIds=shared.tolist(),
        retainedSourceFaceIds=keep.tolist(),touchedSourceFaceIds=touched.tolist(),
        midpoint_m=target.tolist(),maximumTwoWayPiecewiseLinearDistance_m=error,
        linkConditionPassed=True,collapsedTriangleImagesCovered=True,
        sourceVertexCorrespondence='Identity except the two explicit midpoint images',
        boundMeaning='Continuous two-way PL surface-distance bound, not tissue strain or anatomical acceptance',
        anatomicalAcceptance=False)
    return q,ff,proof
