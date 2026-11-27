"""Bound entire source surfaces using distance Lipschitz continuity.

Every original triangle is covered, without nearest-vertex assumptions.
This certifies geometric closeness only, never anatomical identity.
"""
import time
import numpy as np
from .section_course_clearance import closest_triangles

def directed_surface_bound(first,second,*,maximum_distance_m,wall_seconds=6.,maximum_points=200000,maximum_levels=9):
    a,b=map(lambda v:np.asarray(v,float),(first,second))
    if any(t.ndim!=3 or t.shape[1:]!=(3,3) or not len(t) or not np.isfinite(t).all() for t in [a,b]):raise ValueError('Two finite nonempty triangle surfaces required')
    if not np.isfinite(maximum_distance_m) or maximum_distance_m<=0 or wall_seconds<=0 or maximum_levels<1:raise ValueError('Positive finite source coverage limits required')
    start=time.monotonic();world=max(float(np.abs(a).max()),float(np.abs(b).max()));origin=(a.mean((0,1))+b.mean((0,1)))*.5;a=a-origin;b=b-origin
    arithmetic=1024*np.finfo(float).eps*(world+float(np.abs(origin).max())+float(np.abs(a).max())+float(np.abs(b).max()))
    pending=a;upper=0.;lower=0.;points=0;levels=[]
    for level in range(maximum_levels):
        if time.monotonic()-start>wall_seconds:raise TimeoutError('Finite surface coverage deadline')
        centers=pending.mean(1);points+=len(centers)
        if points>maximum_points:raise ValueError('Finite surface coverage sample budget exceeded')
        distances=[]
        for offset in range(0,len(centers),512):
            if time.monotonic()-start>wall_seconds:raise TimeoutError('Finite surface coverage deadline')
            p=centers[offset:offset+512];q,_=closest_triangles(p,b,batch_size=128);distances.append(np.linalg.norm(p-q,axis=1))
        distance=np.concatenate(distances);radius=np.linalg.norm(pending-centers[:,None],axis=2).max(1);bound=distance+radius+arithmetic;lower=max(lower,float(distance.max()-arithmetic))
        if lower>maximum_distance_m:raise ValueError(f'Source surfaces exceed the declared bound: witnessed {lower:g}m')
        accepted=bound<=maximum_distance_m
        if accepted.any():upper=max(upper,float(bound[accepted].max()))
        levels.append(dict(level=level,triangles=len(pending),accepted=int(accepted.sum())))
        if accepted.all():return dict(complete=True,passed=True,directedDistanceUpperBound_m=upper,witnessedDistanceLowerBound_m=max(lower,0.),sourceTriangles=len(a),targetTriangles=len(b),queriedPoints=points,levels=levels,arithmeticBound_m=arithmetic,nearestPointMethod='Exhaustive Float64 actual-triangle closest points',surfaceCover='Distance to a fixed closed triangle set is1-Lipschitz; centroid distance plus the complete subtriangle radius covers every point.',elapsedSeconds=time.monotonic()-start,anatomicalIdentityInferred=False)
        t=pending[~accepted];x=(t[:,0]+t[:,1])*.5;y=(t[:,1]+t[:,2])*.5;z=(t[:,2]+t[:,0])*.5
        pending=np.concatenate([np.stack([t[:,0],x,z],axis=1),np.stack([x,t[:,1],y],axis=1),np.stack([z,y,t[:,2]],axis=1),np.stack([x,y,z],axis=1)])
    raise ValueError('Source coverage unresolved at finite subdivision limit')

def symmetric_surface_bound(first,second,**limits):
    a=directed_surface_bound(first,second,**limits);b=directed_surface_bound(second,first,**limits)
    return dict(complete=True,passed=True,symmetricDistanceUpperBound_m=max(a['directedDistanceUpperBound_m'],b['directedDistanceUpperBound_m']),directions=[a,b],anatomicalIdentityInferred=False)
