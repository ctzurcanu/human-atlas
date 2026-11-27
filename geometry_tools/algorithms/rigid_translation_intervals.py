"""Conservative SAT time intervals for two rigidly translating triangles.

All corners enter each support bound. Axes stay fixed under a common proper
translation, so their overlap conditions are linear in travelled distance.
The complement of the returned interval has a separating-plane witness.
This does not establish an anatomical corridor, closed volume or deforming
collar clearance. Floating-point bounds remain explicit.
"""
import numpy as np


def triangle_overlap_interval(first, second, direction, *, maximum_distance_m,
                              clearance_m=0.):
    a,b,d=(np.asarray(x,dtype=float) for x in (first,second,direction))
    if a.shape!=(3,3) or b.shape!=(3,3) or d.shape!=(3,) or not all(np.isfinite(x).all() for x in (a,b,d)):
        raise ValueError('Two finite actual triangles and one direction required')
    if not np.isfinite(maximum_distance_m) or maximum_distance_m<=0 or not np.isfinite(clearance_m) or clearance_m<0:
        raise ValueError('Finite positive travel cap and nonnegative clearance required')
    length=float(np.linalg.norm(d))
    if not length:raise ValueError('Nonzero direction required')
    d=d/length;origin=(a.mean(0)+b.mean(0))*.5
    world_scale=max(float(np.abs(a).max()),float(np.abs(b).max()))
    a=a-origin;b=b-origin
    ea=np.roll(a,-1,axis=0)-a;eb=np.roll(b,-1,axis=0)-b
    na=np.cross(ea[0],ea[1]);nb=np.cross(eb[0],eb[1])
    if np.linalg.norm(na)==0 or np.linalg.norm(nb)==0:
        raise ValueError('Degenerate triangles require separate source classification')
    axes=np.vstack([na,nb,np.cross(ea[:,None],eb[None]).reshape(-1,3),
                    np.cross(na,ea),np.cross(nb,eb)])
    lengths=np.linalg.norm(axes,axis=1);axes=axes[lengths>0]/lengths[lengths>0,None]
    pa=a@axes.T;pb=b@axes.T;speed=axes@d
    arithmetic=256*np.finfo(float).eps*(world_scale+float(np.abs(origin).max())+
        float(np.abs(a).max())+float(np.abs(b).max())+maximum_distance_m+clearance_m)
    lower=pb.min(0)-pa.max(0)-clearance_m-arithmetic
    upper=pb.max(0)-pa.min(0)+clearance_m+arithmetic
    lo,hi=0.,float(maximum_distance_m)
    for low,high,slope in zip(lower,upper,speed):
        if slope>0:lo=max(lo,float(low/slope));hi=min(hi,float(high/slope))
        elif slope<0:lo=max(lo,float(high/slope));hi=min(hi,float(low/slope))
        elif low>0 or high<0:lo,hi=1.,0.;break
        if lo>hi:break
    return dict(forbidden_interval_m=None if lo>hi else [lo,hi],
        maximum_distance_m=float(maximum_distance_m),unit_direction=d.tolist(),
        requested_projected_clearance_m=float(clearance_m),support_arithmetic_bound_m=float(arithmetic),
        actual_corner_supports=True,axes_tested=len(axes),rigid_translation_only=True,
        anatomical_corridor_established=False)


def clear_translation_intervals(intervals, *, maximum_distance_m):
    """Return open interiors outside the union of conservative blocked ranges."""
    if not np.isfinite(maximum_distance_m) or maximum_distance_m<=0:raise ValueError('Positive finite travel cap required')
    blocked=[]
    for interval in intervals:
        if interval is None:continue
        values=np.asarray(interval,float)
        if values.shape!=(2,) or not np.isfinite(values).all() or values[0]>values[1]:
            raise ValueError('Finite ordered intervals required')
        lo=max(0.,float(values[0]));hi=min(float(maximum_distance_m),float(values[1]))
        if lo<=hi:blocked.append([lo,hi])
    merged=[]
    for lo,hi in sorted(blocked):
        if merged and lo<=merged[-1][1]:merged[-1][1]=max(hi,merged[-1][1])
        else:merged.append([lo,hi])
    result=[];start=0.
    for lo,hi in merged:
        if start<lo:result.append([start,lo])
        start=max(start,hi)
    if start<maximum_distance_m:result.append([start,float(maximum_distance_m)])
    return dict(clear_open_interiors_m=result,blocked_closed_intervals_m=merged,
        interval_endpoints_are_not_clearance_certificates=True)
