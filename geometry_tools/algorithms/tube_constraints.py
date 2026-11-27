"""Monotone section translation with hard per-span strain caps.

Only scalar translation along a fixed direction is redistributed. A caller
must derive caps from actual wall/packing predicates and recheck saved geometry.
"""
import numpy as np
from .tube_sections import TubeGeometryError


def capped_translation_transition(desired_offsets, maximum_span_change):
    desired=np.asarray(desired_offsets,float);caps=np.asarray(maximum_span_change,float)
    if desired.ndim!=1 or not len(desired) or caps.shape!=(len(desired)-1,) or not np.isfinite(desired).all() or not np.isfinite(caps).all() or np.any(caps<0):
        raise TubeGeometryError('finite offsets and nonnegative per-span caps required')
    demand=-np.diff(desired);total=desired[0]-desired[-1]
    if np.any(demand<0) or total<0:raise TubeGeometryError('offset transition must be monotone decreasing')
    if not np.isfinite(total) or not np.isfinite(caps.sum()) or caps.sum()<total:raise TubeGeometryError('available span caps cannot meet the pinned endpoint translation')
    if total==0:return desired.copy()
    weights=np.maximum(demand,1e-12)
    lo,hi=0.,1.
    while np.minimum(caps,hi*weights).sum()<total:hi*=2
    for _ in range(80):
        mid=(lo+hi)/2
        if np.minimum(caps,mid*weights).sum()<total:lo=mid
        else:hi=mid
    steps=np.minimum(caps,hi*weights);offsets=np.r_[desired[0],desired[0]-np.cumsum(steps)]
    offsets[-1]=desired[-1]
    return offsets
