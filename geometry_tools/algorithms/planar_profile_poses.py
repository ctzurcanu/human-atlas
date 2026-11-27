"""Whole-profile proper-pose tests against explicitly typed planar geometry.

Solid envelopes, open-sheet barriers and a containing enclosure have different
predicates. No outline is repaired, closed or promoted into anatomical tissue.
"""
import numpy as np
import shapely


def profile_pose_mask(points, poses, *, solid=None, barrier=None, container=None,
                      side_exclusion=None, margin_m=.0001):
    p,t=np.asarray(points,float),np.asarray(poses,float)
    if(p.ndim!=2 or p.shape[1]!=2 or len(p)<3 or t.ndim!=2 or t.shape[1]!=3
       or not np.isfinite(np.r_[p.ravel(),t.ravel()]).all()
       or not np.isfinite(margin_m)or margin_m<0):
        raise ValueError('Finite source outline, X/Z/roll poses and nonnegative margin required')
    source=shapely.polygons(p)
    if not shapely.is_valid(source)or shapely.area(source)<=0:raise ValueError('Actual simple positive-area source profile required')
    c,s=np.cos(t[:,2]),np.sin(t[:,2])
    placed=np.stack([c[:,None]*p[None,:,0]+s[:,None]*p[None,:,1],
                     -s[:,None]*p[None,:,0]+c[:,None]*p[None,:,1]],axis=2)+t[:,None,:2]
    shapes=shapely.polygons(placed);clear=np.ones(len(t),bool);gap=np.full(len(t),np.inf)
    for obstacle in [solid,barrier,side_exclusion]:
        if obstacle is None or shapely.is_empty(obstacle):continue
        distance=shapely.distance(shapes,obstacle);clear&=shapely.disjoint(shapes,obstacle)&(distance>=margin_m);gap=np.minimum(gap,distance)
    if container is not None:
        if not shapely.is_valid(container)or shapely.is_empty(container):raise ValueError('Actual qualified enclosure required')
        distance=shapely.distance(shapes,shapely.boundary(container));clear&=shapely.contains(container,shapes)&(distance>=margin_m);gap=np.minimum(gap,distance)
    relative=np.abs(shapely.area(shapes)/shapely.area(source)-1)
    if np.any(relative>1e-10):raise ValueError('Proper-pose source area replay failed')
    return dict(clear=clear,minimum_gap_m=gap,placed_points=placed,area_relative_error=relative)
