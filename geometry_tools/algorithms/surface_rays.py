"""Bounded Float64 first-hit rays on original triangles through the mesh BVH.

This is a geometric measurement, not an anatomical wall-thickness inference.
No normalized Float32 copy is substituted for the original source surface.
"""
import time
import numpy as np

from ..spatial import BVH
from ..distance import barycentric


def first_surface_hits(mesh,origins,directions,*,index=None,minimum_m=0.,maximum_m=np.inf,
                       time_limit=5.,max_triangle_tests=1000000):
    origins=np.asarray(origins,float);directions=np.asarray(directions,float)
    if (origins.ndim!=2 or origins.shape[1:]!=(3,)or directions.shape!=origins.shape
        or not np.isfinite(origins).all()or not np.isfinite(directions).all()
        or minimum_m<0 or maximum_m<=minimum_m or time_limit<=0 or max_triangle_tests<1):
        raise ValueError('Finite ray arrays and positive, bounded query controls are required')
    lengths=np.linalg.norm(directions,axis=1)
    if np.any(lengths==0):raise ValueError('A surface ray requires a nonzero direction')
    directions=directions/lengths[:,None]
    started=time.monotonic();deadline=started+time_limit
    index=BVH(mesh)if index is None else index
    if index.mesh.geometry_sha!=mesh.geometry_sha or index.mesh.frame!=mesh.frame:
        raise ValueError('Ray BVH does not bind the original mesh and frame')
    distances=np.full(len(origins),np.inf);face_ids=np.full(len(origins),-1,np.int64)
    weights=np.full((len(origins),3),np.nan);points=np.full((len(origins),3),np.nan)
    tests=0;eps=np.finfo(float).eps
    scale=max(1.,float(np.abs(mesh.vertices).max()))
    broad_bound=64*eps*scale

    def entry(origin,direction,node,limit):
        low=index.low[node]-broad_bound;high=index.high[node]+broad_bound
        near,far=minimum_m,limit
        for axis in range(3):
            if direction[axis]==0.:
                if origin[axis]<low[axis]or origin[axis]>high[axis]:return np.inf
            else:
                a=(low[axis]-origin[axis])/direction[axis];b=(high[axis]-origin[axis])/direction[axis]
                near=max(near,min(a,b));far=min(far,max(a,b))
                if far<near:return np.inf
        return near

    for ray,(origin,direction)in enumerate(zip(origins,directions)):
        best=maximum_m;stack=[(entry(origin,direction,0,best),0)]
        while stack:
            if time.monotonic()>=deadline:raise TimeoutError('Original-triangle ray query reached its declared deadline')
            distance,node=stack.pop()
            if not np.isfinite(distance)or distance>best:continue
            if index.leaf(node):
                ids=index.ids(node);tests+=len(ids)
                if tests>max_triangle_tests:raise RuntimeError('Original-triangle ray query exhausted its declared triangle budget')
                tri=mesh.triangles[ids]
                normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
                denominator=normal@direction
                parallel=np.abs(denominator)<=64*eps*np.linalg.norm(normal,axis=1)
                candidates=np.flatnonzero(~parallel)
                if not len(candidates):continue
                t=np.einsum('ij,ij->i',normal[candidates],tri[candidates,0]-origin)/denominator[candidates]
                keep=(t>minimum_m)&(t<=best)
                candidates,t=candidates[keep],t[keep]
                if not len(candidates):continue
                actual=origin+t[:,None]*direction
                w,nondegenerate=barycentric(actual,tri[candidates])
                keep=nondegenerate&np.isfinite(w).all(1)&(w.min(1)>=-128*eps)&(w.max(1)<=1+128*eps)
                if not keep.any():continue
                candidates,t,w,actual=candidates[keep],t[keep],w[keep],actual[keep]
                nearest=int(np.argmin(t));best=float(t[nearest])
                distances[ray]=best;face_ids[ray]=int(ids[candidates[nearest]])
                weights[ray]=w[nearest];points[ray]=actual[nearest]
            else:
                children=[(entry(origin,direction,int(child),best),int(child))for child in index.children[node]]
                children.sort(reverse=True)
                stack.extend(children)
    return dict(distance_m=distances,face_ids=face_ids,barycentric=weights,points=points,
                hit=face_ids>=0),dict(backend='float64-bvh-original-triangles',frame=mesh.frame,
                geometry_sha256=mesh.geometry_sha,original_coordinates_preserved=True,
                queries=len(origins),tested_triangles=tests,complete_within_budget=True,
                maximum_m=float(maximum_m),minimum_m=float(minimum_m),elapsed_seconds=time.monotonic()-started,
                wall_thickness_or_anatomy_inferred=False)
