"""Bounded axis-parallel rays against literal Float64 source triangles.

Useful for clearance proposals when a Float32 accelerator cannot preserve a
thin mesh. These are numerical witnesses, not exact predicates or a continuous
surface-clearance certificate. Shared-edge/tangent ambiguities stay explicit.
"""
from __future__ import annotations
import numpy as np
from ..spatial import BVH


def cast_axis_rays(mesh, origins, *, axis=1, sign=-1, bvh=None,
                   barycentric_epsilon=1e-12, max_triangle_tests=2000000):
    origins=np.asarray(origins,dtype=np.float64)
    if origins.ndim!=2 or origins.shape[1]!=3 or not np.isfinite(origins).all():
        raise ValueError('Finite origins[N,3] are required')
    if axis not in (0,1,2) or sign not in (-1,1):
        raise ValueError('An axis 0/1/2 and direction sign -1/+1 are required')
    if not 0<=barycentric_epsilon<=1e-8 or max_triangle_tests<1:
        raise ValueError('Invalid numerical margin or work budget')
    bvh=bvh or BVH(mesh)
    if bvh.mesh.geometry_sha!=mesh.geometry_sha:
        raise ValueError('BVH belongs to another source mesh')
    count=len(origins);face_ids=np.full(count,-1,np.int64)
    weights=np.full((count,3),np.nan);points=np.full((count,3),np.nan)
    distances=np.full(count,np.nan);ambiguous=np.zeros(count,bool)
    projected_tangent_candidates=np.zeros(count,np.int64)
    transverse=[k for k in range(3) if k!=axis];tested=0
    endpoint=mesh.bounds[1 if sign>0 else 0,axis]
    for i,origin in enumerate(origins):
        if sign*(endpoint-origin[axis])<0:continue
        low=origin.copy();high=origin.copy()
        low[axis]=min(origin[axis],endpoint);high[axis]=max(origin[axis],endpoint)
        ids=bvh.query_bounds(low,high)
        tested+=len(ids)
        if tested>max_triangle_tests:raise RuntimeError('Axis-ray triangle-test budget exceeded')
        if not len(ids):continue
        triangles=mesh.triangles[ids];a=triangles[:,0]
        u=triangles[:,1]-a;v=triangles[:,2]-a;d=origin-a
        j,k=transverse;det=u[:,j]*v[:,k]-u[:,k]*v[:,j]
        valid=det!=0
        projected_tangent_candidates[i]=int((~valid).sum())
        ids=ids[valid];a=a[valid];u=u[valid];v=v[valid];d=d[valid];det=det[valid]
        s=(d[:,j]*v[:,k]-d[:,k]*v[:,j])/det
        w=(u[:,j]*d[:,k]-u[:,k]*d[:,j])/det
        bary=np.c_[1-s-w,s,w]
        surface=a+s[:,None]*u+w[:,None]*v
        travel=sign*(surface[:,axis]-origin[axis])
        hit=(bary.min(1)>=-barycentric_epsilon)&(travel>=0)
        if not hit.any():continue
        found=np.flatnonzero(hit)
        nearest=found[np.lexsort((ids[found],travel[found]))[0]]
        face_ids[i]=ids[nearest];weights[i]=bary[nearest]
        points[i]=surface[nearest];distances[i]=travel[nearest]
        ambiguous[i]=bool(bary[nearest].min()<=barycentric_epsilon or projected_tangent_candidates[i])
    hit=face_ids>=0
    return dict(hit=hit,face_ids=face_ids,barycentric=weights,points=points,
                distance_m=distances,boundary_ambiguous=ambiguous,
                projected_tangent_candidates=projected_tangent_candidates),dict(
        backend='literal-float64-axis-ray-bvh',geometry_sha256=mesh.geometry_sha,
        queries=count,hits=int(hit.sum()),tested_triangles=tested,
        max_triangle_tests=max_triangle_tests,barycentric_epsilon=barycentric_epsilon,
        source_positions_preserved=True,approximate_query=False,
        precision='Float64 arithmetic, not exact predicates',anatomical_acceptance=False)
