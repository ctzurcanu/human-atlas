"""Bounded clearance-weighted routing around fixed solids inside an oriented skin.

The grid and route are proposals. Full transported triangles still require
independent contact, source-interface and rendered review.
"""
import time
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra
from scipy.spatial import cKDTree
from geometry_tools.open3d_backend import library


def sampled_clearance(points,bone,skin,surface_obstacles=()):
    q=np.asarray(points,float);bh,_=bone.closest_points(q);sh,_=skin.closest_points(q)
    inside=bone.scene.compute_occupancy(library().core.Tensor(((q-bone.origin)/bone.scale).astype(np.float32)),
                                      nthreads=2,nsamples=3).numpy()>.5
    t=skin.mesh.triangles[sh['face_ids']]
    n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    sg=np.sum((sh['points']-q)*n,axis=1)
    gap=np.minimum(bh['distance_m']*np.where(inside,-1.,1.),sg)
    for index in surface_obstacles:
        hit,_=index.closest_points(q);gap=np.minimum(gap,hit['distance_m'])
    return gap


def free_space_route(start,end,bone,skin,guide,*,radius_m,clearance_m=.0002,
                     grid_m=.001,padding_m=.015,maximum_grid_points=500000,wall_seconds=16.,
                     surface_obstacles=()):
    begin=time.monotonic();start,end=np.asarray(start,float),np.asarray(end,float)
    guide=np.asarray(guide,float)
    if not 0<grid_m<=.002 or not 0<radius_m<.005 or not 0<wall_seconds<=20:
        raise ValueError('Finite grid, source radius and native deadline required')
    def deadline():
        if time.monotonic()-begin>wall_seconds:raise TimeoutError('Native free-space route deadline')
    anchors=np.r_[guide,start[None],end[None]]
    low=anchors.min(0)-padding_m;high=anchors.max(0)+padding_m
    shape=np.ceil((high-low)/grid_m).astype(int)+1
    total=int(np.prod(shape))
    if total>maximum_grid_points:raise ValueError('Declared route grid exceeds its finite point budget')
    lattice=np.indices(shape).reshape(3,-1).T;points=low+lattice*grid_m
    gap=np.empty(total)
    for first in range(0,total,32768):
        gap[first:first+32768]=sampled_clearance(points[first:first+32768],bone,skin,surface_obstacles);deadline()
    margin=radius_m+clearance_m+grid_m*.5
    free=gap>=margin
    ids=np.full(total,-1,int);ids[free]=np.arange(int(free.sum()))
    if not free.any():raise ValueError('No sampled tube corridor inside this fixed context')
    live=points[free];tree=cKDTree(live)
    _,ends=tree.query([start,end]);ends=np.asarray(ends,int)
    snap=np.linalg.norm(live[ends]-[start,end],axis=1)
    if (snap>grid_m*2).any():raise ValueError('A fixed mouth lacks a nearby free corridor state')
    distance=cKDTree(guide).query(live)[0]
    cost=1.+.3*(distance/.005)**2+(.002/np.maximum(gap[free],.0001))**2
    cube=ids.reshape(shape);sources=[];destinations=[]
    for axis in range(3):
        left=[slice(None)]*3;right=left.copy();left[axis]=slice(None,-1);right[axis]=slice(1,None)
        a,b=cube[tuple(left)].ravel(),cube[tuple(right)].ravel();valid=(a>=0)&(b>=0)
        sources.append(a[valid]);destinations.append(b[valid])
    a,b=np.concatenate(sources),np.concatenate(destinations)
    weights=grid_m*(cost[a]+cost[b])*.5
    graph=sparse.csr_matrix((np.r_[weights,weights],(np.r_[a,b],np.r_[b,a])),shape=(len(live),len(live)))
    distance,parent=dijkstra(graph,directed=False,indices=int(ends[0]),return_predecessors=True)
    deadline()
    if not np.isfinite(distance[ends[1]]):raise ValueError('Fixed mouths are in different sampled free-space components')
    k=int(ends[1]);path=[k]
    while k!=ends[0]:
        k=int(parent[k])
        if k<0:raise ValueError('Incomplete route predecessor chain')
        path.append(k)
    route=np.r_[start[None],live[path[::-1]],end[None]]
    # Pull straight spans only when their actual context samples are clear.
    pulled=[route[0]];i=0;queries=0
    while i<len(route)-1:
        candidates=[];samples=[];ends_at=[]
        for j in range(i+1,len(route)):
            count=max(3,int(np.ceil(np.linalg.norm(route[j]-route[i])/(grid_m*.25)))+1)
            samples.append(np.linspace(route[i],route[j],count));candidates.append(j);ends_at.append(count)
        values=sampled_clearance(np.concatenate(samples),bone,skin,surface_obstacles);queries+=len(values)
        chunks=np.split(values,np.cumsum(ends_at)[:-1])
        good=[j for j,c in zip(candidates,chunks)if c.min()>=radius_m+clearance_m+grid_m*.25]
        if not good:raise ValueError('Grid route does not replay as a finite-radius span')
        i=good[-1];pulled.append(route[i]);deadline()
    result=np.asarray(pulled)
    return result,dict(gridPoints=total,freeGridPoints=int(free.sum()),gridSpacing_m=grid_m,
        sourceRadius_m=radius_m,clearance_m=clearance_m,gridMargin_m=margin,
        snappedMouthDistances_m=snap.tolist(),rawPathPoints=len(route),pulledPathPoints=len(result),
        replaySamples=queries,skinIsOrientedOpenSurface=True,fullTriangleAuditRequired=True,
        additionalOpenSurfaceObstacles=len(surface_obstacles),
        elapsedSeconds=time.monotonic()-begin,nativeWallSeconds=wall_seconds)
