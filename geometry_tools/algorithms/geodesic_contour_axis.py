"""Source-only tube axes from connected geodesic level contours.

The surface is never contracted, simplified or reconstructed. A declared mouth
seeds edge-geodesic distance; a surface path to the declared distal point picks
one connected contour at each level. Other branches are retained by the caller.
This is a placement guide, not an anatomical or lumen acceptance certificate.
"""
import time
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra, connected_components


def geodesic_contour_axis(vertices, faces, mouth_raw, distal_point, *,
                          spacing_m=.00075, wall_seconds=3., maximum_levels=512):
    started=time.monotonic()
    v=np.asarray(vertices,float);f=np.asarray(faces,int);mouth=np.asarray(mouth_raw,int)
    end=np.asarray(distal_point,float)
    if (v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or
        not np.isfinite(v).all() or end.shape!=(3,) or not np.isfinite(end).all() or
        not len(mouth) or mouth.min()<0 or mouth.max()>=len(v) or
        not 0<spacing_m<=.003 or not 0<wall_seconds<=5 or not 2<=maximum_levels<=2048):
        raise ValueError('Finite source mesh, actual mouth, distal point and bounded axis sampling required')
    points,physical=np.unique(v,axis=0,return_inverse=True);tri=physical[f]
    edges,edge_ids=np.unique(np.sort(tri[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0,return_inverse=True)
    edge_ids=edge_ids.reshape(-1,3)
    lengths=np.linalg.norm(points[edges[:,1]]-points[edges[:,0]],axis=1)
    valid=lengths>1e-15
    a,b=edges[valid].T
    graph=sparse.csr_matrix((np.r_[lengths[valid],lengths[valid]],(np.r_[a,b],np.r_[b,a])),shape=(len(points),len(points)))
    seeds=np.unique(physical[mouth])
    distance,predecessor,_=dijkstra(graph,indices=seeds,min_only=True,return_predecessors=True)
    if not np.isfinite(distance).all():raise ValueError('Axis source must be one connected physical component')
    tip=int(np.linalg.norm(points-end,axis=1).argmin());limit=float(distance[tip])
    if limit<2*spacing_m:raise ValueError('Declared distal point is too close to its mouth')
    path=[tip]
    while distance[path[-1]]>0:
        parent=int(predecessor[path[-1]])
        if parent<0 or len(path)>len(points):raise ValueError('Incomplete source geodesic predecessor chain')
        path.append(parent)
    path=np.array(path[::-1],int)
    lookup={tuple(edge):i for i,edge in enumerate(edges)}
    levels=np.arange(spacing_m*.5,limit-spacing_m*.1,spacing_m)
    if len(levels)>maximum_levels:raise ValueError('Source axis exceeds its explicit level budget')
    centres=[v[mouth].mean(0)];witnesses=[]
    for level in levels:
        if time.monotonic()-started>wall_seconds:raise TimeoutError('Native source-contour axis deadline')
        # Avoid a source vertex exactly on the contour, without moving any
        # source geometry or relying on coincident floating point segments.
        for _ in range(8):
            if np.min(abs(distance-level))>1e-10:break
            level+=3.71e-10
        else:raise ValueError('Unresolved critical geodesic contour level')
        crossing=(distance[edges[:,0]]<level)&(distance[edges[:,1]]>level)
        crossing|=(distance[edges[:,1]]<level)&(distance[edges[:,0]]>level)
        live=np.flatnonzero(crossing);local=np.full(len(edges),-1,int);local[live]=np.arange(len(live))
        in_face=crossing[edge_ids]
        if np.any((in_face.sum(1)!=0)&(in_face.sum(1)!=2)):
            raise ValueError('Nonregular source level intersection')
        face_ids=np.flatnonzero(in_face.sum(1)==2)
        segments=local[edge_ids[face_ids][in_face[face_ids]].reshape(-1,2)]
        adjacency=sparse.csr_matrix((np.ones(len(segments)*2),(segments.ravel(),segments[:,::-1].ravel())),shape=(len(live),len(live)))
        _,labels=connected_components(adjacency,directed=False)
        k=int(np.searchsorted(distance[path],level)-1)
        if not 0<=k<len(path)-1:raise ValueError('Declared distal path misses a contour level')
        path_edge=lookup[tuple(sorted((int(path[k]),int(path[k+1]))))]
        if local[path_edge]<0:raise ValueError('Geodesic path edge does not cross its contour')
        chosen=labels[local[path_edge]];keep=np.all(labels[segments]==chosen,axis=1)
        selected=segments[keep];ids=np.unique(selected)
        degree=np.bincount(selected.ravel(),minlength=len(live))
        if np.any(degree[ids]!=2):raise ValueError('Selected geodesic contour is open or branched')
        e=edges[live];fraction=(level-distance[e[:,0]])/(distance[e[:,1]]-distance[e[:,0]])
        hits=points[e[:,0]]+fraction[:,None]*(points[e[:,1]]-points[e[:,0]])
        spans=hits[selected];weight=np.linalg.norm(spans[:,1]-spans[:,0],axis=1)
        if weight.sum()<1e-12:raise ValueError('Collapsed source contour')
        centre=np.sum(spans.mean(1)*weight[:,None],axis=0)/weight.sum()
        centres.append(centre)
        witnesses.append(dict(level_m=float(level),sourceFaces=face_ids[keep].tolist(),
                              sourceEdges=live[ids].tolist(),perimeter_m=float(weight.sum()),
                              maximumContourRadius_m=float(np.linalg.norm(hits[ids]-centre,axis=1).max())))
    centres.append(end.copy());axis=np.asarray(centres)
    if np.any(np.linalg.norm(np.diff(axis,axis=0),axis=1)<1e-12):
        raise ValueError('Coincident source axis stations require explicit review')
    return axis,dict(method='Connected source edge-geodesic contours selected by a declared distal surface path',
        sourceVertices=len(v),sourceFaces=len(f),physicalVertices=len(points),
        sourceSurfacePathPhysicalIds=path.tolist(),distalGeodesicDistance_m=limit,
        contourWitnesses=witnesses,spacing_m=spacing_m,sourceGeometryModified=False,
        sourceBranchesRemoved=0,anatomicalAcceptance=False,fullTriangleReplayRequired=True,
        nativeWallSeconds=wall_seconds,elapsedSeconds=time.monotonic()-started)
