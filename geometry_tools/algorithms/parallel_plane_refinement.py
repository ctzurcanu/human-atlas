"""Conforming subdivision by parallel planes, without changing the surface.

Each added point belongs to one original edge. Adjacent faces reuse its exact
edge/plane identity; no tolerance welding, repair, anatomical cap or smoothing.
Original controls retain their indices. Face and point ancestry permit replay.
"""
import time
import numpy as np


def refine_parallel_planes(vertices, faces, levels, *, axis=1,
                           time_limit_seconds=10., maximum_faces=200000):
    start=time.monotonic()
    v=np.asarray(vertices,float);f=np.asarray(faces,int);levels=np.asarray(levels,float)
    if (v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3
            or not np.isfinite(v).all() or not np.isfinite(levels).all()
            or levels.ndim!=1 or len(np.unique(levels))!=len(levels)
            or axis not in (0,1,2) or not 0<time_limit_seconds<=20
            or np.any(f<0) or np.any(f>=len(v))):
        raise ValueError('Finite source triangles and distinct finite planes required')
    levels=np.sort(levels);points=list(v.copy());edges=[[i,i]for i in range(len(v))]
    weights=[[1.,0.]for _ in v];plane_ids=[-1]*len(v);lookup={};triangles=[];parents=[]

    def node(bary,tri):
        ids=np.flatnonzero(bary)
        if len(ids)==1:return int(tri[ids[0]])
        if len(ids)!=2:raise ValueError('Refinement point left its original source edge')
        a,b=sorted(map(int,tri[ids]),key=lambda k:(tuple(v[k]),k));y=float(bary@v[tri,axis])
        li=int(np.argmin(abs(levels-y)))
        if abs(levels[li]-y)>1e-12:raise ValueError('Missing source plane identity')
        key=(a,b,li)
        if key not in lookup:
            t=float((levels[li]-v[a,axis])/(v[b,axis]-v[a,axis]))
            if not 0<t<1:raise ValueError('Added source-edge point is not strictly internal')
            p=(1-t)*v[a]+t*v[b];p[axis]=levels[li]
            lookup[key]=len(points);points.append(p);edges.append([a,b]);weights.append([1-t,t]);plane_ids.append(li)
        return lookup[key]

    def clip(poly,tri,y,keep_above):
        result=[]
        for a,b in zip(poly,poly[1:]+poly[:1]):
            sa=float(a@v[tri,axis]-y);sb=float(b@v[tri,axis]-y)
            ia=sa>=0 if keep_above else sa<=0;ib=sb>=0 if keep_above else sb<=0
            if ia:result.append(a)
            if ia!=ib:
                t=sa/(sa-sb);result.append((1-t)*a+t*b)
        return result

    for fi,tri in enumerate(f):
        lo,hi=float(v[tri,axis].min()),float(v[tri,axis].max())
        cuts=levels[(levels>lo)&(levels<hi)]
        if not len(cuts):triangles.append(tri.tolist());parents.append(fi);continue
        bounds=np.r_[lo,cuts,hi]
        for lower,upper in zip(bounds,bounds[1:]):
            poly=clip(clip(list(np.eye(3)),tri,lower,True),tri,upper,False)
            ids=[node(p,tri)for p in poly]
            ids=[p for j,p in enumerate(ids)if j==0 or p!=ids[j-1]]
            if len(ids)>1 and ids[-1]==ids[0]:ids.pop()
            for j in range(1,len(ids)-1):
                face=[ids[0],ids[j],ids[j+1]]
                if len(set(face))==3:triangles.append(face);parents.append(fi)
        if len(triangles)>maximum_faces:raise ValueError('Explicit refinement face budget exceeded')
        if fi%64==0 and time.monotonic()-start>time_limit_seconds:
            raise TimeoutError('Native parallel-plane refinement deadline')
    result=dict(vertices=np.asarray(points),faces=np.asarray(triangles,int),
                parent_faces=np.asarray(parents,int),vertex_source_edges=np.asarray(edges,int),
                vertex_source_weights=np.asarray(weights),vertex_plane_ids=np.asarray(plane_ids,int),
                levels=levels,axis=axis,original_vertex_count=len(v))
    q=result['vertices'][result['faces']];cross=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
    p=v[f];original_cross=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
    if np.any(np.einsum('ij,ij->i',cross,original_cross[result['parent_faces']])<=0):
        raise ValueError('Degenerate source/subtriangle or reversed refinement orientation')
    area=np.bincount(result['parent_faces'],weights=np.linalg.norm(cross,axis=1),minlength=len(f))
    if not np.allclose(area,np.linalg.norm(original_cross,axis=1),rtol=1e-10,atol=1e-16):
        raise ValueError('Source face partition failed area replay')
    result['elapsed_seconds']=time.monotonic()-start
    return result
