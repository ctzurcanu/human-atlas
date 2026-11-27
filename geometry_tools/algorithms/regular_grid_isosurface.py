"""Conforming marching tetrahedra on an explicit bounded regular SDF grid.

Negative values denote occupied tissue. Every cut edge is evaluated once and
shared by its incident tetrahedra. This adds no physiological opening caps.
Geometry/topology/source fidelity still require independent audits.
"""
import numpy as np

def extract(field,origin,spacing,*,maximum_active_cells=1200000):
    d=np.asarray(field);shape=np.asarray(d.shape,int);origin=np.asarray(origin,float);h=float(spacing)
    if d.ndim!=3 or min(shape)<2 or not np.isfinite(d).all()or origin.shape!=(3,)or h<=0:raise ValueError('Finite 3D regular field, origin and positive spacing required')
    roundoff=128*np.finfo(float).eps*max(h,float(np.max(np.abs(d))))
    near_zero=(d!=0)&(np.abs(d)<roundoff)
    if np.any(near_zero):d=d.copy();d[near_zero]=0
    offsets=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]])
    lo=np.full(tuple(shape-1),np.inf,dtype='<f4');hi=-lo.copy()
    for off in offsets:
        q=d[tuple(slice(int(k),int(k+s-1))for k,s in zip(off,shape))];np.minimum(lo,q,out=lo);np.maximum(hi,q,out=hi)
    cells=np.column_stack(np.nonzero((lo<0)&(hi>=0)));del lo,hi
    if len(cells)>maximum_active_cells:raise ValueError('Isosurface active-cell budget reached')
    strides=np.array([shape[1]*shape[2],shape[2],1],np.int64);corner_ids=(cells@strides)[:,None]+(offsets@strides)[None,:]
    tetrahedra=np.array([[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]])
    flat=d.reshape(-1);tri_edges=[];inside_points=[];outside_points=[]
    for tet in tetrahedra:
        ids=corner_ids[:,tet];values=flat[ids];mask=np.sum((values<0)*(1<<np.arange(4)),axis=1)
        for case in range(1,15):
            rows=np.flatnonzero(mask==case)
            if not len(rows):continue
            inner=[i for i in range(4)if case&(1<<i)];outer=[i for i in range(4)if i not in inner];source=ids[rows]
            if len(inner)==1:patterns=[[(inner[0],j)for j in outer]]
            elif len(inner)==3:patterns=[[(j,outer[0])for j in inner]]
            else:
                a,b=inner;c,e=outer;patterns=[[(a,c),(b,c),(b,e)],[(a,c),(b,e),(a,e)]]
            for pattern in patterns:
                tri_edges.append(np.sort(source[:,np.asarray(pattern)],axis=2))
                inside_points.append(source[:,inner[0]]);outside_points.append(source[:,outer[0]])
    edges=np.concatenate(tri_edges);pairs,inverse=np.unique(edges.reshape(-1,2),axis=0,return_inverse=True);faces=inverse.reshape(-1,3)
    ijk=np.column_stack(np.unravel_index(pairs.ravel(),tuple(shape))).reshape(-1,2,3);values=flat[pairs].astype(float)
    # A zero grid value needs one consistent vertex key instead of several
    # coincident keys from incident edges.
    t=values[:,0]/(values[:,0]-values[:,1]);vertices=origin+h*(ijk[:,0]+t[:,None]*(ijk[:,1]-ijk[:,0]))
    _,first,canonical=np.unique(vertices,axis=0,return_index=True,return_inverse=True);vertices=vertices[first];faces=canonical[faces]
    inside=np.column_stack(np.unravel_index(np.concatenate(inside_points),tuple(shape)));outside=np.column_stack(np.unravel_index(np.concatenate(outside_points),tuple(shape)))
    q=vertices[faces];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);flip=np.einsum('ij,ij->i',n,outside-inside)<0;faces[flip]=faces[flip][:,[0,2,1]]
    keep=(faces[:,0]!=faces[:,1])&(faces[:,1]!=faces[:,2])&(faces[:,2]!=faces[:,0]);faces=faces[keep]
    return vertices,faces,dict(active_cells=len(cells),cut_edges=len(pairs),vertices=len(vertices),triangles=len(faces),coincident_zero_controls_welded=True)
