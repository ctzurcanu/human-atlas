"""Conforming eight-child tetrahedral refinement with shared edge controls.

The central octahedron uses its shortest diagonal; every exterior face has
the same four-triangle subdivision on both neighboring tetrahedra.
"""
import numpy as np


def refine_tetrahedral_edges(positions, tetrahedra):
    p=np.asarray(positions,float);t=np.asarray(tetrahedra,dtype=np.int64)
    if p.ndim!=2 or p.shape[1]!=3 or t.ndim!=2 or t.shape[1]!=4 or not np.isfinite(p).all() or t.min()<0 or t.max()>=len(p):
        raise ValueError('Finite indexed tetrahedral source required')
    edge_corners=np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
    edges,inv=np.unique(np.sort(t[:,edge_corners].reshape(-1,2),axis=1),axis=0,return_inverse=True)
    midpoint=p[edges].mean(1);ids=inv.reshape(-1,6)+len(p);v=np.r_[p,midpoint]
    a,b,c,d=t.T;ab,ac,ad,bc,bd,cd=ids.T
    children=[np.c_[a,ab,ac,ad],np.c_[b,ab,bc,bd],np.c_[c,ac,bc,cd],np.c_[d,ad,bd,cd]]
    diagonals=np.stack([np.c_[ab,cd],np.c_[ac,bd],np.c_[ad,bc]],axis=1)
    chosen=np.argmin(np.linalg.norm(v[diagonals[:,:,0]]-v[diagonals[:,:,1]],axis=2),axis=1)
    cycles=[np.c_[ac,ad,bd,bc],np.c_[ab,ad,cd,bc],np.c_[ab,ac,cd,bd]]
    for k in range(3):
        selected=chosen==k;x,y=diagonals[selected,k].T;cycle=cycles[k][selected]
        for j in range(4):children.append(np.c_[x,y,cycle[:,j],cycle[:,(j+1)%4]])
    result=np.concatenate(children);q=v[result];det=np.linalg.det(np.transpose(q[:,1:]-q[:,:1],(0,2,1)))
    if np.any(det==0):raise ValueError('Source refinement produces a collapsed tetrahedron')
    neg=det<0;result[neg]=result[neg][:,[0,2,1,3]]
    return dict(positions=v,tetrahedra=result,edges=edges,midpoint_ids=np.arange(len(edges))+len(p))
