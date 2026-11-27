"""Interval derivatives with a declared rotating plane, not a frozen world plane.

These are local proposals only. A finite pose requires original-pair contact
replay and new linearization; a derivative is never clearance evidence.
"""
import numpy as np
from scipy import sparse

def corotating_vertex_pair_rows(first_points,second_points,first_groups,second_groups,
        centers,levers,axes,plane_groups,plane_weights):
    a,b,n=[np.asarray(v,float) for v in [first_points,second_points,axes]]
    ga,gb=[np.asarray(v) for v in [first_groups,second_groups]]
    c,l=np.asarray(centers,float),np.asarray(levers,float)
    pg,pw=np.asarray(plane_groups),np.asarray(plane_weights,float)
    if a.ndim!=2 or a.shape[1]!=3 or b.shape!=a.shape or n.shape!=a.shape or ga.shape!=(len(a),) or gb.shape!=ga.shape or c.shape!=(len(l),3) or pg.shape!=pw.shape or pg.ndim!=2 or len(pg)!=len(a):
        raise ValueError('Explicit matching source vertex pairs, pose cells and plane ownership required')
    if any(g.dtype.kind not in 'iu' or np.any(g < -1) or np.any(g>=len(c)) for g in [ga,gb,pg]) or np.any(l<=0) or np.any(pw<0) or np.any(pw.sum(1)>1+1e-12) or not all(np.isfinite(v).all() for v in [a,b,n,c,l,pw]):
        raise ValueError('Finite proper pose cells, nonnegative plane weights and fixed-cell sentinel required')
    if np.any(abs(np.linalg.norm(n,axis=1)-1)>1e-10):raise ValueError('Unit separating-plane axes required')
    count=len(c);rows=[];cols=[];values=[];row_ids=np.arange(len(a))
    def add(groups,offset,vectors,mask):
        rows.append(np.repeat(row_ids[mask],3))
        cols.append((3*(offset+groups[mask,None])+np.arange(3)).ravel())
        values.append(vectors[mask].ravel())
    for sign,p,g in [(1.,a,ga),(-1.,b,gb)]:
        live=g>=0;safe=np.maximum(g,0)
        add(g,0,sign*n,live)
        add(g,count,sign*np.cross(p-c[safe],n)/l[safe,None],live)
    rotation=np.cross(n,a-b)
    for k in range(pg.shape[1]):
        g,w=pg[:,k],pw[:,k];safe=np.maximum(g,0)
        add(g,count,w[:,None]*rotation/l[safe,None],(g>=0)&(w>0))
    result=sparse.coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),shape=(len(a),6*count)).tocsr()
    result.eliminate_zeros()
    return result,np.einsum('ij,ij->i',a-b,n)
