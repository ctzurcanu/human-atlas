"""Bound control transitions by literal source-mesh graph distance.

This supplies proposal weights, not a clearance or injectivity certificate.
Euclidean proximity never creates a graph edge or an anatomical dependency.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

def fixed_interface_weights(vertices, faces, cells, fixed_cells, translations,
                            *, transition_m=.020, displacement_multiple=4.):
 v=np.asarray(vertices,float);f=np.asarray(faces,int);cell=np.asarray(cells,int);delta=np.asarray(translations,float)
 if v.ndim!=2 or v.shape[1]!=3 or cell.shape!=(len(v),)or delta.ndim!=2 or delta.shape[1]!=3 or not np.isfinite(v).all()or not np.isfinite(delta).all():raise ValueError('Finite source positions, translation cells and control vectors required')
 n=len(delta);fixed=np.array(sorted(set(map(int,fixed_cells))),int)
 if not len(fixed)or np.any(fixed<0)or np.any(fixed>=n)or np.any(cell<0)or np.any(cell>=n):raise ValueError('Explicit valid fixed source cells required')
 if not np.isfinite(transition_m)or transition_m<=0 or not np.isfinite(displacement_multiple)or displacement_multiple<1:raise ValueError('Positive transition length and displacement multiple >=1 required')
 edges=f[:,[0,1,1,2,2,0]].reshape(-1,2);length=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1);pairs=np.sort(cell[edges],axis=1);keep=(pairs[:,0]!=pairs[:,1])&(length>0);pairs,length=pairs[keep],length[keep]
 unique,inverse=np.unique(pairs,axis=0,return_inverse=True);distance=np.full(len(unique),np.inf);np.minimum.at(distance,inverse,length)
 graph=coo_matrix((np.r_[distance,distance],(np.r_[unique[:,0],unique[:,1]],np.r_[unique[:,1],unique[:,0]])),shape=(n,n)).tocsr()
 distances=dijkstra(graph,directed=False,indices=fixed,min_only=True)
 required=np.maximum(transition_m,displacement_multiple*np.linalg.norm(delta,axis=1));x=np.clip(distances/required,0,1);weight=x*x*x*(10+x*(-15+6*x));weight[fixed]=0
 return weight,dict(distance_m=distances,required_transition_m=required,source_graph_edges=len(unique),explicit_fixed_cells=len(fixed))
