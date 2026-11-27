"""One diagonal for both oriented uses of a corresponding wall edge."""
import numpy as np


def canonical_wall_rims(oriented_edges, inner_offset):
 edges=np.asarray(oriented_edges,dtype=np.int64)
 if edges.ndim!=2 or edges.shape[1]!=2 or not isinstance(inner_offset,(int,np.integer))or inner_offset<=0 or (len(edges)and(edges.min()<0 or edges.max()>=inner_offset or np.any(edges[:,0]==edges[:,1]))):
  raise ValueError('Distinct oriented carrier edges and a positive inner-bank offset required')
 a,b=edges.T;first=np.c_[b,a,a+inner_offset];second=np.c_[b,a+inner_offset,b+inner_offset];reverse=a>b;first[reverse]=np.c_[b[reverse],a[reverse],b[reverse]+inner_offset];second[reverse]=np.c_[a[reverse],a[reverse]+inner_offset,b[reverse]+inner_offset]
 return np.r_[first,second]
