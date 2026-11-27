"""Shortest source-edge handle cycles around an explicitly declared conduit.

An angular winding alone can enclose ordinary solid tissue. The face-dual
cut must also remain connected: contractible surface loops are rejected.
The result is a geometric handle reference, not anatomical acceptance.
"""
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra,connected_components


def handle_cycle(vertices,faces,center,axis=(0.,1.,0.),*,half_length=.025,radius=.04,
                 excluded_centers=(),maximum_starts=24,wall_seconds=8.):
 start=time.monotonic();deadline=start+wall_seconds
 if wall_seconds<=0 or maximum_starts<1 or half_length<=0 or radius<=0:raise ValueError('Positive finite cycle-search limits required')
 v=np.asarray(vertices,float);f=np.asarray(faces,int);c=np.asarray(center,float);a=np.asarray(axis,float)
 if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or c.shape!=(3,)or a.shape!=(3,)or not np.isfinite(np.r_[v.ravel(),c,a]).all():raise ValueError('Finite indexed geometry, center and axis required')
 a=a/np.linalg.norm(a);seed=np.eye(3)[np.argmin(abs(a))];u=seed-a*(seed@a);u/=np.linalg.norm(u);w=np.cross(a,u)
 q=v-c;axial=q@a;r=np.linalg.norm(q-axial[:,None]*a,axis=1);angle=np.arctan2(q@w,q@u)
 allowed=(abs(axial)<=half_length)&(r<=radius)&(r>1e-10);edges,owner=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0,return_inverse=True)
 incident=[[]for _ in edges]
 for i,j in enumerate(owner):incident[j].append(i%len(f))
 if any(len(x)!=2 for x in incident):raise ValueError('A closed exact-index two-manifold is required')
 length=np.linalg.norm(v[edges[:,0]]-v[edges[:,1]],axis=1);edge_id={tuple(e):i for i,e in enumerate(edges)}
 take=allowed[edges].all(1);selected=edges[take];cost=length[take];n=len(v);rr=[];cc=[];ww=[]
 for direction in [0,1]:
  first,second=selected[:,direction],selected[:,1-direction];raw=angle[second]-angle[first];wrapped=(raw+np.pi)%(2*np.pi)-np.pi;shift=np.rint((angle[first]+wrapped-angle[second])/(2*np.pi)).astype(int)
  for sheet in range(3):
   target=sheet+shift;ok=(target>=0)&(target<3);rr.extend((sheet*n+first[ok]).tolist());cc.extend((target[ok]*n+second[ok]).tolist());ww.extend(cost[ok].tolist())
 graph=coo_matrix((ww,(rr,cc)),shape=(3*n,3*n)).tocsr()
 candidates=[]
 for left,right in zip(np.linspace(-np.pi,np.pi,maximum_starts+1)[:-1],np.linspace(-np.pi,np.pi,maximum_starts+1)[1:]):
  ids=np.flatnonzero(allowed&(angle>=left)&(angle<right))
  if len(ids):candidates.append(int(ids[np.argmin(r[ids]+.2*abs(axial[ids]))]))
 best=None;rejected=[]
 for anchor in dict.fromkeys(candidates):
  if time.monotonic()>=deadline:raise TimeoutError('Source handle-cycle search reached its finite deadline')
  distance,pred=dijkstra(graph,directed=True,indices=n+anchor,return_predecessors=True,limit=2*np.pi*radius*2)
  target=2*n+anchor
  if not np.isfinite(distance[target]):continue
  path=[target]
  while path[-1]!=n+anchor:
   parent=int(pred[path[-1]]);assert parent>=0;path.append(parent)
  cycle=np.array(path[::-1],int)%n;cycle=cycle[:-1]
  if len(cycle)<3 or len(np.unique(cycle))!=len(cycle):continue
  def winding(point):
   rel=v[cycle]-point;theta=np.arctan2(rel@w,rel@u);delta=(np.roll(theta,-1)-theta+np.pi)%(2*np.pi)-np.pi;return float(delta.sum()/(2*np.pi))
  if abs(winding(c)-1)>1e-8 or any(abs(winding(np.asarray(point,float)))>.5 for point in excluded_centers):continue
  blocked={edge_id[tuple(sorted((int(i),int(j))))]for i,j in zip(cycle,np.roll(cycle,-1))};face_edges=np.array([x for i,x in enumerate(incident)if i not in blocked],int)
  dual=coo_matrix((np.ones(2*len(face_edges)),(np.r_[face_edges[:,0],face_edges[:,1]],np.r_[face_edges[:,1],face_edges[:,0]])),shape=(len(f),len(f))).tocsr()
  components=connected_components(dual,directed=False,return_labels=False)
  if components!=1:rejected.append(dict(reason='separating or contractible cycle',vertices=len(cycle)));continue
  value=float(distance[target])
  if best is None or value<best[0]:best=(value,cycle)
 if best is None:raise ValueError('No nonseparating source-edge handle cycle satisfies the declared conduit constraints')
 return best[1],dict(perimeter_m=best[0],vertices=len(best[1]),winding=1,face_dual_components_after_cut=1,nonseparating_source_cycle=True,excluded_axes=len(excluded_centers),tested_starts=len(candidates),rejected_cycles=rejected,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False)
