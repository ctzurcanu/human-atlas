"""Bounded halfedge collapses with every incident tissue's full link condition.

Topology condition: https://doc.cgal.org/latest/Surface_mesh_simplification/
No arbitrary close-point contractions. Immutable carrier controls are retained.
"""
from collections import defaultdict
import time
import numpy as np

def collapse_shared_stars(positions,cells,reference_normals,fixed,*,bound=2.5e-7,seconds=20.):
 start=time.monotonic();p=np.asarray(positions,float);f=[np.array(q,int,copy=True)for q in cells];live=[np.ones(len(q),bool)for q in f];vf=[]
 for faces in f:
  index=defaultdict(set)
  for i,t in enumerate(faces):
   for v in t:index[int(v)].add(i)
  vf.append(index)
 allfaces=np.concatenate(f);edges=np.unique(np.sort(np.concatenate([allfaces[:,[0,1]],allfaces[:,[1,2]],allfaces[:,[2,0]]]),axis=1),axis=0);length=np.linalg.norm(p[edges[:,0]]-p[edges[:,1]],axis=1);order=np.flatnonzero(length<=bound);order=order[np.argsort(length[order])];parent=np.arange(len(p));members={};moves=[];rejected=0
 def root(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return int(i)
 for edge in edges[order]:
  if time.monotonic()-start>=seconds:raise TimeoutError('Explicit shared-star conditioning deadline')
  a,b=map(root,edge)
  if a==b or fixed[a]and fixed[b]:continue
  if fixed[b]or not fixed[a]and b<a:a,b=b,a
  merged=members.get(a,[a])+members.get(b,[b])
  if np.linalg.norm(p[merged]-p[a],axis=1).max()>bound:continue
  allowed=True;updates=[]
  for k,faces in enumerate(f):
   fa=vf[k].get(a,set());fb=vf[k].get(b,set())
   if not fb:continue
   if fa:
    common=fa&fb
    if len(common)!=2:allowed=False;break
    def link(vertex,incident):
     pairs={tuple(sorted(int(v)for v in faces[i]if v!=vertex))for i in incident}
     return {v for pair in pairs for v in pair},pairs
    va,ea=link(a,fa);vb,eb=link(b,fb);edge_link={int(v)for i in common for v in faces[i]if v not in [a,b]}
    if va&vb!=edge_link or ea&eb:allowed=False;break
   ids=np.asarray(sorted(fb),int);old=faces[ids];new=old.copy();new[new==b]=a;keep=np.array([len(set(t))==3 for t in new]);q=p[old];on=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);q=p[new];nn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);reference=reference_normals[k][ids];od=np.einsum('ij,ij->i',on,reference);nd=np.einsum('ij,ij->i',nn,reference)
   if np.any(keep&(od>0)&(nd<=0)):allowed=False;break
   updates.append((k,ids,new,keep))
  if not allowed:rejected+=1;continue
  for k,ids,new,keep in updates:
   for i,t,active in zip(ids,new,keep):
    for v in f[k][i]:vf[k][int(v)].remove(int(i))
    f[k][i]=t;live[k][i]=active
    if active:
     for v in t:vf[k][int(v)].add(int(i))
  parent[b]=a;members[a]=merged;members.pop(b,None);moves.append(dict(source=b,target=a,maximum_displacement_m=float(np.linalg.norm(p[merged]-p[a],axis=1).max())))
 return [q[mask]for q,mask in zip(f,live)],dict(collapses=moves,rejected_link_or_orientation=rejected,maximum_allowed_displacement_m=bound,carrier_controls_literal=True,elapsed_seconds=time.monotonic()-start)
