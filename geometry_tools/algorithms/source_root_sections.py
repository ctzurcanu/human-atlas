"""Recover source neck sections from literal mouth-to-mouth annular incidence.

Unlike exporter-quad proposals, this handles a source neck changing from twelve
to six controls. Recovery stops at the first branch or non-annular layer.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra,connected_components
from .source_tube_topology import _annulus


def recover_root_sections(vertices,faces,root_raw_ids,*,maximum_layers=64):
 p=np.asarray(vertices,float);f=np.asarray(faces,int);root=np.asarray(root_raw_ids,int)
 if not isinstance(maximum_layers,int)or not 1<=maximum_layers<=256:raise ValueError('Finite source layer budget required')
 xyz,representative,inv=np.unique(p,axis=0,return_index=True,return_inverse=True);ff=inv[f];e=np.unique(np.sort(ff[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0)
 graph=coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(xyz),len(xyz))).tocsr();distance=dijkstra(graph,directed=False,indices=np.unique(inv[root]),unweighted=True,min_only=True)
 previous=np.unique(inv[root]);rows=[]
 for level in range(1,maximum_layers+1):
  ring=np.flatnonzero(distance==level)
  if len(ring)<3:break
  induced=graph[ring][:,ring]
  if not np.all(np.diff(induced.indptr)==2)or connected_components(induced,directed=False)[0]!=1:break
  order=[int(ring[0])]
  while len(order)<len(ring):
   neighbors=graph.indices[graph.indptr[order[-1]]:graph.indptr[order[-1]+1]];choices=sorted(set(map(int,neighbors))&set(map(int,ring))-set(order))
   if not choices:raise ValueError('A source neck layer does not close into one simple cycle')
   order.append(choices[0])
  if order[0]not in graph.indices[graph.indptr[order[-1]]:graph.indptr[order[-1]+1]]:raise ValueError('Unclosed source neck ring')
  band=np.flatnonzero(np.all((distance[ff]>=level-1)&(distance[ff]<=level),axis=1))
  if not len(band)or not _annulus(ff[band],previous,ring):break
  rows.append(dict(layer=level,rawIds=representative[np.array(order)].tolist(),sourceFaceIds=band.tolist(),certifiedManifoldAnnulus=True));previous=ring
 return rows


def source_annular_stem(vertices,faces,root_raw_ids,profile_candidates):
 """Find the longest branch-free source sleeve ending on a declared profile.

 Face adjacency is cut only at the actual candidate ring edges. Every accepted
 sleeve must have exactly its two declared source boundaries and Euler zero.
 This avoids assuming that triangulation graph-distance layers are sections.
 """
 from collections import defaultdict
 from .tube_continuity import boundary_loops,continuity_audit
 p=np.asarray(vertices,float);f=np.asarray(faces,int);physical,representative,inv=np.unique(p,axis=0,return_index=True,return_inverse=True);ff=inv[f];root=set(map(int,inv[np.asarray(root_raw_ids,int)]));edge_faces=defaultdict(list)
 for i,tri in enumerate(ff):
  for a,b in zip(tri,np.roll(tri,-1)):edge_faces[tuple(sorted((int(a),int(b))))].append(i)
 if any(len(ids)>2 for ids in edge_faces.values()):raise ValueError('Source neck has nonmanifold incidence')
 root_edges={edge for edge,owners in edge_faces.items()if len(owners)==1 and set(edge)<=root};starts={i for edge in root_edges for i in edge_faces[edge]};choices=[]
 if len(root_edges)!=len(root):raise ValueError('Declared root is not an actual complete boundary')
 for raw in profile_candidates:
  ring=inv[np.asarray(raw,int)];cut=set(map(int,ring))
  if cut&root or len(cut)!=len(ring):continue
  barrier={tuple(sorted((int(a),int(b))))for a,b in zip(ring,np.roll(ring,-1))}
  if not barrier<=set(edge_faces):continue
  adjacency=[[]for _ in range(len(f))]
  for edge,owners in edge_faces.items():
   if len(owners)==2 and edge not in barrier:
    a,b=owners;adjacency[a].append(b);adjacency[b].append(a)
  seen=set();stack=list(starts)
  while stack:
   i=stack.pop()
   if i not in seen:seen.add(i);stack.extend(adjacency[i])
  ids=np.array(sorted(seen),int)
  if not len(ids)or len(ids)==len(f):continue
  try:loops=boundary_loops(ff[ids])
  except ValueError:continue
  if len(loops)!=2 or {frozenset(x)for x in loops}!={frozenset(root),frozenset(cut)}:continue
  audit=continuity_audit(physical,ff[ids],expected_mouths=2,expected_components=1,expected_euler=0)
  if not audit['topological_continuity_pass']:continue
  kept=np.setdiff1d(np.arange(len(f)),ids);kept_loops=boundary_loops(ff[kept]);distal=next((x for x in kept_loops if set(x)==cut),None)
  if distal is None:continue
  choices.append(dict(sourceFaceIds=ids.tolist(),cutRawIds=representative[distal].tolist(),annularAudit=audit))
 if not choices:raise ValueError('No branch-free source sleeve ends at a declared profile')
 best=max(choices,key=lambda x:len(x['sourceFaceIds']));inside=set(best['sourceFaceIds']);nested=sorted([x for x in choices if set(x['sourceFaceIds'])<=inside],key=lambda x:len(x['sourceFaceIds']));previous=set()
 for row in nested:
  current=set(row['sourceFaceIds'])
  if not previous<=current:raise ValueError('Declared source sections do not form one nested annular sleeve')
  previous=current
 best['sourceSectionProfiles']=[dict(rawIds=x['cutRawIds'],sourcePrefixFaces=len(x['sourceFaceIds']))for x in nested]
 return best
