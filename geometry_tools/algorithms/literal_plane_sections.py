"""Uncapped plane/surface contours joined only by exact source incidence.

An open source may still have valid local section loops. Open/branched cuts
are reported, never closed or welded by a distance tolerance.
"""
from collections import defaultdict
import numpy as np

def literal_plane_sections(vertices,faces,point,normal):
 v=np.asarray(vertices,float);f=np.asarray(faces,int);p=np.asarray(point,float);n=np.asarray(normal,float)
 if v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or not np.isfinite(v).all()or not np.isfinite(p).all()or not np.isfinite(n).all()or np.linalg.norm(n)==0:raise ValueError('Finite triangle surface and physical plane required')
 n=n/np.linalg.norm(n);physical,representative,inverse=np.unique(v,axis=0,return_index=True,return_inverse=True);f=inverse[f];s=(physical-p)@n;points={};edges={};coplanar=[];unresolved=[]
 def hit(a,b):
  if s[a]==0:key=('v',int(a));q=physical[a]
  elif s[b]==0:key=('v',int(b));q=physical[b]
  else:
   a,b=sorted((int(a),int(b)));key=('e',a,b);q=physical[a]+s[a]/(s[a]-s[b])*(physical[b]-physical[a])
  if key not in points:points[key]=q
  return key
 indices=np.flatnonzero((s[f].min(1)<=0)&(s[f].max(1)>=0))
 for fi in indices:
  tri=f[fi]
  if np.all(s[tri]==0):coplanar.append(int(fi));continue
  keys=[]
  for a,b in zip(tri,np.roll(tri,-1)):
   if s[a]*s[b]<0 or s[a]==0 or s[b]==0:
    key=hit(a,b)
    if key not in keys:keys.append(key)
  if len(keys)==2 and keys[0]!=keys[1]:edges.setdefault(tuple(sorted(keys)),[]).append(int(fi))
  elif len(keys)>2:unresolved.append(int(fi))
 graph=defaultdict(set)
 for a,b in edges:graph[a].add(b);graph[b].add(a)
 unseen=set(graph);loops=[];loop_source_faces=[];loop_raw_edges=[];loop_edge_weights=[];open_components=[]
 while unseen:
  todo=[min(unseen)];component=set()
  while todo:
   k=todo.pop()
   if k in component:continue
   component.add(k);todo.extend(graph[k]-component)
  unseen-=component
  if any(len(graph[k])!=2 for k in component):open_components.append(dict(points=np.array([points[k]for k in sorted(component)]),degrees=[len(graph[k])for k in sorted(component)]));continue
  first=min(component);order=[first];previous=None;current=first
  while True:
   candidates=graph[current]-({previous}if previous is not None else set());nxt=min(candidates)
   if nxt==first:break
   if nxt in order:raise ValueError('Non-simple source section incidence')
   order.append(nxt);previous,current=current,nxt
  loops.append(np.array([points[k]for k in order]))
  raw=[];weights=[]
  for key in order:
   if key[0]=='v':raw.append([representative[key[1]],representative[key[1]]]);weights.append([1.,0.])
   else:
    a,b=key[1:];t=s[a]/(s[a]-s[b]);raw.append([representative[a],representative[b]]);weights.append([1-t,t])
  loop_raw_edges.append(np.asarray(raw,int));loop_edge_weights.append(np.asarray(weights,float))
  loop_source_faces.append(sorted({fi for a,b in zip(order,order[1:]+order[:1])for fi in edges[tuple(sorted((a,b)))]}))
 axis=np.eye(3)[np.argmin(abs(n))];u=np.cross(n,axis);u/=np.linalg.norm(u);basis=np.column_stack((u,np.cross(n,u)))
 segments=np.array([[points[a],points[b]]for a,b in edges]).reshape(-1,2,3)
 return dict(point_m=p,normal=n,basis=basis,loops=loops,loop_source_faces=loop_source_faces,loop_raw_edges=loop_raw_edges,loop_edge_weights=loop_edge_weights,plane_loops=[(q-p)@basis for q in loops],segments=segments,plane_segments=(segments-p)@basis,segment_source_faces=list(edges.values()),open_components=open_components,coplanar_source_faces=coplanar,unresolved_source_faces=unresolved,section_source_faces=indices,closed_local_sections=not open_components and not coplanar and not unresolved,source_surface_closed_not_assumed=True,caps_created=0,source_geometry_changed=False,join_policy='exact source-coordinate incidence; no near-point welding')
