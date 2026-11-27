"""Remove near-collinear facets by flipping only numerically planar quads.

Controls stay literal. A replacement must preserve the convex patch's
orientation and materially improve both triangle altitudes. Full geometry
and source-attribute review remain the caller's responsibility.
"""
from collections import defaultdict
import time
import numpy as np

def flip_planar_slivers(vertices,faces,*,minimum_altitude_m,plane_tolerance_m=2e-8,maximum_flips=200,seconds=10.):
 v=np.asarray(vertices,float);f=np.array(faces,int,copy=True);start=time.monotonic();edges=defaultdict(set)
 def keys(t):return [tuple(sorted((int(a),int(b))))for a,b in zip(t,np.roll(t,-1))]
 for i,t in enumerate(f):
  for e in keys(t):edges[e].add(i)
 def altitude(t):
  p=v[t];return np.linalg.norm(np.cross(p[1]-p[0],p[2]-p[0]))/max(np.linalg.norm(p-np.roll(p,1,axis=0),axis=1))
 changes=[];blocked=set()
 for _ in range(maximum_flips):
  if time.monotonic()-start>seconds:raise TimeoutError('Planar diagonal work reached its explicit deadline')
  p=v[f];heights=np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)/np.linalg.norm(p-np.roll(p,1,axis=1),axis=2).max(1);bad=[int(i)for i in np.flatnonzero(heights<minimum_altitude_m)if int(i)not in blocked]
  if not bad:break
  i=min(bad,key=lambda k:heights[k]);best=None
  for a,b,c in [f[i],np.roll(f[i],-1),np.roll(f[i],-2)]:
   peers=edges[tuple(sorted((int(a),int(b))))]-{i}
   if len(peers)!=1:continue
   j=next(iter(peers));other=f[j];d=next((int(q)for q in other if q not in [a,b]),None)
   if d is None or d==c or tuple(sorted((int(c),d)))in edges:continue
   if not any(np.array_equal(np.roll(other,-k)[:2],[b,a])for k in range(3)):continue
   old=v[np.r_[f[i],f[j]].reshape(2,3)];norm=np.cross(old[:,1]-old[:,0],old[:,2]-old[:,0]);reference=norm[np.argmax(np.linalg.norm(norm,axis=1))];reference/=np.linalg.norm(reference)
   residual=float(abs((v[[a,b,c,d]]-v[a])@reference).max())
   if residual>plane_tolerance_m:continue
   proposal=np.array([[c,a,d],[c,d,b]],int);new=v[proposal];normal=np.cross(new[:,1]-new[:,0],new[:,2]-new[:,0])
   if np.any(normal@reference<=0):continue
   height=min(map(altitude,proposal))
   if height<=min(heights[i],heights[j])*1.5 or height<=heights[i]*1.5:continue
   choice=(height,j,proposal,residual)
   if best is None or height>best[0]:best=choice
  if best is None:blocked.add(i);continue
  height,j,proposal,residual=best
  for k in [i,j]:
   for e in keys(f[k]):
    edges[e].remove(k)
    if not edges[e]:del edges[e]
  f[[i,j]]=proposal
  for k in [i,j]:
   for e in keys(f[k]):edges[e].add(k)
  changes.append(dict(faces=[i,j],minimum_new_altitude_m=float(height),maximum_plane_residual_m=residual));blocked.discard(j)
 p=v[f];heights=np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1)/np.linalg.norm(p-np.roll(p,1,axis=1),axis=2).max(1)
 return f,dict(controls_literal=True,flips=changes,remaining_below_threshold=int(np.sum(heights<minimum_altitude_m)),minimum_altitude_m=float(heights.min()),plane_tolerance_m=plane_tolerance_m,elapsed_seconds=time.monotonic()-start,full_audit_pending=True)
