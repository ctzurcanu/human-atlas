"""Conforming local subtraction of a convex prismatic opening from a surface.

The cutter is a projected polygon plus a bounded axis neighborhood. Original
controls are kept literally. Edge intersections are calculated from canonical
source-edge/plane identities; adjacent triangles therefore share one control.
This constructs an opening, not a solid or anatomical identity.
"""
from __future__ import annotations
import collections
import numpy as np

class SurfaceChannelError(ValueError):pass

def area(p):
 p=np.asarray(p);return float(np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1))/2) if len(p)>=3 else 0

def clean(poly):
 out=[]
 for p in poly:
  if not out or np.linalg.norm(np.asarray(p)-out[-1])>1e-14:out.append(np.asarray(p))
 if len(out)>1 and np.linalg.norm(out[-1]-out[0])<=1e-14:out.pop()
 return np.asarray(out).reshape(-1,2)

def clip(poly,c):
 p=np.asarray(poly);d=p@c[:2]+c[2];out=[]
 for i,j in zip(range(len(p)),np.roll(np.arange(len(p)),-1)):
  if d[i]>=-1e-14:out.append(p[i])
  if (d[i]>=-1e-14)!=(d[j]>=-1e-14):out.append(p[i]+(p[j]-p[i])*d[i]/(d[i]-d[j]))
 return clean(out)

def subtract(poly,planes):
 inside=np.asarray(poly);pieces=[]
 for c in planes:
  outside=clip(inside,-c)
  if len(outside)>=3 and abs(area(outside))>1e-20:pieces.append(outside)
  inside=clip(inside,c)
  if len(inside)<3:break
 return pieces

class PointBank:
 def __init__(self,positions,normals=None):
  self.source=np.asarray(positions,dtype=np.float32).copy();self.positions=self.source.tolist();self.precise=self.source.astype(float).tolist();self.normals=(np.zeros_like(self.source) if normals is None else np.asarray(normals,dtype=np.float32)).tolist();self.lookup={tuple(p):i for i,p in enumerate(self.source)};self.ancestry={};self.edge_points=collections.defaultdict(set)
 def add(self,p,normal=None,ancestry=None):
  p=np.asarray(p,dtype=np.float32);key=tuple(p);idx=self.lookup.get(key)
  if idx is not None:return idx
  idx=len(self.positions);self.lookup[key]=idx;self.positions.append(p.tolist());self.precise.append(np.asarray(ancestry.get('precisePositionMeters',p) if isinstance(ancestry,dict) else p,dtype=float).tolist());self.normals.append(([0,0,0] if normal is None else np.asarray(normal,dtype=np.float32).tolist()));self.ancestry[idx]=ancestry;return idx


def earcut(poly,coordinates):
 # Retain collinear boundary controls: a diagonal may not skip any such point.
 ids=[]
 for j in poly:
  if not ids or j!=ids[-1]:ids.append(int(j))
 if len(ids)>1 and ids[0]==ids[-1]:ids.pop()
 if len(ids)<3:return []
 p=np.array([coordinates[j] for j in ids]);sign=1 if area(p)>0 else -1;remaining=list(range(len(ids)));out=[]
 def cross(a,b,c):q=b-a;r=c-a;return q[0]*r[1]-q[1]*r[0]
 while len(remaining)>3:
  found=False
  for k in range(len(remaining)):
   a,b,c=remaining[k-1],remaining[k],remaining[(k+1)%len(remaining)]
   if sign*cross(p[a],p[b],p[c])<=1e-22:continue
   blocked=False
   for j in remaining:
    if j in [a,b,c]:continue
    if min(sign*cross(p[a],p[b],p[j]),sign*cross(p[b],p[c],p[j]),sign*cross(p[c],p[a],p[j]))>=-1e-22:blocked=True;break
   if not blocked:out.append([ids[a],ids[b],ids[c]]);remaining.pop(k);found=True;break
  if not found:raise SurfaceChannelError('Local conforming polygon cannot be triangulated')
 if abs(cross(p[remaining[0]],p[remaining[1]],p[remaining[2]]))>1e-22:out.append([ids[j] for j in remaining])
 return out


def subtract_surface(bank,faces,center,basis,cutter,axis_interval,port_id):
 faces=np.asarray(faces,dtype=np.uint32);source=np.array(bank.positions,dtype=float);uv=(source-center)@basis[:,:2];depth=(source-center)@basis[:,2];cutter=np.asarray(cutter,dtype=float)
 if area(cutter)<0:cutter=cutter[::-1]
 planes=[]
 for a,b in zip(cutter,np.roll(cutter,-1,axis=0)):
  e=b-a;c=np.array([-e[1],e[0],e[1]*a[0]-e[0]*a[1]]);c/=np.linalg.norm(c[:2]);planes.append(c)
 planes=np.array(planes);polygons={};changed=[];canonical={};coordinates={i:uv[i] for i in range(len(source))}
 def add_on_triangle(fi,q):
  tri=faces[fi];xy=uv[tri];ab=xy[1]-xy[0];ac=xy[2]-xy[0];w=np.linalg.solve(np.column_stack([ab,ac]),q-xy[0]);weights=np.array([1-w.sum(),*w]);near=np.flatnonzero(abs(weights)<1e-10)
  exact=np.flatnonzero(np.linalg.norm(xy-q,axis=1)<1e-13)
  if len(exact):return int(tri[exact[0]])
  edge=None
  if len(near):
   ids=[j for j in range(3) if j!=int(near[0])];a,b=sorted(int(tri[j]) for j in ids);line_dist=abs(planes[:,:2]@q+planes[:,2]);plane=int(line_dist.argmin());c=planes[plane];t=-(uv[a]@c[:2]+c[2])/((uv[b]-uv[a])@c[:2]);key=(a,b,port_id,plane)
   if key in canonical:return canonical[key]
   p=source[a]+(source[b]-source[a])*t;n=np.asarray(bank.normals[a])*(1-t)+np.asarray(bank.normals[b])*t;edge=(a,b);ancestry={'sourceEdge':[a,b],'fraction':float(t),'port':port_id,'cutterPlane':plane}
  else:
   # A remaining interior point is an exact polygon-corner ray/face crossing.
   distances=abs(planes[:,:2]@q+planes[:,2]);pair=sorted(np.argsort(distances)[:2].tolist());assert distances[pair].max()<1e-10
   q=np.linalg.solve(planes[pair,:2],-planes[pair,2]);w=np.linalg.solve(np.column_stack([ab,ac]),q-xy[0]);weights=np.array([1-w.sum(),*w]);key=(fi,port_id,'planeIntersection',*pair);p=weights@source[tri];n=weights@np.asarray(bank.normals)[tri];ancestry={'sourceFace':int(fi),'barycentric':weights.tolist(),'port':port_id,'cutterPlaneIntersection':pair}
  ancestry['precisePositionMeters']=np.asarray(p,dtype=float).tolist();idx=bank.add(p,n,ancestry);canonical[key]=idx;coordinates[idx]=q
  if edge:bank.edge_points[edge].add(idx)
  return idx
 for fi,tri in enumerate(faces):
  xy=uv[tri]
  if np.max(depth[tri])<axis_interval[0] or np.min(depth[tri])>axis_interval[1]:continue
  if np.any(xy.max(0)<cutter.min(0)-1e-14) or np.any(xy.min(0)>cutter.max(0)+1e-14):continue
  if abs(area(xy))<1e-18:continue
  overlap=xy.copy()
  for c in planes:
   overlap=clip(overlap,c)
   if len(overlap)<3:break
  if len(overlap)<3 or abs(area(overlap))<=1e-20:continue
  pieces=subtract(xy,planes);polygons[fi]=[[add_on_triangle(fi,q) for q in poly] for poly in pieces];changed.append(fi)
 # Split adjacent source faces and all local piece edges at EVERY shared
 # control. Only source faces touching such a local split are retessellated.
 output=[];parents=[];retess=[]
 for fi,tri in enumerate(faces):
  pieces=polygons.get(fi,[tri.tolist()]);points=set(j for poly in pieces for j in poly)
  for a,b in zip(tri,np.roll(tri,-1)):
   points.update(bank.edge_points.get(tuple(sorted((int(a),int(b)))),set()))
  localcoords={j:coordinates.get(j,(np.array(bank.positions[j])-center)@basis[:,:2]) for j in points}
  refined=[]
  for poly in pieces:
   new=[]
   for a,b in zip(poly,np.roll(poly,-1)):
    qa,qb=localcoords[a],localcoords[b];e=qb-qa;ll=e@e;new.append(a)
    if ll==0:continue
    knots=[]
    for j in points:
     if j in [a,b]:continue
     t=np.dot(localcoords[j]-qa,e)/ll
     if 1e-12<t<1-1e-12 and np.linalg.norm(localcoords[j]-(qa+t*e))<1e-11:knots.append((t,j))
    new.extend(j for t,j in sorted(knots))
   refined.append(new)
  if fi not in polygons and len(refined)==1 and refined[0]==tri.tolist():output.append(tri.tolist());parents.append(fi);continue
  if fi not in polygons:retess.append(fi)
  for poly in refined:
   for t in earcut(poly,localcoords):output.append(t);parents.append(fi)
 return np.asarray(output,dtype=np.uint32),np.asarray(parents,dtype=np.uint32),{'sourceFacesCut':changed,'adjacentSourceFacesRetessellated':retess,'unchangedSourceFaces':len(faces)-len(changed)-len(retess),'port':port_id,'openingOnly':True}
