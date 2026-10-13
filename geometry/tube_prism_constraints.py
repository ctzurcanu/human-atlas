"""Exact continuous prism inequalities for shared station translations.

Source faults require explicit diagnostic opt-in and remain constrained to their
original maximum; this is never a geometry-acceptance certificate.
Translations are equal at each corresponding layer control. Tissue depth
vectors remain fixed; each face spans precisely two adjacent stations.

Example (metres, outward carrier faces, inward wall boundaries):
    gate = StationPrismConstraints(fields, faces, station_ids, pinned, .05)
    maxima, gradient = gate.evaluate(station_translations)
    # gradient[k, face] differentiates with respect to high-minus-low motion.
    feasible = np.all(maxima <= gate.limit)
    candidate = fields + station_translations[station_ids][None]

Recheck the serialized candidate, source section measures, actual intersections
and adjoining ports before use. A feasible local determinant does not establish
whole-organ packing, lumen continuity or anatomical acceptance.
"""
import numpy as np
from .tube_sections import TubeGeometryError

class StationPrismConstraints:
 def __init__(self,fields,faces,station,protected,fraction=.1,*,allow_source_faults=False):
  fields=np.asarray(fields,dtype=float);faces=np.asarray(faces);station=np.asarray(station);protected=np.asarray(protected,dtype=bool)
  if fields.ndim!=3 or fields.shape[-1]!=3 or len(fields)<2 or not np.isfinite(fields).all():
   raise TubeGeometryError('finite wall fields must have shape [interfaces, controls, 3]')
  if faces.ndim!=2 or faces.shape[1]!=3 or not len(faces) or not np.issubdtype(faces.dtype,np.integer) or faces.min()<0 or faces.max()>=fields.shape[1]:
   raise TubeGeometryError('nonempty integer triangles must index the wall controls')
  if station.shape!=(fields.shape[1],) or not np.issubdtype(station.dtype,np.integer) or station.min()<0 or not np.array_equal(np.unique(station),np.arange(station.max()+1)):
   raise TubeGeometryError('one contiguous nonnegative station index per control required')
  if protected.shape!=(station.max()+1,):
   raise TubeGeometryError('one protected flag per station required')
  if not np.isfinite(fraction) or not 0<fraction<1:
   raise TubeGeometryError('source reserve fraction must lie between zero and one')
  self.faces=faces;self.station=station;self.protected=protected
  face_station=station[faces];self.low=face_station.min(1);self.high=face_station.max(1)
  if not np.all(self.high-self.low==1):
   raise TubeGeometryError('each carrier triangle must span exactly two adjacent stations')
  flag=(face_station-self.low[:,None]).astype(float);alpha=flag[:,1]-flag[:,0];beta=flag[:,2]-flag[:,0]
  aa=fields[:-1,faces].astype(float);bb=fields[1:,faces].astype(float);depth=bb-aa
  u=aa[:,:,1]-aa[:,:,0];v=aa[:,:,2]-aa[:,:,0];du=depth[:,:,1]-depth[:,:,0];dv=depth[:,:,2]-depth[:,:,0]
  self.c0=np.einsum('kfi,kfvi->kfv',np.cross(u,v),depth)
  self.c1=np.einsum('kfi,kfvi->kfv',np.cross(du,v)+np.cross(u,dv),depth)
  self.c2=np.einsum('kfi,kfvi->kfv',np.cross(du,dv),depth)
  self.g0=alpha[None,:,None,None]*np.cross(v[:,:,None],depth)+beta[None,:,None,None]*np.cross(depth,u[:,:,None])
  self.g1=alpha[None,:,None,None]*np.cross(dv[:,:,None],depth)+beta[None,:,None,None]*np.cross(depth,du[:,:,None])
  base,_=self.evaluate(np.zeros((int(station.max())+1,3)))
  self.sourceHealthy=base<0
  self.sourceInvalidCounts=(~self.sourceHealthy).sum(axis=1)
  if not allow_source_faults and not np.all(self.sourceHealthy):
   raise TubeGeometryError("source wall has invalid continuous prisms; no implicit acceptance")
  self.limit=np.where(self.sourceHealthy,base*fraction,base)
  self.sourceMaximum=base
 def evaluate(self,translation):
  translation=np.asarray(translation,dtype=float)
  if translation.shape!=(len(self.protected),3) or not np.isfinite(translation).all():
   raise TubeGeometryError('one finite 3D translation per station required')
  relative=translation[self.high]-translation[self.low]
  c0=self.c0+np.einsum('kfvi,fi->kfv',self.g0,relative)
  c1=self.c1+np.einsum('kfvi,fi->kfv',self.g1,relative)
  t=np.divide(-c1,2*self.c2,out=np.zeros_like(c1),where=self.c2!=0)
  middle=(self.c2!=0)&(t>0)&(t<1)
  start=c0;end=c0+c1+self.c2;interior=np.where(middle,c0+c1*t+self.c2*t*t,-np.inf)
  values=np.stack((start,end,interior),axis=-1)
  depth_choice=values.argmax(-1);corner_value=values.max(-1)
  corner=corner_value.argmax(-1)
  k,f=np.indices(corner.shape);maximum=corner_value[k,f,corner]
  active_t=np.where(depth_choice==0,0,np.where(depth_choice==1,1,t))[k,f,corner]
  gradient=self.g0[k,f,corner]+self.g1[k,f,corner]*active_t[:,:,None]
  return maximum,gradient
 def project(self,translation,sweeps=20):
  """Bounded diagnostic relaxation; callers must recheck all inequalities."""
  self.evaluate(translation)
  if not isinstance(sweeps,(int,np.integer)) or sweeps<0:raise TubeGeometryError('nonnegative integer sweep count required')
  out=np.asarray(translation,dtype=float).copy();out[self.protected]=0
  for _ in range(sweeps):
   value,g=self.evaluate(out);violation=value-self.limit;coats,faces=np.where(violation>0)
   if not len(faces):break
   low,high=self.low[faces],self.high[faces];gradient=g[coats,faces]
   movable_low=(~self.protected[low]).astype(float);movable_high=(~self.protected[high]).astype(float)
   denominator=np.sum(gradient*gradient,axis=1)*(movable_low+movable_high)
   use=denominator>1e-30;low,high,gradient,denominator,coats,faces,movable_low,movable_high=[x[use] for x in [low,high,gradient,denominator,coats,faces,movable_low,movable_high]]
   correction=gradient*(violation[coats,faces]/denominator)[:,None]
   delta=np.zeros_like(out);weights=np.zeros(len(out))
   np.add.at(delta,low,correction*movable_low[:,None]);np.add.at(delta,high,-correction*movable_high[:,None])
   np.add.at(weights,low,movable_low);np.add.at(weights,high,movable_high)
   out+=delta/np.maximum(weights[:,None],1);out[self.protected]=0
  return out
