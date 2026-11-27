"""Conservative surface-guarded transport of an existing spatial field.

A finite triangle cover yields a globally 1-Lipschitz lower bound on distance
to fixed surfaces. Its zero set contains every point of those triangles.
Every bounded Euler map is bi-Lipschitz when its explicit step bound is <1;
thus the continuous ambient map fixes those surfaces. Discrete interpolated
source faces still require self/contact checks; this does not establish lumen
patency, closed tissue occupancy, or anatomical acceptance.
"""
import numpy as np
from scipy.spatial import cKDTree

class TriangleSurfaceGuard:
 def __init__(self, points, radius_m):
  self.points=np.asarray(points,float)
  self.radius=float(radius_m)
  if self.points.ndim!=2 or self.points.shape[1]!=3 or not len(self.points) or not np.isfinite(self.points).all() or not np.isfinite(self.radius) or self.radius<=0:raise ValueError('Finite nonempty triangle cover and positive radius required')
  self.index=cKDTree(self.points)

 @classmethod
 def from_triangles(cls, triangles, *, radius_m, maximum_samples=2000000, anisotropic=False):
  triangles=np.asarray(triangles,float);radius=float(radius_m)
  if triangles.ndim!=3 or triangles.shape[1:]!=(3,3) or not len(triangles) or not np.isfinite(triangles).all() or not np.isfinite(radius) or radius<=0:raise ValueError('Finite nonempty source triangles and positive cover radius required')
  edge=np.linalg.norm(triangles-np.roll(triangles,1,axis=1),axis=2).max(1)
  divisions=np.maximum(1,np.ceil(edge/radius).astype(np.int64))
  costs=(divisions+1)*(divisions+2)//2;strips=np.zeros(len(triangles),bool);ordered=triangles;row_divisions=np.ones(len(triangles),int)
  if anisotropic:
   lengths=np.linalg.norm(triangles[:,[1,2,0]]-triangles[:,[0,1,2]],axis=2);longest=lengths.argmax(1);orders=np.array([[0,1,2],[1,2,0],[2,0,1]]);ordered=triangles[np.arange(len(triangles))[:,None],orders[longest]]
   ab=ordered[:,1]-ordered[:,0];ac=ordered[:,2]-ordered[:,0];den=np.einsum('ij,ij->i',ab,ab);fraction=np.divide(np.einsum('ij,ij->i',ab,ac),den,out=np.zeros(len(den)),where=den>0);height=np.linalg.norm(ac-fraction[:,None]*ab,axis=1);spacing=radius/np.sqrt(2);row_divisions=np.maximum(1,np.ceil(height/spacing).astype(int))
   strip_upper=np.ceil((row_divisions+1)*(edge/(2*spacing)+2)).astype(np.int64);strips=strip_upper<costs;costs=np.where(strips,strip_upper,costs)
  count=int(costs.sum())
  if count>maximum_samples:raise ValueError(f'Triangle cover needs at most {count} samples; explicit cap is {maximum_samples}')
  patches=[]
  for n in np.unique(divisions[~strips]):
   t=triangles[(divisions==n)&~strips];i,j=np.triu_indices(int(n)+1);u=i/n;v=(j-i)/n
   # Each small lattice triangle has maximum edge <= radius.
   patches.append((t[:,None,0]+u[None,:,None]*(t[:,None,1]-t[:,None,0])+v[None,:,None]*(t[:,None,2]-t[:,None,0])).reshape(-1,3))
  for index in np.flatnonzero(strips):
   t=ordered[index];n=int(row_divisions[index]);levels=np.arange(n+1)/n;length=edge[index]*(1-levels);steps=np.ceil(length/spacing).astype(int);counts=steps+1;total=int(counts.sum());assert total<=costs[index];levels=np.repeat(levels,counts);starts=np.repeat(np.cumsum(counts)-counts,counts);local=np.arange(total)-starts;steps=np.repeat(steps,counts);fraction=np.divide(local,steps,out=np.zeros(total),where=steps>0)
   patches.append(t[0]+levels[:,None]*(t[2]-t[0])+((1-levels)*fraction)[:,None]*(t[1]-t[0]))
   # Longest-edge coordinates put the apex above the base interval. A
   # point projects into the next lower row's segment. Vertical distance
   # <= spacing and nearest-row horizontal distance <= spacing/2 imply
   # nearest cover distance <= sqrt(1.25)*spacing < radius.
  points=np.unique(np.concatenate(patches),axis=0)
  margin=radius*1e-10+1e-12
  guard=cls(points,radius+margin)
  guard.evidence=dict(sourceTriangles=len(triangles),coverPoints=len(points),maximumSubdivisionEdge_m=float((edge/divisions).max()),coverRadius_m=radius,subtractionRadius_m=guard.radius,roundoffCoverMargin_m=margin,maximumSamples=maximum_samples,zeroSetContainsAllSourceTriangles=True,distanceLowerBoundLipschitz=1.,sourceTrianglesUnmodified=True,anisotropicStripTriangles=int(strips.sum()),maximumRawCoverPoints=count)
  return guard

 def lower_distance_and_gradient(self, points):
  q=np.asarray(points,float)
  if q.ndim!=2 or q.shape[1]!=3 or not np.isfinite(q).all():raise ValueError('Finite N by 3 query controls required')
  distance,ids=self.index.query(q,k=1,eps=0,workers=2)
  lower=np.maximum(0.,distance-self.radius)
  grad=np.zeros_like(q);take=lower>0
  grad[take]=(q[take]-self.points[ids[take]])/distance[take,None]
  return lower,grad


def guarded_translation_flow(points, fields, guard, *, fade_m, steps, skin_band_guard=None, band_inner_m=None, band_outer_m=None):
 q=np.asarray(points,float).copy();fade=float(fade_m)
 if q.ndim!=2 or q.shape[1]!=3 or not np.isfinite(q).all() or not np.isfinite(fade) or fade<=0 or isinstance(steps,bool) or int(steps)!=steps or not 1<=steps<=256 or not fields:raise ValueError('Finite controls, positive fade and 1..256 steps required')
 steps=int(steps);velocity_L=sum(f.lipschitz_bound+np.linalg.norm(f.displacement)/fade for f in fields)
 if skin_band_guard is not None:
  if not np.isfinite([band_inner_m,band_outer_m]).all()or not 0<=band_inner_m<band_outer_m:raise ValueError('Ordered finite skin-band distances required')
  velocity_L+=sum(np.linalg.norm(f.displacement)for f in fields)*2.109375/(band_outer_m-band_inner_m)

 if velocity_L/steps>=1:raise ValueError('Euler step cannot certify injectivity')
 J=np.broadcast_to(np.eye(3),(len(q),3,3)).copy()
 for _ in range(steps):
  velocity=np.zeros_like(q);dv=np.zeros_like(J)
  for field in fields:
   if hasattr(field,'map_and_jacobians'):
    mapped,derivative=field.map_and_jacobians(q)
   else:mapped,derivative=field.map(q),field.jacobians(q)
   velocity+=mapped-q;dv+=derivative-np.eye(3)
  active=np.any(velocity!=0,axis=1)
  if not active.any():break
  ids=np.flatnonzero(active);v=velocity[ids];derivative=dv[ids]
  if skin_band_guard is not None:
   depth,normal=skin_band_guard.lower_distance_and_gradient(q[ids]);t=np.clip((depth-band_inner_m)/(band_outer_m-band_inner_m),0,1);band=(1-t)**4*(4*t+1);dband=np.where((depth>band_inner_m)&(depth<band_outer_m),-20*t*(1-t)**3/(band_outer_m-band_inner_m),0.)[:,None]*normal
   derivative=band[:,None,None]*derivative+v[:,:,None]*dband[:,None,:];v=band[:,None]*v;keep=band>0;ids=ids[keep];v=v[keep];derivative=derivative[keep]
   if not len(ids):break
  d,g=guard.lower_distance_and_gradient(q[ids]);ex=np.exp(-d/fade);weight=-np.expm1(-d/fade);dw=ex[:,None]*g/fade
  step_J=np.eye(3)+(weight[:,None,None]*derivative+v[:,:,None]*dw[:,None,:])/steps
  q[ids]+=weight[:,None]*v/steps
  J[ids]=step_J@J[ids]
 proof=dict(kind='surface-guarded-positive-translation-flow',steps=steps,fade_m=fade,velocityLipschitzBound=float(velocity_L),stepLipschitzBound=float(velocity_L/steps),minimumStretchLowerBound=float((1-velocity_L/steps)**steps),analytic_determinant_lower_bound=float((1-velocity_L/steps)**(3*steps)),fixedGuardSurfaces=True,barrierUsesConservativeTriangleCover=True,discreteMeshChecksPending=True,skinBandRestricted=skin_band_guard is not None,skinBandInner_m=band_inner_m,skinBandOuter_m=band_outer_m)
 return q,J,proof


class EllipsoidalTranslationField:
 """Rigid translation core with a finite C2 ellipsoidal transition."""
 def __init__(self, center_m, support_axes_m, core_fraction, displacement_m):
  from .compact_positive_fields import compact_ellipsoidal_translation
  self.center=np.asarray(center_m,float);self.axes=np.asarray(support_axes_m,float)
  self.core=float(core_fraction);self.displacement=np.asarray(displacement_m,float)
  self._apply=compact_ellipsoidal_translation
  self._apply(np.empty((0,3)),center=self.center,support_axes_m=self.axes,core_fraction=self.core,translation_m=self.displacement)
  self.lipschitz_bound=float(2.109375*np.linalg.norm(self.displacement)*np.linalg.norm(1/self.axes)/(1-self.core))
 def map_and_jacobians(self,points):
  q,J,_=self._apply(points,center=self.center,support_axes_m=self.axes,core_fraction=self.core,translation_m=self.displacement)
  return q,J
 def map(self,points):return self.map_and_jacobians(points)[0]
 def jacobians(self,points):return self.map_and_jacobians(points)[1]


class SeparableCorridorField:
 """Constant straight-corridor translation with a conservative full gradient bound."""
 def __init__(self, parameters, permutation):
  from .compact_positive_fields import separable_corridor_translation
  self.parameters=parameters;self.permutation=np.asarray(permutation,int);self.inverse=np.argsort(self.permutation)
  if sorted(self.permutation.tolist())!=[0,1,2]:raise ValueError('A coordinate permutation is required')
  stations=np.asarray([[s['y_m'],s['x_m'],s['z_m']]for s in parameters['centerline_stations']],float)
  if np.ptp(stations[:,1:],axis=0).max()!=0 or 'translation_stations'in parameters:raise ValueError('Only constant straight-corridor translations have this bound')
  core=np.asarray(parameters['core_axes_m']);support=np.asarray(parameters['support_axes_m']);width=support-core
  if 'y_core_interval_m'in parameters:
   yc=parameters['y_core_interval_m'];ys=parameters['y_support_interval_m'];width[1]=min(yc[0]-ys[0],ys[1]-yc[1])
  if np.any(width<=0):raise ValueError('Positive transition widths required')
  self.displacement=np.asarray(parameters['translation_m'])[self.inverse]
  self.lipschitz_bound=float(np.linalg.norm(self.displacement)*2.109375*np.linalg.norm(1/width))
  self._apply=separable_corridor_translation
 def map_and_jacobians(self, points):
  q,J,_=self._apply(np.asarray(points)[:,self.permutation],**self.parameters)
  return q[:,self.inverse],J[:,self.inverse,:][:,:,self.inverse]
 def map(self, points):
  return self.map_and_jacobians(points)[0]
 def jacobians(self, points):
  return self._apply(np.asarray(points)[:,self.permutation],**self.parameters)[1][:,self.inverse,:][:,:,self.inverse]


def swept_surface_guard_mask(protected_triangles, moving_triangles, displacement_bounds_m, *, maximum_proposals=12000000):
 """Retain every fixed triangle that could meet a bounded moving source face.

 Query boxes enclose each entire source triangle plus the componentwise
 maximum movement of all its vertices. They are broad-phase volumes only.
 This does not infer anatomical attachment or declare a collision.
 """
 from .triangle_cross_candidates import TriangleCrossIndex
 protected=np.asarray(protected_triangles,float);moving=np.asarray(moving_triangles,float);delta=np.asarray(displacement_bounds_m,float)
 if delta.shape!=(3,)or np.any(delta<0)or not np.isfinite(delta).all():raise ValueError('Finite nonnegative componentwise movement bound required')
 if moving.ndim!=3 or moving.shape[1:]!=(3,3)or not len(moving):raise ValueError('Nonempty actual moving source faces required')
 padding=1e-12;low=moving.min(1)-delta-padding;high=moving.max(1)+delta+padding;center=(low+high)/2
 boxes=np.stack([low,high,center],axis=1);index=TriangleCrossIndex(protected);keep=np.zeros(len(protected),bool);count=0
 for pairs in index.batches(boxes,maximum_proposals=maximum_proposals):
  keep[pairs[:,1]]=True;count+=len(pairs)
 return keep,dict(complete=True,movingSourceFaces=len(moving),protectedSourceFaces=len(protected),retainedProtectedFaces=int(keep.sum()),displacementBounds_m=delta.tolist(),roundoffBoundsPadding_m=padding,proposalPairs=count,maximumProposals=maximum_proposals,allPotentialMovingFaceContactsCovered=True,sourceTrianglesUnmodified=True)
