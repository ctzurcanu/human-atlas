"""Proper common section poses for explicitly named vascular bundles.

The source correspondence and anchor ownership must be supplied by the caller.
No radius scale, inferred attachment or clearance guarantee is introduced.
"""
import numpy as np
from scipy.interpolate import CubicSpline

def quintic(x):
 x=np.clip(x,0,1);return x*x*x*(10+x*(-15+6*x))

def member_translation(y, corrections):
 result=np.zeros(3)
 for correction in corrections:
  support=np.asarray(correction['ySupport_m'],float);delta=np.asarray(correction['translation_m'],float)
  if support.shape!=(4,)or not np.all(np.diff(support)>0)or delta.shape!=(3,)or not np.isfinite(np.r_[support,delta]).all():raise ValueError('Finite local source-member correction and ordered support required')
  a,b,c,d=support;result+=quintic((y-a)/(b-a))*(1-quintic((y-c)/(d-c)))*delta
 return result

class AxialBundlePose:
 def __init__(self, controls, support, anchors, *, transition_m=.030):
  ys=np.asarray(controls['y_m'],float);source=np.asarray(controls['sourcePairMidpointXZ_m'],float);target=np.asarray(controls['targetPairMidpointXZ_m'],float);angle=np.asarray(controls['angleRadians'],float);self.support=np.asarray(support,float);self.anchors={k:np.asarray(v,float).reshape(-1,2)for k,v in anchors.items()};self.transition=float(transition_m)
  if ys.ndim!=1 or len(ys)<2 or not np.all(np.diff(ys)>0)or source.shape!=(len(ys),2)or target.shape!=source.shape or angle.shape!=ys.shape:raise ValueError('Corresponding ordered source bone-frame stations required')
  if self.support.shape!=(4,)or not np.all(np.diff(self.support)>0)or self.transition<=0 or not np.isfinite(self.transition):raise ValueError('Ordered support and positive fixed-interface transition required')
  if not np.isfinite(np.r_[ys,source.ravel(),target.ravel(),angle,self.support]).all()or any(not np.isfinite(v).all()or np.any(v[:,1]<v[:,0])for v in self.anchors.values()):raise ValueError('Finite ordered anchor intervals and frame controls required')
  self.ys=ys;self.sc=CubicSpline(ys,source,bc_type='natural');self.tc=CubicSpline(ys,target,bc_type='natural');self.ac=CubicSpline(ys,angle,bc_type='natural')
 def weight(self,y,bundles):
  if not bundles or any(k not in self.anchors for k in bundles):raise ValueError('Explicit known bundle ownership required')
  a,b,c,d=self.support;base=float(quintic((y-a)/(b-a))*(1-quintic((y-c)/(d-c))));weights=[]
  for key in bundles:
   intervals=self.anchors[key];distance=np.maximum(np.maximum(intervals[:,0]-y,y-intervals[:,1]),0);weights.append(float(np.prod(quintic(distance/self.transition))))
  return base*min(weights)
 def pose(self,y,bundles):
  w=self.weight(float(y),bundles);station=np.clip(y,self.ys[0],self.ys[-1]);a,b=self.sc(station),self.tc(station);angle=float(self.ac(station))*w;c,s=np.cos(angle),np.sin(angle);R=np.array([[c,0,-s],[0,1,0],[s,0,c]]);origin=np.array([a[0],y,a[1]]);target=origin+np.array([w*(b[0]-a[0]),0,w*(b[1]-a[1])]);return R,target-R@origin,w
