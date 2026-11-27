"""Bounded cubic Hermite stem proposals with explicit physical bend limits.

Curvature extrema are replayed at every real stationary polynomial root.
Float64 polynomial roots are numerical evidence, not an exact-arithmetic proof
or a substitute for complete emitted-surface and anatomical checks.
"""
import time
import numpy as np
from numpy.polynomial import polynomial as poly
from scipy.optimize import minimize_scalar
from .tube_sections import unit


def _critical(coefficients):
 c=poly.polytrim(np.asarray(coefficients,float),tol=1e-24)
 if len(c)<2:return np.empty(0)
 roots=poly.polyroots(c);real=roots[np.abs(roots.imag)<1e-7].real
 return real[(real>0)&(real<1)]


def cubic_metrics(coefficients):
 c=np.asarray(coefficients,float);v=np.array([poly.polyder(x)for x in c.T]).T;a=np.array([poly.polyder(x)for x in v.T]).T
 def add_squares(values):
  result=np.zeros(1)
  for x in values:result=poly.polyadd(result,poly.polymul(x,x))
  return result
 speed=add_squares(v.T)
 cross=np.zeros((4,3))
 for k,(i,j)in enumerate([(1,2),(2,0),(0,1)]):
  value=poly.polysub(poly.polymul(v[:,i],a[:,j]),poly.polymul(v[:,j],a[:,i]));cross[:len(value),k]=value
 numerator=add_squares(cross.T)
 stationary=poly.polysub(poly.polymul(poly.polyder(numerator),speed),3*poly.polymul(numerator,poly.polyder(speed)))
 samples=np.r_[0.,_critical(stationary),1.];s=poly.polyval(samples,speed);n=poly.polyval(samples,numerator)
 minimum_speed2=float(poly.polyval(np.r_[0.,_critical(poly.polyder(speed)),1.],speed).min())
 if minimum_speed2<=1e-24:return dict(minimumRadius_m=0.,minimumSpeed2=minimum_speed2,stationaryParameters=samples.tolist())
 k2=np.maximum(n,0)/s**3;maximum=float(np.sqrt(k2.max(initial=0)));radius=1/maximum if maximum else 1e30
 dense=np.linspace(0,1,257);length=float(np.trapezoid(np.sqrt(np.maximum(poly.polyval(dense,speed),0)),dense))
 return dict(minimumRadius_m=radius,minimumSpeed2=minimum_speed2,stationaryParameters=samples.tolist(),length_m=length)


def port_hermite_course(start,start_tangent,end,end_tangent,*,minimum_radius_m,
                       maximum_length_ratio=3.,wall_seconds=.5,maximum_evaluations=48):
 p0,p1=np.asarray(start,float),np.asarray(end,float);t0,t1=unit(start_tangent),unit(end_tangent);distance=float(np.linalg.norm(p1-p0))
 if p0.shape!=(3,)or p1.shape!=(3,)or not np.isfinite(np.r_[p0,p1,t0,t1,minimum_radius_m,maximum_length_ratio,wall_seconds]).all()or distance<1e-10 or minimum_radius_m<=0 or maximum_length_ratio<1 or not 0<wall_seconds<=3:raise ValueError('Finite declared port dimensions and bounded course limits required')
 began=time.monotonic();evaluations=[]
 def evaluate(relative):
  if time.monotonic()-began>wall_seconds:raise TimeoutError('Native port-course proposal deadline')
  handle=relative*distance;v0=t0*handle;v1=t1*handle;c=np.array([p0,v0,3*(p1-p0)-2*v0-v1,2*(p0-p1)+v0+v1]);metrics=cubic_metrics(c)
  allowed=metrics['length_m']<=maximum_length_ratio*distance
  score=-metrics['minimumRadius_m']/distance if allowed else 1.+metrics['length_m']/distance
  evaluations.append((score,c,metrics,relative));return score
 result=minimize_scalar(evaluate,bounds=(.1,2.),method='bounded',options=dict(maxiter=maximum_evaluations,xatol=1e-5))
 _,c,metrics,relative=min(evaluations,key=lambda x:x[0]);accepted=metrics['minimumRadius_m']>=minimum_radius_m*1.01 and metrics['length_m']<=maximum_length_ratio*distance
 evidence=dict(complete=bool(result.success),accepted=bool(result.success and accepted),**metrics,
  minimumRequestedRadius_m=float(minimum_radius_m),maximumLengthRatio=maximum_length_ratio,chordLength_m=distance,handleRatio=float(relative),evaluations=len(evaluations),
  criterion='Numerical stationary polynomial replay with one-percent bend-radius margin; emitted surface checks pending',
  anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-began)
 return c,evidence
