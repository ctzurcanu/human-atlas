"""Minimum-volume axis-aligned support for a complete material-point cohort.

In normalized coordinates, minimize -sum(log(diagonal A)) subject to
||A*p+b|| <= radius_limit. This convex formulation fixes the symmetry-axis
center and may constrain the support's inferior extent. It fits a support
domain only, never an anatomical mesh or a clinical correspondence.
"""
import time
import numpy as np
from scipy.optimize import minimize
from scipy.spatial import ConvexHull

def fit_support(points,center,axes,*,radius_limit=.8,minimum_y_m=None,wall_seconds=2.):
 start=time.monotonic();p=np.asarray(points,float);center=np.asarray(center,float);axes=np.asarray(axes,float)
 if p.ndim!=2 or p.shape[1]!=3 or len(p)<4 or center.shape!=(3,)or axes.shape!=(3,)or not np.isfinite(np.r_[p.ravel(),center,axes]).all()or np.any(axes<=0)or not 0<radius_limit<1:raise ValueError('Finite 3D points and positive enclosing support required')
 q=(p-center)/axes;hull=ConvexHull(q);q=q[hull.vertices]
 def guard():
  if time.monotonic()-start>wall_seconds:raise TimeoutError('Finite convex support-fit deadline')
 def values(x):guard();return q*x[:3]+np.array([0.,x[3],x[4]])
 def constraints(x):return radius_limit**2-(values(x)**2).sum(1)
 def jac(x):
  v=values(x);return np.c_[-2*v*q,-2*v[:,1],-2*v[:,2]]
 cons=[dict(type='ineq',fun=constraints,jac=jac)]
 if minimum_y_m is not None:
  k=(center[1]-minimum_y_m)/axes[1];cons.append(dict(type='ineq',fun=lambda x:k*x[1]-x[3]-1,jac=lambda x:np.array([0.,k,0.,-1.,0.])))
 result=minimize(lambda x:-np.log(x[:3]).sum(),np.array([1.,1.,1.,0.,0.]),jac=lambda x:np.r_[-1/x[:3],0.,0.],constraints=cons,bounds=[(.05,20.)]*3+[(None,None)]*2,method='SLSQP',options=dict(maxiter=200,ftol=1e-12));guard()
 if not result.success or constraints(result.x).min()<-1e-8:raise ValueError('Convex support enclosure did not converge: '+result.message)
 x=result.x;new_axes=axes/x[:3];new_center=center-np.array([0.,x[3],x[4]])*new_axes;r=np.linalg.norm((p-new_center)/new_axes,axis=1)
 if r.max()>radius_limit+1e-8 or minimum_y_m is not None and new_center[1]-new_axes[1]<minimum_y_m-1e-9:raise ValueError('Fitted support failed original-point/inferior-bound replay')
 return new_center,new_axes,dict(complete=True,algorithm='convex diagonal inverse-ellipsoid enclosure',pointCount=len(p),hullPointCount=len(q),normalizedMaximumRadius=float(r.max()),radiusLimit=radius_limit,volumeRatio=float(np.prod(new_axes)/np.prod(axes)),iterations=int(result.nit),minimumY_m=float(new_center[1]-new_axes[1]),centerXFixed=True,anatomicalAcceptance=False,elapsedSeconds=time.monotonic()-start)
