"""C2 compact translation field with a certified continuous Jacobian bound."""
import numpy as np


def compact_axis_field(points,centre,radii,amplitude,*,axis=2,gate=None):
 p,c,r=np.asarray(points,float),np.asarray(centre,float),np.asarray(radii,float)
 if p.ndim!=2 or p.shape[1]!=3 or c.shape!=(3,)or r.shape!=(3,)or axis not in [0,1,2]or not np.isfinite(np.r_[p.ravel(),c,r,amplitude]).all()or np.any(r<=0):raise ValueError('Finite points, centre, positive radii and a coordinate axis required')
 # max |d(1-u.u)^3/du_axis| = 96/(25 sqrt(5)); rank-one J.
 lower=1-abs(amplitude)/r[axis]*96/(25*np.sqrt(5))
 if lower<=0:raise ValueError('Requested field lacks a positive continuous Jacobian bound')
 u=(p-c)/r;s=np.sum(u*u,axis=1);a=np.maximum(0,1-s);weight=a**3;gradient=-6*a[:,None]**2*u/r
 if gate is not None:
  ga,sign,start,width=gate['axis'],gate['sign'],gate['start_m'],gate['width_m']
  if ga not in [0,1,2]or ga==axis or sign not in [-1,1]or not np.isfinite([start,width]).all()or width<=0:raise ValueError('A finite smooth gate perpendicular to the displacement axis is required')
  t=np.clip((sign*p[:,ga]-start)/width,0,1);g=t**3*(10+t*(-15+6*t));dg=30*t**2*(1-t)**2*sign/width
  gradient*=g[:,None];gradient[:,ga]+=weight*dg;weight*=g
  # Gate has no derivative along the displacement axis, so the same rank-one
  # determinant bound remains valid despite the additional transverse shear.
 q=p.copy();q[:,axis]+=amplitude*weight;J=np.repeat(np.eye(3)[None],len(p),axis=0);J[:,axis,:]+=amplitude*gradient
 return dict(points=q,jacobians=J,supported=weight>0,continuous_determinant_lower_bound=float(lower))
