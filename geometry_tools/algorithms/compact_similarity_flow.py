"""Compact proper-similarity flow with an analytic global orientation bound.

Each small step is identity plus a globally Lipschitz displacement with
constant below one in ellipsoid coordinates. Its Jacobian is nonsingular
along the identity homotopy and has positive determinant. Composition keeps
that orientation; neither this certificate nor an ICP fit accepts anatomy.
"""
import numpy as np
from scipy.linalg import logm,expm

class CompactSimilarityFlow:
 def __init__(self,matrix,*,center,support_axes_m,core_fraction=.7,steps=64):
  m=np.asarray(matrix,float);c=np.asarray(center,float);axes=np.asarray(support_axes_m,float);core=float(core_fraction)
  if m.shape!=(4,4)or c.shape!=(3,)or axes.shape!=(3,)or not np.isfinite(np.r_[m.ravel(),c,axes,core]).all()or np.any(axes<=0)or not 0<=core<1 or not isinstance(steps,int)or not 1<=steps<=256:
   raise ValueError('Finite proper similarity, positive support axes, core fraction and 1–256 steps required')
  if not np.array_equal(m[3],[0,0,0,1]):raise ValueError('Homogeneous affine matrix required')
  a=m[:3,:3];scale=np.cbrt(np.linalg.det(a))
  if scale<=0 or not np.allclose(a.T@a,np.eye(3)*scale**2,rtol=1e-10,atol=1e-12):raise ValueError('Only positive uniform scale and proper rotation are accepted')
  centered=m.copy();centered[:3,3]=a@c+m[:3,3]-c;generator=logm(centered)
  if np.max(np.abs(np.imag(generator)))>1e-10:raise ValueError('Similarity logarithm is not numerically real')
  small=expm(np.real(generator)/steps);self.linear=small[:3,:3];self.c=(self.linear*axes[None,:])/axes[:,None]-np.eye(3);self.d=small[:3,3]/axes
  # For |q| <= 1, |Cq+d| <= ||C|| + |d|. Outside the support
  # displacement and derivative are exactly zero. Wendland max |w'|=135/64.
  norm=float(np.linalg.norm(self.c,2));lipschitz=norm+2.109375/(1-core)*(norm+float(np.linalg.norm(self.d)))
  if not lipschitz<1-1e-8:raise ValueError('Declared step count cannot certify the compact similarity')
  self.center=c.copy();self.axes=axes.copy();self.core=core;self.steps=steps;self.matrix=m.copy();self.step_lipschitz=lipschitz;self.log_lower=3*steps*np.log1p(-lipschitz);self.lower=float(np.exp(self.log_lower))
  if not self.lower>0:raise ValueError('Analytic determinant bound underflows; use a better-conditioned support')

 def apply(self,points):
  p=np.asarray(points,float)
  if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite N by 3 points required')
  out=p.copy();jac=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();q=(p-self.center)/self.axes;ids=np.flatnonzero(np.linalg.norm(q,axis=1)<1)
  if len(ids):
   q=q[ids].copy();j=np.broadcast_to(np.eye(3),(len(ids),3,3)).copy()
   for _ in range(self.steps):
    radius=np.linalg.norm(q,axis=1);t=np.clip((radius-self.core)/(1-self.core),0,1);w=(1-t)**4*(4*t+1);derivative=np.where((radius>self.core)&(radius<1),-20*t*(1-t)**3/(1-self.core),0.);gradient=np.divide(q,radius[:,None],out=np.zeros_like(q),where=radius[:,None]>0)*derivative[:,None];d=q@self.c.T+self.d;local=np.eye(3)[None]+w[:,None,None]*self.c[None]+d[:,:,None]*gradient[:,None,:];q+=w[:,None]*d;j=local@j
   out[ids]=self.center+q*self.axes;jac[ids]=j*self.axes[None,:,None]/self.axes[None,None,:]
  return out,jac,dict(kind='compact-proper-similarity-flow',steps=self.steps,step_displacement_lipschitz_bound=self.step_lipschitz,analytic_determinant_lower_bound=self.lower,log_determinant_lower_bound=float(self.log_lower),identity_outside_support=True,source_faces_changed=False,discrete_mesh_checks_pending=True,anatomicalAcceptance=False)
