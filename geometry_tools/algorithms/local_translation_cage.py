"""Compact, globally injective translation field with explicit stretch bound."""
from dataclasses import dataclass
import numpy as np

@dataclass
class LocalTranslationCage:
    windows:np.ndarray
    translation:np.ndarray
    minimum_jacobian:float=.2

    def __post_init__(self):
        self.windows=np.asarray(self.windows,float);self.translation=np.asarray(self.translation,float)
        if self.windows.shape!=(3,4)or self.translation.shape!=(3,)or not np.isfinite(np.r_[self.windows.ravel(),self.translation]).all():raise ValueError('Finite 3×4 windows and metre translation required')
        if np.any(np.diff(self.windows,axis=1)<=0):raise ValueError('Strictly ordered support/plateau bounds required')
        # Positive translation can reduce the determinant only on the falling
        # side of its support; negative translation uses the rising side.
        spans=np.where(self.translation>=0,self.windows[:,3]-self.windows[:,2],self.windows[:,1]-self.windows[:,0])
        self.determinant_lower_bound=float(1-np.sum(abs(self.translation)*1.5/spans))
        if not 0<self.minimum_jacobian<=1 or self.determinant_lower_bound<self.minimum_jacobian:raise ValueError('Translation exceeds the declared positive-Jacobian bound')

    def weights(self,points):
        p=np.asarray(points,float).reshape(-1,3);lo,a,b,hi=self.windows.T
        u=np.clip((p-lo)/(a-lo),0,1);v=np.clip((hi-p)/(hi-b),0,1)
        first=u*u*(3-2*u);second=v*v*(3-2*v);value=first*second
        derivative=(6*u*(1-u)/(a-lo))*second-first*(6*v*(1-v)/(hi-b))
        weight=np.prod(value,axis=1);gradient=np.empty_like(p)
        for k in range(3):gradient[:,k]=derivative[:,k]*np.prod(value[:,[j for j in range(3)if j!=k]],axis=1)
        return weight,gradient

    def map(self,points):
        p=np.asarray(points,float);weight,_=self.weights(p);return p+weight[:,None]*self.translation

    def jacobians(self,points):
        _,gradient=self.weights(points);j=np.broadcast_to(np.eye(3),(len(gradient),3,3)).copy();j+=self.translation[None,:,None]*gradient[:,None,:]
        if np.linalg.det(j).min(initial=1)<self.determinant_lower_bound-1e-12:raise ValueError('Analytic Jacobian bound failed')
        return j
