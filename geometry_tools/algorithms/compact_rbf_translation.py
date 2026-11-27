"""Small shared C2 corrections with an explicit global injectivity bound.

The whole assembly receives one sum of compact translations. Literal shared
coordinates stay shared. Continuous injectivity does not certify triangulated
surfaces, tissue contact or anatomy; independently check the emitted arrays.
"""
import numpy as np
from .compact_translation import CompactTranslation


class CompactRBFTranslation:
    def __init__(self,centres,displacements,radius_m,frame,*,maximum_lipschitz=.1):
        p=np.asarray(centres,float);d=np.asarray(displacements,float)
        if p.ndim!=2 or p.shape[1]!=3 or d.shape!=p.shape or not 1<=len(p)<=256:
            raise ValueError('One to 256 explicit paired controls required')
        if not np.isfinite(p).all()or not np.isfinite(d).all()or not 0<radius_m or not 0<maximum_lipschitz<1:
            raise ValueError('Finite controls, positive support and injective bound required')
        basis=[CompactTranslation(c,np.full(3,radius_m),np.zeros(3),frame)for c in p]
        W=np.column_stack([f.weights_and_gradients(p)[0]for f in basis])
        if np.linalg.cond(W)>1e10:raise ValueError('Ill-conditioned or duplicate correction controls')
        coefficients=np.linalg.solve(W,d)
        self.fields=[CompactTranslation(c,np.full(3,radius_m),q,frame)for c,q in zip(p,coefficients)]
        self.lipschitz_bound=sum(f.lipschitz_bound for f in self.fields)
        if self.lipschitz_bound>=maximum_lipschitz:raise ValueError('Shared correction exceeds its global distortion bound')
        self.frame=frame
        if not np.allclose(self.map(p),p+d,rtol=0,atol=1e-14):raise ValueError('Correction controls did not replay')

    def map(self,points):
        p=np.asarray(points,float);delta=np.zeros_like(p)
        for f in self.fields:delta+=f.weights_and_gradients(p)[0][:,None]*f.displacement
        return p+delta

    def jacobians(self,points):
        p=np.asarray(points,float);J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy()
        for f in self.fields:J+=f.displacement[None,:,None]*f.weights_and_gradients(p)[1][:,None,:]
        return J

    def transform_normals(self,points,normals):
        n=np.asarray(normals,float);q=np.linalg.solve(self.jacobians(points).transpose(0,2,1),n[...,None])[...,0]
        length=np.linalg.norm(q,axis=1)
        if np.any(length==0):raise ValueError('Zero normals cannot be transported')
        return q/length[:,None]

    def evidence(self):
        return dict(frame=self.frame,fields=[f.evidence()for f in self.fields],
            sumDisplacementLipschitzBound=self.lipschitz_bound,
            continuousMinimumStretchBound=1-self.lipschitz_bound,
            continuousMapGloballyInjective=True,discreteGeometryAndContactsVerified=False)
