"""Compact positive-Jacobian translation of one convex protected core.

Unit outward planes describe the polytope (inside <= 0). A log-mean-exp
smooth maximum retains that entire polytope in the flat core. The core
advects with the translation; finite rank-one maps transport covectors.
The analytic certificate does not certify a discrete anatomical mesh.
"""
from dataclasses import dataclass
import numpy as np


def swept_face_separation(triangles, normals, offsets, translation, batch_size=256):
    """Common-plane lower bound for filled triangles versus C+[0,D]."""
    t=np.asarray(triangles,float);n=np.asarray(normals,float)
    b=np.asarray(offsets,float);d=np.asarray(translation,float)
    shifts=np.maximum(0.,n@d);values=[];planes=[]
    for start in range(0,len(t),batch_size):
        q=(t[start:start+batch_size]@n.T+b).min(1)-shifts;pick=q.argmax(1)
        values.extend(q[np.arange(len(q)),pick]);planes.extend(pick)
    return np.asarray(values),np.asarray(planes,dtype=np.int64)


@dataclass(frozen=True)
class ConvexCoreTransport:
    normals: np.ndarray
    offsets: np.ndarray
    translation: np.ndarray
    margin_m: float
    tau_m: float=2e-5
    steps: int=64

    def __post_init__(self):
        n=np.asarray(self.normals,float);b=np.asarray(self.offsets,float);d=np.asarray(self.translation,float)
        if n.ndim!=2 or n.shape[1]!=3 or b.shape!=(len(n),) or d.shape!=(3,) or len(n)<4:
            raise ValueError('Unit planes, matching offsets and one three-vector required')
        if not np.isfinite(np.r_[n.ravel(),b,d]).all() or not np.allclose(np.linalg.norm(n,axis=1),1.,rtol=0,atol=2e-14):
            raise ValueError('Finite unit outward planes required')
        if not np.isfinite(self.margin_m) or self.margin_m<=0 or not np.isfinite(self.tau_m) or self.tau_m<=0:
            raise ValueError('Positive finite falloff margin and smoothing scale required')
        if not isinstance(self.steps,int) or self.steps<1 or self.step_determinant_lower_bound<=0:
            raise ValueError('Step count must certify a strictly positive determinant')
        for key,value in [('normals',n),('offsets',b),('translation',d)]:
            value=value.copy();value.flags.writeable=False;object.__setattr__(self,key,value)

    @property
    def soft_core_expansion_upper_m(self):
        return self.tau_m*np.log(len(self.normals))

    @property
    def step_determinant_lower_bound(self):
        return 1.-np.linalg.norm(self.translation)*1.875/(self.steps*self.margin_m)

    def field(self,points):
        """Return weight, analytic gradient and smooth convex-core distance."""
        p=np.asarray(points,float);q=p@self.normals.T+self.offsets;maximum=q.max(1)
        e=np.exp((q-maximum[:,None])/self.tau_m);total=e.sum(1)
        distance=maximum+self.tau_m*np.log(total)-self.soft_core_expansion_upper_m
        u=np.clip(distance/self.margin_m,0.,1.)
        weight=np.clip(1.-u*u*u*(10.+u*(-15.+6.*u)),0.,1.)
        derivative=-30.*u*u*(1.-u)*(1.-u)/self.margin_m
        gradient=(e@self.normals)/total[:,None]*derivative[:,None]
        weight[distance<=0.]=1.;weight[distance>=self.margin_m]=0.
        gradient[(distance<=0.)|(distance>=self.margin_m)]=0.
        return weight,gradient,distance

    def apply(self,points,normals=None,batch_size=2048):
        """Apply the same field to all points, with no source identity mask.

        Initial flat-core points take their mathematically exact rigid result
        in one addition. Others use the declared finite maps. Covectors use
        the exact rank-one inverse transpose at every finite step. Numerically
        unchanged normal vectors retain their original representation.
        """
        p=np.asarray(points,float)
        if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite [N,3] positions required')
        if not isinstance(batch_size,int) or batch_size<1:raise ValueError('A positive batch size required')
        result=p.copy();original_normals=None if normals is None else np.asarray(normals,float)
        if original_normals is not None and (original_normals.shape!=p.shape or not np.isfinite(original_normals).all()):
            raise ValueError('One finite normal per coordinate required')
        mapped_normals=None if original_normals is None else original_normals.copy();delta=self.translation/self.steps
        stats={'points':len(p),'sweptSupportCandidates':0,'initialFlatCorePoints':0,'transitionPoints':0,
               'stepDeterminantMinimumObserved':1.,'stepDeterminantLowerBound':float(self.step_determinant_lower_bound),
               'steps':self.steps,'softCoreExpansionUpper_m':float(self.soft_core_expansion_upper_m)}
        shifts=np.maximum(0.,self.normals@self.translation)
        for start in range(0,len(p),batch_size):
            stop=min(start+batch_size,len(p));source=p[start:stop]
            plane=source@self.normals.T+self.offsets-shifts
            ids=np.flatnonzero(plane.max(1)<self.margin_m+self.soft_core_expansion_upper_m)+start
            if not len(ids):continue
            stats['sweptSupportCandidates']+=len(ids)
            _,_,distance=self.field(p[ids]);flat=distance<=0.
            result[ids[flat]]=p[ids[flat]]+self.translation;stats['initialFlatCorePoints']+=int(flat.sum())
            ids=ids[~flat]
            if not len(ids):continue
            x=p[ids].copy();cov=None if mapped_normals is None else original_normals[ids].copy()
            touched=np.zeros(len(ids),bool);gradient_touched=touched.copy()
            for step in range(self.steps):
                w,g,_=self.field(x-(step/self.steps)*self.translation);determinant=1.+g@delta
                if not np.isfinite(determinant).all() or determinant.min()<=0.:raise ValueError('A finite map lost its positive determinant')
                stats['stepDeterminantMinimumObserved']=min(stats['stepDeterminantMinimumObserved'],float(determinant.min()))
                touched|=w>0.;gradient_touched|=np.any(g!=0.,axis=1)
                if cov is not None:cov-=g*((cov@delta)/determinant)[:,None]
                x+=w[:,None]*delta
            result[ids]=x;stats['transitionPoints']+=int(touched.sum())
            if cov is not None:
                scale=np.linalg.norm(cov[gradient_touched],axis=1)
                if np.any(scale<=0.) or not np.isfinite(scale).all():raise ValueError('Transported normal is invalid')
                cov[gradient_touched]/=scale[:,None];mapped_normals[ids]=cov
        if not np.isfinite(result).all():raise ValueError('Transport returned nonfinite coordinates')
        return result,mapped_normals,stats
