"""Extend observed proper profile poses within an explicit anatomical group.

These are source-indexed registration observations, not inferred interfaces
or collision acceptance. The grouping prevents neighbouring digits/sides
from contaminating each other's pose fields.
"""
import numpy as np
from scipy.spatial import cKDTree

from .dual_quaternion_carriage import blend_proper_poses


def rigid_profile_observations(source,target):
    a,b=np.asarray(source,float),np.asarray(target,float)
    if a.shape!=b.shape or a.ndim!=3 or a.shape[2]!=3 or a.shape[1]<3 or not len(a):
        raise ValueError('Corresponding nonempty complete source profiles required')
    if not np.isfinite(a).all()or not np.isfinite(b).all():raise ValueError('Finite profile observations required')
    c,d=a.mean(1),b.mean(1);x,y=a-c[:,None],b-d[:,None]
    u,s,vt=np.linalg.svd(np.einsum('nki,nkj->nij',x,y))
    if np.any(s[:,1]<1e-20):raise ValueError('A collapsed profile cannot determine a proper pose')
    v=vt.transpose(0,2,1);sign=np.ones((len(a),3));sign[:,-1]=np.linalg.det(v@u.transpose(0,2,1))
    R=(v*sign[:,None,:])@u.transpose(0,2,1);t=d-np.einsum('nij,nj->ni',R,c)
    fit=np.einsum('nij,nkj->nki',R,a)+t[:,None]
    return c,R,t,np.linalg.norm(fit-b,axis=2).max(1)


class ObservedProfileCarrier:
    def __init__(self,centres,rotations,translations,groups):
        self.centres=np.asarray(centres,float);self.R=np.asarray(rotations,float);self.t=np.asarray(translations,float);self.groups=np.asarray(groups,int)
        n=len(self.centres)
        if self.centres.shape!=(n,3)or self.R.shape!=(n,3,3)or self.t.shape!=(n,3)or self.groups.shape!=(n,)or not n:
            raise ValueError('One explicit proper pose and anatomical group per observation required')
        if not np.isfinite(np.r_[self.centres.ravel(),self.R.ravel(),self.t.ravel()]).all():raise ValueError('Finite poses required')
        if not np.allclose(self.R.transpose(0,2,1)@self.R,np.eye(3),atol=1e-12,rtol=0)or np.any(np.linalg.det(self.R)<1-1e-12):
            raise ValueError('Reflection/shear is not a proper source profile pose')

    def poses(self,points,groups,*,neighbours=12,maximum_source_distance_m=.02):
        p=np.asarray(points,float);labels=np.asarray(groups,int)
        if p.ndim!=2 or p.shape[1]!=3 or labels.shape!=(len(p),)or not np.isfinite(p).all():raise ValueError('Finite owned source queries required')
        if not 1<=neighbours<=64 or not 0<maximum_source_distance_m<=.05:raise ValueError('Finite explicit neighbourhood limits required')
        q=np.empty_like(p);R=np.empty((len(p),3,3));proof=[]
        for label in np.unique(labels):
            queries=np.flatnonzero(labels==label);obs=np.flatnonzero(self.groups==label)
            if not len(obs):raise ValueError('No observations for the requested anatomical group')
            k=min(neighbours,len(obs));distance,index=cKDTree(self.centres[obs]).query(p[queries],k=k)
            distance,index=np.asarray(distance).reshape(len(queries),k),np.asarray(index).reshape(len(queries),k)
            if np.any(distance[:,0]>maximum_source_distance_m):raise ValueError('Source query outside its declared observation domain')
            w=1/np.maximum(distance,1e-12)**2;exact=distance==0
            if exact.any():w[exact.any(1)]=exact[exact.any(1)].astype(float)
            w/=w.sum(1)[:,None];weights=np.zeros((len(queries),len(obs)))
            weights[np.arange(len(queries))[:,None],index]=w
            q[queries],R[queries]=blend_proper_poses(p[queries],weights,self.R[obs],self.t[obs])
            proof.append(dict(group=int(label),queries=len(queries),observations=len(obs),nearestSourceDistanceMaximum_m=float(distance[:,0].max(initial=0))))
        return q,R,dict(groups=proof,anatomicalGroupsKeptSeparate=True,properRotationsOnly=True,anatomicalAcceptance=False)
