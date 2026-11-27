"""Explicit three-dimensional source-patch proposal with actual metric bounds.

All fixed controls remain literal. This numerical construction establishes
neither anatomical interfaces nor acceptance; whole-mesh checks stay external.
"""
import time
import numpy as np
from scipy.optimize import minimize


def separate_patch(vertices, faces, first, second, free, direction, *,
                   gap_m, maximum_motion_m, maximum_strain, wall_seconds):
    start=time.monotonic();p=np.asarray(vertices,float);f=np.asarray(faces,int)
    free=np.asarray(sorted(set(map(int,free))),int);n=np.asarray(direction,float)
    if not len(free) or not np.isfinite(p).all() or abs(np.linalg.norm(n)-1)>1e-10:
        raise ValueError('Explicit patch, finite positions and unit direction required')
    if min(gap_m,maximum_motion_m,maximum_strain,wall_seconds)<=0 or maximum_strain>=1:
        raise ValueError('Positive bounds/deadline and strain less than one required')
    scale=1e6;x=(p-p[np.r_[first,second]].mean(0))*scale
    index=np.full(len(p),-1,int);index[free]=np.arange(len(free));K=len(free)
    touched=np.flatnonzero(np.any(np.isin(f,free),axis=1));tri=x[f[touched]]
    E=np.swapaxes(tri[:,1:]-tri[:,:1],1,2)
    U,s,_=np.linalg.svd(E,full_matrices=False)
    if np.any(s[:,-1]<=1e-10):raise ValueError('Nondegenerate source facets required')
    G=np.linalg.pinv(E)@U
    D=np.zeros((len(touched),2,K))
    for corner,factor in [(0,-G.sum(1)),(1,G[:,0,:]),(2,G[:,1,:])]:
        for j,raw in enumerate(f[touched,corner]):
            if index[raw]>=0:D[j,:,index[raw]]+=factor[j]
    C=[];rhs=[]
    for a in first:
        for b in second:
            row=np.zeros((K,3))
            if index[a]>=0:row[index[a]]+=n
            if index[b]>=0:row[index[b]]-=n
            C.append(row.ravel());rhs.append(gap_m*scale-(x[a]-x[b])@n)
    C=np.asarray(C);rhs=np.asarray(rhs)
    def deadline(_=None):
        if time.monotonic()-start>wall_seconds:raise TimeoutError('Metric patch native deadline')
    def objective(u):deadline();return .5*float(u@u)
    def tensors(u):
        F=U+np.einsum('fjk,kc->fcj',D,u.reshape(K,3))
        eigen,vectors=np.linalg.eigh(np.swapaxes(F,1,2)@F)
        return F,eigen,vectors
    def metrics(u):
        _,e,_=tensors(u)
        return np.r_[(e[:,0]-(1-maximum_strain)**2),((1+maximum_strain)**2-e[:,1])]
    def metrics_jac(u):
        F,_,v=tensors(u);rows=[]
        for eig,sign in [(0,1),(1,-1)]:
            vv=v[:,:,eig];physical=np.einsum('fcj,fj->fc',F,vv)
            coefficient=np.einsum('fjk,fj->fk',D,vv)
            rows.append((sign*2*coefficient[:,:,None]*physical[:,None,:]).reshape(len(F),-1))
        return np.concatenate(rows)
    def motion(u):return (maximum_motion_m*scale)**2-np.sum(u.reshape(K,3)**2,axis=1)
    def motion_jac(u):
        j=np.zeros((K,K,3));j[np.arange(K),np.arange(K)]=-2*u.reshape(K,3)
        return j.reshape(K,-1)
    r=minimize(objective,np.zeros(K*3),jac=lambda u:u,method='SLSQP',
        bounds=[(-maximum_motion_m*scale,maximum_motion_m*scale)]*(K*3),
        constraints=[dict(type='ineq',fun=lambda u:C@u-rhs,jac=lambda u:C),
            dict(type='ineq',fun=metrics,jac=metrics_jac),
            dict(type='ineq',fun=motion,jac=motion_jac)],
        callback=deadline,options=dict(maxiter=150,ftol=1e-12,disp=False))
    deadline();q=p.copy();q[free]+=r.x.reshape(K,3)/scale
    _,e,_=tensors(r.x);strain=float(np.abs(np.sqrt(e)-1).max())
    gap=float(((q[first]-p[first[0]])@n).min()-((q[second]-p[first[0]])@n).max())
    motion=float(np.linalg.norm(q-p,axis=1).max())
    proof=dict(success=bool(r.success),message=str(r.message),iterations=int(r.nit),
        maximumIntrinsicFacetStrain=strain,projectedGap_m=gap,maximumMotion_m=motion,
        elapsedSeconds=time.monotonic()-start,anatomicalAcceptance=False)
    if not r.success or strain>maximum_strain*(1+1e-7) or gap<gap_m-1e-12 or motion>maximum_motion_m*(1+1e-7):
        raise ValueError('Inadmissible actual-metric patch: '+str(proof))
    return q,proof
