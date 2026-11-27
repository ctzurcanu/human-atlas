"""Bounded translation fit for explicitly declared rigid source profiles.

Surface samples and oriented closest facets propose clearance constraints;
they are not continuous intersection evidence. Callers must independently
check every emitted triangle, fixed owner, source interface and packed field.
"""
import time

import numpy as np
from scipy.optimize import minimize

from ..open3d_backend import SurfaceIndex


def fit_profile_translations(reference,base,weights,faces,obstacle,obstacle_normals,*,
                             maximum_motion_m,minimum_clearance_m,
                             maximum_facet_rotation_degrees,
                             wall_seconds=12.,maximum_iterations=80):
    tick=time.monotonic();old=np.asarray(reference,float);base=np.asarray(base,float);w=np.asarray(weights,float);f=np.asarray(faces,int)
    normals=np.asarray(obstacle_normals,float);k=w.shape[1];scale=.001
    if old.shape!=base.shape or w.shape[0]!=len(old)or np.any(w<0)or np.any(w.sum(1)>1):raise ValueError('Corresponding explicit rigid-profile weights required')
    if not 0<wall_seconds<=20 or not 1<=maximum_iterations<=100:raise ValueError('Finite fitting stage required')
    index=SurfaceIndex(obstacle);n=6
    bary=np.array([[i/n,j/n,1-(i+j)/n]for i in range(n+1)for j in range(n+1-i)])
    sw=np.einsum('bi,tik->tbk',bary,w[f]).reshape(-1,k)
    old_samples=np.einsum('bi,tij->tbj',bary,old[f]).reshape(-1,3)
    sample_base=np.einsum('bi,tij->tbj',bary,base[f]).reshape(-1,3)
    active=sw.sum(1)>0;sw=sw[active];old_samples=old_samples[active];sample_base=sample_base[active]
    def signed(q):
        a,_=index.closest_points(q)
        normal=normals[a['face_ids']]
        return np.einsum('ij,ij->i',q-a['points'],normal),normal
    initial_signed,_=signed(old_samples)
    required=np.minimum(minimum_clearance_m,initial_signed)
    tri=old[f];ea=tri[:,1]-tri[:,0];eb=tri[:,2]-tri[:,0];original=np.cross(ea,eb);original/=np.linalg.norm(original,axis=1)[:,None]
    dw1=w[f[:,1]]-w[f[:,0]];dw2=w[f[:,2]]-w[f[:,0]]
    live=np.flatnonzero(w.sum(1)>0);counter=0;cached=None
    def evaluate(x):
        nonlocal counter,cached
        if cached is not None and np.array_equal(cached[0],x):return cached[1:]
        if time.monotonic()-tick>wall_seconds:raise TimeoutError('Profile fitting native stage deadline')
        counter+=1;v=base+w@x.reshape(k,3)*scale;q=sample_base+sw@x.reshape(k,3)*scale
        distance,normal=signed(q);clearance=(distance-required)/scale
        jclear=(sw[:,:,None]*normal[:,None,:]).reshape(len(sw),3*k)
        tri=v[f];e1=tri[:,1]-tri[:,0];e2=tri[:,2]-tri[:,0];a=np.cross(e1,e2);length=np.linalg.norm(a,axis=1);unit=a/length[:,None]
        cosine=np.einsum('ij,ij->i',unit,original);rotation=cosine-np.cos(np.radians(maximum_facet_rotation_degrees))
        tangent=(original-cosine[:,None]*unit)/length[:,None]
        jrotation=np.empty((len(f),k,3))
        for axis in range(3):
            basis=np.eye(3)[axis]
            derivative=np.cross(dw1[:,:,None]*basis[None,None,:],e2[:,None,:])+np.cross(e1[:,None,:],dw2[:,:,None]*basis[None,None,:])
            jrotation[:,:,axis]=np.einsum('tki,ti->tk',derivative,tangent)*scale
        delta=(v[live]-old[live])/scale;radius=np.linalg.norm(delta,axis=1);motion=maximum_motion_m/scale-radius
        direction=np.divide(delta,radius[:,None],out=np.zeros_like(delta),where=radius[:,None]>0)
        jmotion=-(w[live,:,None]*direction[:,None,:]).reshape(len(live),3*k)
        objective=float(np.mean(np.sum(delta*delta,axis=1)))
        grad=(2*np.einsum('ni,nk->ki',delta,w[live])/len(live)).ravel()
        constraints=np.r_[clearance,rotation,motion];jac=np.vstack([jclear,jrotation.reshape(len(f),3*k),jmotion])
        cached=(x.copy(),objective,grad,constraints,jac,v)
        return cached[1:]
    start=np.zeros(k*3)
    result=minimize(lambda x:evaluate(x)[0],start,jac=lambda x:evaluate(x)[1],method='SLSQP',
        constraints=[dict(type='ineq',fun=lambda x:evaluate(x)[2],jac=lambda x:evaluate(x)[3])],
        bounds=[(-maximum_motion_m/scale,maximum_motion_m/scale)]*(k*3),
        options=dict(maxiter=maximum_iterations,ftol=1e-10,disp=False))
    objective,gradient,constraints,jac,positions=evaluate(result.x)
    proof=dict(complete=True,solverSuccess=bool(result.success),solverMessage=str(result.message),iterations=int(result.nit),evaluations=counter,
        minimumProposalConstraint=float(constraints.min()),profileTranslations_mm=result.x.reshape(k,3).tolist(),
        sampledPoints=len(sw),samplesPerTriangle=len(bary),proposalOnlyNotIntersectionProof=True,
        maximumMotion_m=float(np.linalg.norm(positions-old,axis=1).max()),objective=objective,elapsedSeconds=time.monotonic()-tick)
    if not result.success or constraints.min() < -1e-7:return None,proof
    return positions,proof
