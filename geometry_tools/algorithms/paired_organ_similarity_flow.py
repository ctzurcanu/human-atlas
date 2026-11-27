"""Compact paired-organ registration with a monotone midsagittal transition.

Similarity cores preserve local calibre; the smooth halo leaves distant roots
literal. The caller selects the dependency closure and independently protects
other assemblies. This is a geometric carriage, not anatomical acceptance.
"""
import time
import numpy as np
from .feature_normals import FeatureNormals,classified_oriented_distance

def carry_points(points,source_mid,target_mid,scale,source_halfspan,target_halfspan,index,*,steps=24,core_m=.006,outer_m=.035,midline_m=.0025,wall_seconds=12.):
    started=time.monotonic();q=np.array(points,float,copy=True);c0=np.asarray(source_mid,float);c1=np.asarray(target_mid,float);s=float(scale);k=float(target_halfspan-s*source_halfspan);offset=c1-s*c0;normal=FeatureNormals(index.mesh)
    if q.ndim!=2 or q.shape[1]!=3 or not np.isfinite(q).all()or not 0<core_m<outer_m or s<=0 or k<0 or steps<1 or wall_seconds<=0:raise ValueError('Finite monotone paired similarity and explicit finite limits required')
    J=np.broadcast_to(np.eye(3),(len(q),3,3)).copy();I=np.eye(3)
    def rhs(t,x,j):
        if time.monotonic()-started>wall_seconds:raise TimeoutError('Paired-organ local flow deadline')
        a=1+t*(s-1);material=(x-t*offset)/a;mx=material[:,0].copy();lo=(x[:,0]-t*offset[0]-t*k)/a;hi=(x[:,0]-t*offset[0]+t*k)/a
        for _ in range(24):
            th=np.tanh((mx-c0[0])/midline_m);f=a*mx+t*offset[0]+t*k*th-x[:,0];derivative=a+t*k/midline_m*(1-th*th);done=abs(f)<1e-13;take=f>0;hi=np.where(take&~done,mx,hi);lo=np.where(~take&~done,mx,lo);nx=mx-f/derivative;outside=(nx<lo)|(nx>hi);mx=np.where(done,mx,np.where(outside,(lo+hi)/2,nx))
            if done.all():break
        material[:,0]=mx;th=np.tanh((mx-c0[0])/midline_m);assert np.max(abs(a*mx+t*offset[0]+t*k*th-x[:,0]),initial=0)<1e-8
        hit,_=index.closest_points(material);nv,rim=normal.at(hit);signed,known=classified_oriented_distance(material,hit,nv,rim);assert known.all();distance=np.maximum(signed,0)
        u=np.clip((outer_m-distance)/(outer_m-core_m),0,1);w=u**3*(10-15*u+6*u*u);dw=-30*u*u*(1-u)**2/(outer_m-core_m);direction=np.divide(material-hit['points'],hit['distance_m'][:,None],out=np.zeros_like(material),where=hit['distance_m'][:,None]>1e-12);gradient=dw[:,None]*direction;gradient[signed<0]=0
        delta=(s-1)*material+offset;delta[:,0]+=k*th;diagonal=np.full_like(x,s);diagonal[:,0]+=k/midline_m*(1-th*th);inverse=1/(1+t*(diagonal-1));gradient*=inverse
        D=delta[:,:,None]*gradient[:,None,:];D+=w[:,None,None]*(diagonal-1)[:,:,None]*inverse[:,:,None]*I
        return w[:,None]*delta,D@j
    dt=1/steps
    for step in range(steps):
        t=step*dt;v1,j1=rhs(t,q,J);v2,j2=rhs(t+dt/2,q+dt*v1/2,J+dt*j1/2);v3,j3=rhs(t+dt/2,q+dt*v2/2,J+dt*j2/2);v4,j4=rhs(t+dt,q+dt*v3,J+dt*j3);q+=dt*(v1+2*v2+2*v3+v4)/6;J+=dt*(j1+2*j2+2*j3+j4)/6
        if np.linalg.det(J).min(initial=1)<=0 or not np.isfinite(q).all():raise ValueError('Nonpositive paired-flow Jacobian')
    return q,J
