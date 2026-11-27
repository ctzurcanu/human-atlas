"""Shared compact carriage about two measured source/recipient courses."""
import time
import numpy as np
from scipy.interpolate import PchipInterpolator

def carry_points(points,zgrid,source_centers,recipient_centers,*,steps=24,wall_seconds=4.):
    start=time.monotonic();q=np.array(points,float,copy=True);zgrid=np.asarray(zgrid,float);paths=[PchipInterpolator(zgrid,np.asarray(c))for c in source_centers];deltas=[PchipInterpolator(zgrid,np.asarray(t)-np.asarray(s))for t,s in zip(recipient_centers,source_centers)];J=np.broadcast_to(np.eye(3),(len(q),3,3)).copy()
    def evaluate(t,x,j):
        if time.monotonic()-start>wall_seconds:raise TimeoutError('Compact paired course point carriage')
        z=x[:,2];zc=np.clip(z,zgrid[0],zgrid[-1]);inside=(z>zgrid[0])&(z<zgrid[-1]);u=np.clip((z+.06)/(.06+zgrid[0]),0,1);fade=u**3*(10-15*u+6*u*u);fd=30*u*u*(1-u)**2/(.06+zgrid[0]);v=np.zeros_like(x);D=np.zeros_like(j)
        for path,delta in zip(paths,deltas):
            c0=path(zc);d=delta(zc)*fade[:,None];dd=delta.derivative()(zc)*inside[:,None]*fade[:,None]+delta(zc)*fd[:,None];c=c0+t*d;cd=path.derivative()(zc)*inside[:,None]+t*dd;r=x[:,:2]-c;dist=np.linalg.norm(r,axis=1)
            u=np.clip((.013-dist)/.01,0,1);w=u**3*(10-15*u+6*u*u);dw=-30*u*u*(1-u)**2/.01;grad=np.zeros_like(x);grad[:,:2]=np.divide(dw[:,None]*r,dist[:,None],out=np.zeros_like(r),where=dist[:,None]>1e-12);grad[:,2]=-np.einsum('ni,ni->n',grad[:,:2],cd)
            full=np.c_[d,np.zeros(len(x))];v+=w[:,None]*full;D+=full[:,:,None]*grad[:,None,:];D[:,:2,2]+=w[:,None]*dd
        return v,D@j
    dt=1/steps
    for k in range(steps):
        t=k*dt;v1,j1=evaluate(t,q,J);v2,j2=evaluate(t+dt/2,q+dt*v1/2,J+dt*j1/2);v3,j3=evaluate(t+dt/2,q+dt*v2/2,J+dt*j2/2);v4,j4=evaluate(t+dt,q+dt*v3,J+dt*j3);q+=dt*(v1+2*v2+2*v3+v4)/6;J+=dt*(j1+2*j2+2*j3+j4)/6
    if not np.isfinite(q).all()or np.linalg.det(J).min(initial=1)<=0:raise ValueError('Nonpositive mapped point Jacobian')
    return q,J
