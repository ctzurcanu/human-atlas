"""Smooth compact translation of explicitly declared neighbouring assemblies."""
import time
import numpy as np
def carry_points(points,centres,deltas,*,core_m=.02,outer_m=.045,steps=24,wall_seconds=6.,upper_station=None):
    start=time.monotonic();q=np.array(points,float,copy=True);centres=np.asarray(centres,float);deltas=np.asarray(deltas,float);J=np.broadcast_to(np.eye(3),(len(q),3,3)).copy()
    if not 0<core_m<outer_m or centres.shape!=deltas.shape or centres.shape[1:]!=(3,):raise ValueError('Declared compact translation controls required')
    if upper_station is not None:
        axis,full,stop=upper_station
        if axis not in (0,1,2) or not full<stop:raise ValueError('Increasing station interval required')
    def rhs(t,x,j):
        if time.monotonic()-start>wall_seconds:raise TimeoutError('Compact translation flow deadline')
        v=np.zeros_like(x);D=np.zeros_like(j);weight=np.zeros(len(x));gradient=np.zeros_like(x)
        for c,d in zip(centres,deltas):
            r=x-(c+t*d);distance=np.linalg.norm(r,axis=1);u=np.clip((outer_m-distance)/(outer_m-core_m),0,1);w=u**3*(10-15*u+6*u*u);dw=-30*u*u*(1-u)**2/(outer_m-core_m);g=np.divide(dw[:,None]*r,distance[:,None],out=np.zeros_like(r),where=distance[:,None]>1e-12);v+=w[:,None]*d;D+=d[None,:,None]*g[:,None,:];weight+=w;gradient+=g
        den=np.maximum(weight,1);D=D/den[:,None,None]-v[:,:,None]*np.where(weight>1,1,0)[:,None,None]*gradient[:,None,:]/den[:,None,None]**2;v/=den[:,None]
        if upper_station is not None:
            u=np.clip((stop-x[:,axis])/(stop-full),0,1);w=u**3*(10-15*u+6*u*u);g=np.zeros_like(x);g[:,axis]=-30*u*u*(1-u)**2/(stop-full)
            D=w[:,None,None]*D+v[:,:,None]*g[:,None,:];v*=w[:,None]
        return v,D@j
    dt=1/steps
    for k in range(steps):
        t=k*dt;v1,j1=rhs(t,q,J);v2,j2=rhs(t+dt/2,q+dt*v1/2,J+dt*j1/2);v3,j3=rhs(t+dt/2,q+dt*v2/2,J+dt*j2/2);v4,j4=rhs(t+dt,q+dt*v3,J+dt*j3);q+=dt*(v1+2*v2+2*v3+v4)/6;J+=dt*(j1+2*j2+2*j3+j4)/6
    if np.linalg.det(J).min(initial=1)<=0:raise ValueError('Translation flow Jacobian failed')
    return q,J
