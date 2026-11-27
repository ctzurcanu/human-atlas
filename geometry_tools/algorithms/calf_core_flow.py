"""Common ambient XZ rigid-core velocity flow with complete Y derivatives.

No mesh, anatomy owner or cage interpolation is used. The exact smooth ODE flow
has positive determinant; finite RK4 positions/Jacobians require their own audit.
"""
from __future__ import annotations
import numpy as np
from scipy.interpolate import PchipInterpolator


def smoothstep5(x):
    x=np.asarray(x,dtype=float)
    return x*x*x*(10+x*(-15+6*x))


def smoothstep5_derivative(x):
    x=np.asarray(x,dtype=float)
    return 30*x*x*(1-x)*(1-x)


class QuinticMajorant:
    """C2 absolute support majorant of endpoint chords at source vertex-Y knots.

    Between consecutive source vertex Y values, all active mesh edges are
    affine. Their maximum support is convex, so its endpoint chord dominates.
    Replacing linear s with h5(s) changes that chord by at most
    max|s-h5(s)| |H_R-H_L| < .15 |H_R-H_L|. A constant per-normal epsilon
    preserves C2 joins and dominates every interval, including short ones.
    """
    def __init__(self,ys,values,epsilon=None):
        self.ys=np.asarray(ys,float);self.values=np.asarray(values,float)
        if self.ys.ndim!=1 or len(self.ys)<2 or np.any(np.diff(self.ys)<=0):raise ValueError('Strictly ordered support knots required')
        if self.values.ndim!=2 or len(self.values)!=len(self.ys)or not np.isfinite(self.values).all():raise ValueError('Finite support matrix required')
        required=.15*np.max(np.abs(np.diff(self.values,axis=0)),axis=0)
        self.epsilon=required if epsilon is None else np.asarray(epsilon,float)
        if self.epsilon.shape!=required.shape or np.any(self.epsilon<required-1e-15):raise ValueError('Support epsilon does not dominate the analytic chord error bound')

    def evaluate(self,y):
        y=np.asarray(y,float);k=np.clip(np.searchsorted(self.ys,y,side='right')-1,0,len(self.ys)-2);width=self.ys[k+1]-self.ys[k];s=np.clip((y-self.ys[k])/width,0,1);delta=self.values[k+1]-self.values[k]
        value=self.values[k]+smoothstep5(s)[...,None]*delta+self.epsilon
        derivative=smoothstep5_derivative(s)[...,None]*delta/width[...,None]
        derivative=np.where(((y>self.ys[0])&(y<self.ys[-1]))[...,None],derivative,0.)
        return value,derivative


def longitudinal_strength(y,core=(.14,.36),support=(.10,.40)):
    y=np.asarray(y,float);left=(y-support[0])/(core[0]-support[0]);right=(support[1]-y)/(support[1]-core[1]);s=np.where(y<core[0],left,np.where(y>core[1],right,1.));s=np.clip(s,0,1);sign=np.where(y<core[0],1/(core[0]-support[0]),np.where(y>core[1],-1/(support[1]-core[1]),0.))
    return smoothstep5(s),smoothstep5_derivative(s)*sign


class RigidCore:
    def __init__(self,spec,normals):
        self.normals=np.asarray(normals,float);ys=np.asarray(spec['stationY_m'],float)
        if self.normals.ndim!=2 or self.normals.shape[1]!=2 or not np.allclose(np.linalg.norm(self.normals,axis=1),1,atol=1e-12):raise ValueError('Unit XZ support normals required')
        if np.any(np.diff(ys)<=0):raise ValueError('Strictly ordered guide stations required')
        self.c=PchipInterpolator(ys,np.asarray(spec['sourceCentresXZ_m'],float),axis=0)
        self.d=PchipInterpolator(ys,np.asarray(spec['targetDeltasXZ_m'],float),axis=0)
        self.theta=PchipInterpolator(ys,np.asarray(spec['rotationRadians'],float))
        self.shape=QuinticMajorant(spec['sourceSupportY_m'],spec['sourceSupportOffsetsWorld_m'],spec['sourceSupportEpsilon_m'])
        self.margin=float(spec['supportMargin_m'])
        if not self.margin>0 or not np.isfinite(self.margin):raise ValueError('A measured positive support margin is required')

    def velocity(self,points,t):
        p=np.asarray(points,float);v=np.zeros_like(p);A=np.zeros((len(p),3,3));sigma,sy=longitudinal_strength(p[:,1]);ids=np.flatnonzero(sigma>0)
        if not len(ids):return v,A
        y=p[ids,1];s=sigma[ids];s_y=sy[ids];q=float(smoothstep5(t));dq=float(smoothstep5_derivative(t));c=self.c(y);cy=self.c(y,1);d=self.d(y);dy=self.d(y,1);theta=self.theta(y);thetay=self.theta(y,1);H,Hy=self.shape.evaluate(y)
        centre=c+(s*q)[:,None]*d;centre_y=cy+q*(s_y[:,None]*d+s[:,None]*dy);angle=s*q*theta;angle_y=q*(s_y*theta+s*thetay);co=np.cos(angle);sn=np.sin(angle);nx,nz=self.normals.T
        nt=np.stack([co[:,None]*nx-sn[:,None]*nz,sn[:,None]*nx+co[:,None]*nz],axis=2);nty=np.stack([-nt[:,:,1],nt[:,:,0]],axis=2)*angle_y[:,None,None];relative=p[ids][:,[0,2]]-centre
        h=H-c@self.normals.T;hy=Hy-cy@self.normals.T;offset=np.einsum('ni,nmi->nm',relative,nt)-h
        active=(offset<self.margin).all(1);ids=ids[active]
        if not len(ids):return v,A
        s,s_y,d,dy,theta,thetay,centre_y,relative,nt,nty,offset,hy=[x[active]for x in [s,s_y,d,dy,theta,thetay,centre_y,relative,nt,nty,offset,hy]]
        u=np.clip(offset/self.margin,0,1);f=(1-u)**4*(4*u+1);fp=np.where((offset>0)&(offset<self.margin),-20*u*(1-u)**3,0.)
        # Prefix/suffix products avoid division by tiny/zero factors near support.
        prefix=np.c_[np.ones(len(f)),np.cumprod(f,axis=1)[:,:-1]];suffix=np.c_[np.cumprod(f[:,::-1],axis=1)[:,::-1][:,1:],np.ones(len(f))];wf=prefix*suffix*fp;weight=np.prod(f,axis=1)
        gradxz=np.einsum('nm,nmi->ni',wf,nt)/self.margin;offset_y=np.einsum('ni,nmi->nm',relative,nty)-np.einsum('ni,nmi->nm',centre_y,nt)-hy;weight_y=np.sum(wf*offset_y,axis=1)/self.margin
        omega=s*dq*theta;omega_y=dq*(s_y*theta+s*thetay);Krelative=np.c_[-relative[:,1],relative[:,0]];Kcentre_y=np.c_[-centre_y[:,1],centre_y[:,0]];rigid=(s*dq)[:,None]*d+omega[:,None]*Krelative;rigid_y=dq*(s_y[:,None]*d+s[:,None]*dy)+omega_y[:,None]*Krelative-omega[:,None]*Kcentre_y
        v[ids,0]=weight*rigid[:,0];v[ids,2]=weight*rigid[:,1];local=rigid[:,:,None]*gradxz[:,None,:];local[:,0,1]-=weight*omega;local[:,1,0]+=weight*omega;vy=rigid*weight_y[:,None]+weight[:,None]*rigid_y
        A[ids,0,0]=local[:,0,0];A[ids,0,2]=local[:,0,1];A[ids,2,0]=local[:,1,0];A[ids,2,2]=local[:,1,1];A[ids,0,1]=vy[:,0];A[ids,2,1]=vy[:,1]
        return v,A


class CalfCoreFlow:
    def __init__(self,spec):
        self.cores=[RigidCore(field,spec['fixedAngularNormalsXZ'])for field in spec['fields']]

    def velocity(self,points,t):
        if not 0<=t<=1:raise ValueError('Flow time must be in [0,1]')
        p=np.asarray(points,float)
        if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite Nx3 source points required')
        v=np.zeros_like(p);A=np.zeros((len(p),3,3))
        for core in self.cores:
            dv,da=core.velocity(p,t);v+=dv;A+=da
        return v,A

    def integrate(self,points,*,steps,normals=None):
        """Direct common RK4 positions and full Dphi; max8192 points per batch."""
        initial=np.asarray(points,float);p=initial.copy()
        if len(p)>8192 or not isinstance(steps,int)or not 1<=steps<=4096:raise ValueError('Explicit1..4096 steps and <=8192 source points per batch required')
        J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();logdet=np.zeros(len(p));dt=1./steps
        def rhs(x,j,t):
            v,A=self.velocity(x,t);return v,A@j,np.trace(A,axis1=1,axis2=2)
        for k in range(steps):
            t=k*dt;p1,j1,l1=rhs(p,J,t);p2,j2,l2=rhs(p+.5*dt*p1,J+.5*dt*j1,t+.5*dt);p3,j3,l3=rhs(p+.5*dt*p2,J+.5*dt*j2,t+.5*dt);p4,j4,l4=rhs(p+dt*p3,J+dt*j3,t+dt)
            p+=dt*(p1+2*p2+2*p3+p4)/6;J+=dt*(j1+2*j2+2*j3+j4)/6;logdet+=dt*(l1+2*l2+2*l3+l4)/6
        p[:,1]=initial[:,1];answer=dict(points=p,jacobian=J,integrated_log_determinant=logdet,finite_jacobian_determinant=np.linalg.det(J),YFixedExactly=np.array_equal(p[:,1],initial[:,1]),steps=steps,exactFlowPositiveDeterminantTheoremIsFiniteIntegrationAcceptance=False)
        if normals is not None:
            n=np.asarray(normals,float)
            if n.shape!=p.shape:raise ValueError('One normal per source point required')
            mapped=np.linalg.solve(J.transpose(0,2,1),n[...,None])[...,0];length=np.linalg.norm(mapped,axis=1);mapped=np.divide(mapped,length[:,None],out=np.zeros_like(mapped),where=length[:,None]>0);answer['normals']=mapped
        return answer
