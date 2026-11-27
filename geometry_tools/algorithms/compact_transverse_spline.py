"""A compact transverse tube-course map with linear source-vertex controls.

The Y coordinate is fixed. Complete transverse sections inside the X/Z core
translate rigidly. Nonnegative B-spline bases bound displacement by the control
convex hull. Continuous positivity never substitutes for discrete/anatomy checks.
"""
import numpy as np
from scipy.interpolate import BSpline,PchipInterpolator


class CompactTransverseSpline:
    def __init__(self,stations,knots,*,degree=3,core_xz_m,support_xz_m):
        values=np.asarray([[s['y_m'],s['x_m'],s['z_m']]for s in stations],float)
        self.knots=np.asarray(knots,float);self.degree=int(degree)
        self.core=np.asarray(core_xz_m,float);self.support=np.asarray(support_xz_m,float)
        if (values.ndim!=2 or values.shape[1]!=3 or len(values)<2 or not np.isfinite(values).all()
            or np.any(np.diff(values[:,0])<=0)or self.degree!=3
            or self.core.shape!=(2,)or self.support.shape!=(2,)
            or not np.isfinite(np.r_[self.core,self.support,self.knots]).all()
            or np.any(self.core<0)or np.any(self.support<=self.core)
            or self.knots.ndim!=1 or len(self.knots)<8 or np.any(np.diff(self.knots)<0)):
            raise ValueError('Ordered finite stations, cubic knots and positive transverse transition widths required')
        self.count=len(self.knots)-self.degree-1
        if self.count<6:raise ValueError('At least six cubic controls are needed for fixed endpoint jets')
        self.domain=self.knots[[self.degree,-self.degree-1]]
        if self.domain[0]>=self.domain[1]:raise ValueError('Positive spline support interval required')
        self.curve=BSpline(self.knots,np.eye(self.count),self.degree,extrapolate=False)
        self.x=PchipInterpolator(values[:,0],values[:,1]);self.z=PchipInterpolator(values[:,0],values[:,2])
        self.station_x_bounds=(float(values[:,1].min()),float(values[:,1].max()))

    def design(self,points,*,derivatives=False):
        p=np.asarray(points,float)
        if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite N by 3 source points required')
        y=p[:,1];cy=np.clip(y,self.x.x[0],self.x.x[-1]);cx=self.x(cy);cz=self.z(cy)
        sx=self.x.derivative()(cy);sz=self.z.derivative()(cy)
        sx[(y<self.x.x[0])|(y>self.x.x[-1])]=0;sz[(y<self.z.x[0])|(y>self.z.x[-1])]=0
        offset=p[:,[0,2]]-np.c_[cx,cz];absolute=abs(offset);width=self.support-self.core
        t=np.clip((absolute-self.core)/width,0,1);w=(1-t)**4*(4*t+1)
        dw=np.where((absolute>self.core)&(absolute<self.support),-20*t*(1-t)**3/width*np.sign(offset),0.)
        weight=w.prod(1);gx=dw[:,0]*w[:,1];gz=w[:,0]*dw[:,1];gy=-sx*gx-sz*gz
        yy=np.clip(y,*self.domain);basis=self.curve(yy);basis_d=self.curve.derivative()(yy)
        active=(y>self.domain[0])&(y<self.domain[1])&np.all(absolute<self.support,axis=1)
        basis[~active]=0;basis_d[~active]=0;matrix=weight[:,None]*basis
        if not derivatives:return matrix
        derivative=basis[:,:,None]*np.c_[gx,gy,gz][:,None,:];derivative[:,:,1]+=weight[:,None]*basis_d
        return matrix,derivative

    def apply(self,points,controls):
        p=np.asarray(points,float);c=np.asarray(controls,float)
        if c.shape!=(self.count,2)or not np.isfinite(c).all():raise ValueError('One finite X/Z translation per spline control required')
        if np.any(c[:3])or np.any(c[-3:]):raise ValueError('First/last three controls must retain the fixed C2 endpoint jets')
        lower=1-2.109375*np.sum(np.max(abs(c),axis=0)/(self.support-self.core))
        if lower<=0:raise ValueError('Source course cannot certify a positive Jacobian')
        matrix,derivative=self.design(p,derivatives=True);delta=matrix@c;q=p.copy();q[:,[0,2]]+=delta
        jac=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();jac[:,0,:]+=np.einsum('nka,k->na',derivative,c[:,0]);jac[:,2,:]+=np.einsum('nka,k->na',derivative,c[:,1])
        return q,jac,dict(kind='compact-transverse-cubic-spline',analytic_determinant_lower_bound=float(lower),Y_fixed=True,
            fixed_outside_support=True,fixed_C2_endpoint_jets=True,core_transverse_sections_rigid=True,discrete_mesh_checks_pending=True)
