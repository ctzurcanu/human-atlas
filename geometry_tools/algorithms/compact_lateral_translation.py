"""An explicitly bounded C2 bilateral lateral translation field.

Y and Z are unchanged. At fixed Y/Z, X remains strictly increasing when the
analytic determinant lower bound is positive, giving a global injective map.
This geometric property does not establish anatomy, skin clearance or calibre.
"""
from dataclasses import dataclass
import numpy as np

def plateau(values,bounds):
    a,b,c,d=map(float,bounds)
    if not a<b<=c<d:raise ValueError('Ordered outer/inner support limits required')
    x=np.asarray(values,float);w=np.zeros_like(x);derivative=np.zeros_like(x)
    for mask,t,scale in [( (x>a)&(x<b),(x-a)/(b-a),1/(b-a)),((x>c)&(x<d),(d-x)/(d-c),-1/(d-c))]:
        q=t[mask];w[mask]=q*q*q*(10+q*(-15+6*q));derivative[mask]=30*q*q*(q-1)*(q-1)*scale
    w[(x>=b)&(x<=c)]=1
    return w,derivative

@dataclass(frozen=True)
class CompactBilateralLateralTranslation:
    distance_m:float
    absolute_x_support:tuple
    y_support:tuple
    z_support:tuple

    @property
    def determinant_lower_bound(self):
        a,b,c,d=self.absolute_x_support
        return 1-abs(self.distance_m)*1.875/min(b-a,d-c)

    def apply(self,points):
        p=np.asarray(points,float)
        if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite Nx3 world-metre positions required')
        if self.absolute_x_support[0]<=0 or self.determinant_lower_bound<=0:raise ValueError('Support must avoid midline and have a positive analytic Jacobian bound')
        wx,dx=plateau(np.abs(p[:,0]),self.absolute_x_support);wy,dy=plateau(p[:,1],self.y_support);wz,dz=plateau(p[:,2],self.z_support);side=np.sign(p[:,0]);q=p.copy();q[:,0]+=side*self.distance_m*wx*wy*wz
        J=np.broadcast_to(np.eye(3),(len(p),3,3)).copy();J[:,0,0]+=self.distance_m*dx*wy*wz;J[:,0,1]=side*self.distance_m*wx*dy*wz;J[:,0,2]=side*self.distance_m*wx*wy*dz
        if (J[:,0,0]<self.determinant_lower_bound-1e-12).any():raise ArithmeticError('Analytic determinant bound violated')
        return q,J

    def normals(self,points,normals):
        _,J=self.apply(points);n=np.asarray(normals,float)
        if n.shape!=np.asarray(points).shape:raise ValueError('Normals must correspond to input positions')
        out=np.linalg.solve(J.swapaxes(1,2),n[...,None])[...,0];length=np.linalg.norm(out,axis=1)
        if (length==0).any():raise ValueError('Nonzero source normals required')
        return out/length[:,None]
