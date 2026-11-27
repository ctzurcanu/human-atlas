"""A globally invertible smooth cage with monotone longitudinal coordinate.

Cross sections undergo a proper rotation and positive principal stretches.
They can follow different organ contours without reversing or collapsing the
continuous coordinate field. Meshed samples still require intersection checks.
Tube lumens should use the separate rigid-section transport, not cage shear.
"""
from dataclasses import dataclass
from collections.abc import Mapping
import numpy as np
from scipy.interpolate import PchipInterpolator


def polygon_moments(polygons):
    """Area, centroid and covariance of disjoint simple filled polygons."""
    pieces=[]
    for polygon in polygons:
        q=np.asarray(polygon,float);origin=q.mean(0);q=q-origin;p=np.roll(q,-1,axis=0);cross=q[:,0]*p[:,1]-p[:,0]*q[:,1];area=cross.sum()/2
        if area<0:q=q[::-1];p=np.roll(q,-1,axis=0);cross=q[:,0]*p[:,1]-p[:,0]*q[:,1];area=cross.sum()/2
        if area<=0:raise ValueError('Positive finite section area required')
        center=((q+p)*cross[:,None]).sum(0)/(6*area)
        xx=np.sum((q[:,0]**2+q[:,0]*p[:,0]+p[:,0]**2)*cross)/(12*area)
        yy=np.sum((q[:,1]**2+q[:,1]*p[:,1]+p[:,1]**2)*cross)/(12*area)
        xy=np.sum((2*q[:,0]*q[:,1]+q[:,0]*p[:,1]+p[:,0]*q[:,1]+2*p[:,0]*p[:,1])*cross)/(24*area)
        cov=np.array([[xx,xy],[xy,yy]])-np.outer(center,center);pieces.append((area,center+origin,cov))
    total=sum(a for a,c,v in pieces);center=sum(a*c for a,c,v in pieces)/total;cov=sum(a*(v+np.outer(c-center,c-center))for a,c,v in pieces)/total
    if np.linalg.eigvalsh(cov).min()<=0:raise ValueError('A nondegenerate finite cross-section is required')
    return total,center,cov


def principal_chart(covariances):
    cov=np.asarray(covariances,float);eig,vec=np.linalg.eigh(cov);eig=eig[:,::-1];vec=vec[:,:,::-1]
    angle=np.unwrap(2*np.arctan2(vec[:,1,0],vec[:,0,0]))/2
    return eig,angle


@dataclass
class LongitudinalOrganCage:
    source_bounds:np.ndarray
    target_bounds:np.ndarray
    fractions:np.ndarray
    source_centers:np.ndarray
    target_centers:np.ndarray
    source_covariances:np.ndarray
    target_covariances:np.ndarray
    minimum_stretch:float=0.
    protected_basal_lift:dict|None=None

    def __post_init__(self):
        for key in self.__dataclass_fields__:
            if key!='protected_basal_lift':setattr(self,key,np.asarray(getattr(self,key),float))
        if self.source_bounds.shape!=(2,)or self.target_bounds.shape!=(2,)or np.any(np.diff(self.source_bounds)<=0)or np.any(np.diff(self.target_bounds)<=0):raise ValueError('Strictly increasing source/target longitudinal bounds required')
        self._basal_lift=None
        if self.protected_basal_lift is not None:self._declare_protected_basal_lift()
        sm,sa=principal_chart(self.source_covariances);tm,ta=principal_chart(self.target_covariances)
        # PCA axes are unoriented; choose the continuous proper target chart
        # closest to the donor instead of introducing arbitrary half turns.
        ta=ta+np.pi*np.rint((sa[0]-ta[0])/np.pi)
        self.scale=(self.target_bounds[1]-self.target_bounds[0])/(self.source_bounds[1]-self.source_bounds[0])
        self._source=PchipInterpolator(self.fractions,self.source_centers)
        self._target=PchipInterpolator(self.fractions,self.target_centers)
        self._angles=PchipInterpolator(self.fractions,np.c_[sa,ta])
        stretches=np.sqrt(tm/sm)
        if self.minimum_stretch<0:raise ValueError('A nonnegative lower stretch bound is required')
        stretches=np.maximum(stretches,self.minimum_stretch)
        self._stretch=PchipInterpolator(self.fractions,np.log(stretches))

    def _require_identity_base(self):
        for first,second in [('source_bounds','target_bounds'),('source_centers','target_centers'),('source_covariances','target_covariances')]:
            a,b=getattr(self,first),getattr(self,second)
            if not np.isfinite(a).all()or not np.array_equal(a,b):
                raise ValueError('Protected basal lift requires identical finite source/target bounds, centers and covariances')
        if self.minimum_stretch.shape!=()or not np.isfinite(self.minimum_stretch)or not 0<=self.minimum_stretch<=1:
            raise ValueError('Protected basal lift requires an identity base; minimum_stretch must be between zero and one')

    def _declare_protected_basal_lift(self):
        """Explicit physical Y lift in an already rigidly rotated Z/X/Y chart.

        Source/target sections remain identical. The added field is triangular:
        Z and X stay fixed and Y is strictly increasing, with limits at both
        infinities. Thus the analytic determinant lower bound also establishes
        global invertibility, not merely a positive sampled determinant.
        """
        spec=self.protected_basal_lift
        required={'amplitude_m','chart_axes','z_fall_bounds_m','y_fall_bounds_m'}
        if not isinstance(spec,Mapping)or set(spec)!=required:
            raise ValueError('Protected basal lift requires amplitude_m, chart_axes, z_fall_bounds_m and y_fall_bounds_m only')
        if not isinstance(spec['chart_axes'],(list,tuple))or tuple(spec['chart_axes'])!=('Z','X','Y'):
            raise ValueError('Protected basal lift supports only the declared rigid Z/X/Y chart')
        amplitude=np.asarray(spec['amplitude_m'],float)
        if amplitude.shape!=()or not np.isfinite(amplitude)or amplitude<0:
            raise ValueError('Protected basal lift amplitude must be a finite nonnegative scalar in metres')
        ranges=[]
        for key in ['z_fall_bounds_m','y_fall_bounds_m']:
            value=np.asarray(spec[key],float)
            if value.shape!=(2,)or not np.isfinite(value).all()or value[1]<=value[0]:
                raise ValueError('Protected basal lift fall bounds must be two finite strictly increasing metre coordinates')
            ranges.append(tuple(map(float,value)))
        self._require_identity_base()
        lower=1-1.5*float(amplitude)/(ranges[1][1]-ranges[1][0])
        if lower<=0:
            raise ValueError('Protected basal lift has a nonpositive global Jacobian determinant lower bound')
        self._basal_lift=(float(amplitude),ranges[0],ranges[1])
        self.protected_basal_lift_jacobian_lower_bound=lower

    @staticmethod
    def _c1_fall(values,bounds):
        u=np.clip((values-bounds[0])/(bounds[1]-bounds[0]),0.,1.)
        fall=1-u*u*(3-2*u)
        derivative=-6*u*(1-u)/(bounds[1]-bounds[0])
        return fall,derivative

    def _basal_terms(self,points):
        self._require_identity_base()
        if not np.isfinite(points).all():raise ValueError('Finite points required for the protected basal lift')
        amplitude,z_bounds,y_bounds=self._basal_lift
        z,dz=self._c1_fall(points[:,0],z_bounds)
        y,dy=self._c1_fall(points[:,2],y_bounds)
        return amplitude*z*y,amplitude*dz*y,amplitude*z*dy

    def sections(self,x):
        u=np.clip((np.asarray(x,float)-self.source_bounds[0])/np.diff(self.source_bounds)[0],self.fractions[0],self.fractions[-1]);angles=self._angles(u);stretches=np.exp(self._stretch(u))
        def rot(angle):
            c,s=np.cos(angle),np.sin(angle);return np.stack([c,-s,s,c],axis=-1).reshape(-1,2,2)
        before,after=rot(angles[:,0]),rot(angles[:,1]);matrix=(after*stretches[:,None,:])@before.transpose(0,2,1)
        return self._source(u),self._target(u),matrix

    def map(self,points):
        points=np.asarray(points,float);shape=points.shape;points=points.reshape(-1,3)
        if self._basal_lift is not None:
            lift,_,_=self._basal_terms(points);out=points.copy();out[:,2]+=lift
            return out.reshape(shape)
        source,target,matrix=self.sections(points[:,0]);out=np.empty_like(points);out[:,0]=self.target_bounds[0]+self.scale*(points[:,0]-self.source_bounds[0]);out[:,1:]=target+np.einsum('pij,pj->pi',matrix,points[:,1:]-source);return out.reshape(shape)

    def jacobians(self,points,step=1e-6):
        points=np.asarray(points,float).reshape(-1,3)
        if self._basal_lift is not None:
            _,dz,dy=self._basal_terms(points);j=np.broadcast_to(np.eye(3),(len(points),3,3)).copy();j[:,2,0]=dz;j[:,2,2]+=dy
        else:
            _,_,matrix=self.sections(points[:,0]);delta=np.array([step,0.,0.]);first=(self.map(points+delta)-self.map(points-delta))/(2*step);j=np.zeros((len(points),3,3));j[:,:,0]=first;j[:,1:,1:]=matrix
        if np.any(np.linalg.det(j)<=0):raise ValueError('Cage reversal or collapse')
        return j

    def tube_field(self,field,section_scale):
        """Move centres with the cage, rigidly rotate every paired section.

        The section scale is independent of organ contour stretches. All
        inner/outer polygon areas therefore change only by section_scale².
        """
        from .tube_sections import path_tangents,rotation_minimizing_frames
        field=np.asarray(field,float);old=field[0].mean(1);new=self.map(old);ot,nt=path_tangents(old),path_tangents(new);before=rotation_minimizing_frames(old,ot,initial_reference=field[1,0,0]-old[0]);initial=self.jacobians(old[:1])[0]@before[0,:,0];after=rotation_minimizing_frames(new,nt,initial_reference=initial);rotations=after@before.transpose(0,2,1)
        if not np.allclose(np.linalg.det(rotations),1,atol=1e-12):raise ValueError('Proper section rotations required')
        return new[None,:,None]+section_scale*np.einsum('sij,ksvj->ksvi',rotations,field-old[None,:,None])
