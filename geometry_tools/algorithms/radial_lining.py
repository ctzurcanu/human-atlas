"""Source-witnessed radial apposition of an existing finite tube wall.

R_new² = R_old² - fraction * window * taper * R_lumen².
Without the outer taper this preserves each radial tissue sector's area.
A monotone outward taper keeps the continuous radial Jacobian positive.
These facts do not establish a valid emitted mesh or a physiological state.
"""
import numpy as np


class RadialLining:
    def __init__(self,mesh):
        self.mesh=mesh;self.triangles=mesh.triangles
        self.e1=self.triangles[:,1]-self.triangles[:,0]
        self.e2=self.triangles[:,2]-self.triangles[:,0]
        self.normals=np.cross(self.e1,self.e2)
        radial=np.einsum('ti,tki->tk',self.normals[:,[0,2]],self.triangles[:,:,[0,2]])
        if not len(radial)or np.any(radial>=0):
            raise ValueError('The supplied wetted lining is not strictly radial about this explicit axis')

    def radii(self,axial,angles,batch_size=128):
        axial=np.asarray(axial,float);angles=np.asarray(angles,float)
        if axial.ndim!=1 or angles.shape!=axial.shape or not np.isfinite(np.r_[axial,angles]).all():
            raise ValueError('Finite corresponding axial coordinates and angles required')
        if type(batch_size)is not int or not 0<batch_size<=256:raise ValueError('A bounded ray batch is required')
        output=[];face_ids=[];weights=[];maximum_spread=0.
        edge_scale=np.linalg.norm(self.e1,axis=1)*np.linalg.norm(self.e2,axis=1)
        for start in range(0,len(axial),batch_size):
            a=axial[start:start+batch_size];theta=angles[start:start+batch_size]
            direction=np.c_[np.cos(theta),np.zeros(len(theta)),np.sin(theta)]
            origins=np.c_[np.zeros(len(a)),a,np.zeros(len(a))]
            pvec=np.cross(direction[:,None,:],self.e2[None]);det=np.einsum('ti,bti->bt',self.e1,pvec)
            delta=origins[:,None,:]-self.triangles[None,:,0,:]
            qvec=np.cross(delta,self.e1[None])
            with np.errstate(divide='ignore',invalid='ignore'):
                u=np.einsum('bti,bti->bt',delta,pvec)/det
                v=np.einsum('bi,bti->bt',direction,qvec)/det
                distance=np.einsum('ti,bti->bt',self.e2,qvec)/det
                error=64*np.finfo(float).eps*np.maximum(1.,edge_scale[None]/abs(det))
            valid=np.isfinite(distance)&(distance>0)&(u>=-error)&(v>=-error)&(u+v<=1+error)&(det!=0)
            candidates=np.where(valid,distance,np.inf);selected=np.argmin(candidates,axis=1)
            nearest=candidates[np.arange(len(a)),selected]
            if not np.isfinite(nearest).all():raise ValueError('A radial query has no literal source-face intersection')
            farthest=np.max(np.where(valid,distance,-np.inf),axis=1)
            spread=farthest-nearest;bound=256*np.finfo(float).eps*np.maximum(1.,abs(nearest))
            if np.any(spread>bound):raise ValueError('A radial query intersects distinct source lining depths')
            maximum_spread=max(maximum_spread,float(spread.max(initial=0.)))
            uv=np.c_[u[np.arange(len(a)),selected],v[np.arange(len(a)),selected]]
            bary=np.c_[1-uv.sum(1),uv];bary=np.maximum(bary,0.);bary/=bary.sum(1)[:,None]
            actual=np.einsum('bi,bij->bj',bary,self.triangles[selected])
            radius=np.einsum('bi,bi->b',actual-origins,direction)
            residual=np.linalg.norm(actual-origins-radius[:,None]*direction,axis=1)
            if np.any(residual>bound):raise ValueError('A source-face witness drifted off its radial ray')
            output.extend(radius);face_ids.extend(selected);weights.extend(bary)
        return np.asarray(output),dict(face_ids=np.asarray(face_ids,int),barycentric=np.asarray(weights,float).reshape(-1,3),
                                      maximum_shared_hit_spread_m=maximum_spread,source_queries=len(axial),
                                      float64_source_witnesses=True)


def radial_apposition(local_points,source_radii,window,fraction,taper_distance,*,core_fraction=.25):
    p=np.asarray(local_points,float);r0=np.asarray(source_radii,float);window=np.asarray(window,float)
    if p.ndim!=2 or p.shape[1]!=3 or r0.shape!=(len(p),)or window.shape!=(len(p),):raise ValueError('Corresponding local points, source radii and windows required')
    if not 0<=fraction<=1 or not np.isfinite(taper_distance)or taper_distance<=0:raise ValueError('Explicit bounded closure and positive outer support required')
    if not np.isfinite(np.r_[p.ravel(),r0,window]).all()or np.any(r0<=0)or np.any((window<0)|(window>1)):
        raise ValueError('Finite positive source radii and unit windows required')
    radius=np.linalg.norm(p[:,[0,2]],axis=1);gap=radius-r0;bound=256*np.finfo(float).eps*np.maximum(1.,radius)
    if np.any(gap < -bound):raise ValueError('An occupied tissue control lies inside the source fluid envelope')
    if not np.isfinite(core_fraction)or not 0<=core_fraction<1:raise ValueError('A finite inner core inside the outer support is required')
    t=np.clip((np.maximum(gap,0.)/taper_distance-core_fraction)/(1-core_fraction),0.,1.);taper=(1-t)**3*(1+3*t+6*t*t)
    squared=radius**2-fraction*window*taper*r0**2
    if np.any(squared < -bound*np.maximum(radius,r0)):raise ValueError('The closure would reverse a radial tissue sector')
    target=np.sqrt(np.maximum(squared,0.));out=p.copy();scale=target/radius
    out[:,0]*=scale;out[:,2]*=scale
    return out,dict(radius_before=radius,radius_after=target,taper=taper,radial_gap=gap,
                    continuous_sector_area_preserved_where_taper_is_one=True,
                    emitted_geometry_verified=False,anatomical_acceptance=False)


def radial_jacobians(local_points,source_radii,source_gradients,window,window_derivative,fraction,taper_distance,*,core_fraction=.25):
    p=np.asarray(local_points,float);r0=np.asarray(source_radii,float);gradient=np.asarray(source_gradients,float)
    window=np.asarray(window,float);window_derivative=np.asarray(window_derivative,float)
    if gradient.shape!=(len(p),2)or window_derivative.shape!=(len(p),)or not np.isfinite(np.r_[gradient.ravel(),window_derivative]).all():
        raise ValueError('Explicit finite source and window derivatives required')
    q,field=radial_apposition(p,r0,window,fraction,taper_distance,core_fraction=core_fraction)
    radius=field['radius_before'];target=field['radius_after'];taper=field['taper']
    if np.any(target<=0):raise ValueError('Full coaptation needs remeshed contact normals, not a singular Jacobian')
    t=np.clip((np.maximum(field['radial_gap'],0.)/taper_distance-core_fraction)/(1-core_fraction),0,1)
    taper_prime=-30*t*t*(1-t)**2/(taper_distance*(1-core_fraction))
    squared_r=2*radius-fraction*window*taper_prime*r0*r0
    squared_a=-fraction*(window_derivative*taper*r0*r0+window*gradient[:,0]*(2*taper*r0-taper_prime*r0*r0))
    squared_theta=-fraction*window*gradient[:,1]*(2*taper*r0-taper_prime*r0*r0)
    dr=squared_r/(2*target);da=squared_a/(2*target);dt=squared_theta/(2*target)
    c=p[:,0]/radius;s=p[:,2]/radius;ratio=target/radius
    j=np.zeros((len(p),3,3));j[:,1,1]=1.
    j[:,0,0]=c*c*dr-c*s*dt/radius+ratio*s*s
    j[:,0,2]=c*s*dr+c*c*dt/radius-ratio*s*c
    j[:,0,1]=c*da
    j[:,2,0]=s*c*dr-s*s*dt/radius-ratio*c*s
    j[:,2,2]=s*s*dr+s*c*dt/radius+ratio*c*c
    j[:,2,1]=s*da
    determinant=squared_r/(2*radius)
    if np.any(determinant<=0)or not np.isfinite(j).all():raise ValueError('Nonpositive or nonfinite continuous radial Jacobian')
    return j,dict(continuous_determinant=determinant,emitted_geometry_verified=False)
