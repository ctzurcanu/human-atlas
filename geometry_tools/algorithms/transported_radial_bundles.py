"""Finite positive transverse bundle maps with exact rigid flat sections.

This module is independent of the application dispatcher. Every ambient point
receives the same maps; object IDs, masks and anatomical installation are absent.
"""
from functools import lru_cache
import numpy as np

MAX_WEIGHT_SLOPE=135./64.
_FIELDS=('source_x_m','source_z_m','target_x_m','target_z_m','angle_radians','core_m','width_m')

@lru_cache(maxsize=64)
def _profiles(stations,support):
    from scipy.interpolate import PchipInterpolator
    a=np.asarray(stations,float);ys=np.asarray(support,float)
    # Constant endcaps make the profiles C1 at the seven-guide endpoints.
    values=np.column_stack([a[:,1:3],a[:,3:5]-a[:,1:3],a[:,5:]])
    y=np.r_[ys[0],a[:,0],ys[1]];v=np.r_[values[:1],values,values[-1:]]
    curves=PchipInterpolator(y,v,axis=0,extrapolate=False)
    return curves,curves.derivative()

def _controls(spec):
    try:
        stations=tuple(tuple(float(s[k])for k in ('y_m',*_FIELDS))for s in spec['stations'])
        yc=np.asarray(spec['y_core_interval_m'],float);ys=np.asarray(spec['y_support_interval_m'],float)
        steps=spec['steps']
    except (KeyError,TypeError,ValueError)as exc:raise ValueError('Explicit finite bundle controls required')from exc
    a=np.asarray(stations,float)
    if (a.ndim!=2 or a.shape[1]!=8 or len(a)<2 or not np.isfinite(a).all()or np.any(np.diff(a[:,0])<=0)
            or yc.shape!=(2,)or ys.shape!=(2,)or not np.isfinite(np.r_[yc,ys]).all()
            or not ys[0]<yc[0]<=a[0,0]<a[-1,0]<=yc[1]<ys[1]
            or np.any(a[:,6]<0)or np.any(a[:,7]<=0)
            or isinstance(steps,bool)or not isinstance(steps,(int,np.integer))or not 1<=steps<=512):
        raise ValueError('Ordered guides, positive widths and finite 1..512 steps required')
    delta=a[:,3:5]-a[:,1:3]
    # PCHIP preserves the range of each delta component; use its full rectangle.
    maximum_delta=float(np.linalg.norm(np.max(np.abs(delta),axis=0)))
    width_min=float(a[:,7].min());lower_step=1-MAX_WEIGHT_SLOPE*maximum_delta/(steps*width_min)
    if lower_step<=0:raise ValueError('Each finite radial translation must certify positive determinant')
    return stations,yc,ys,int(steps),dict(analytic_each_translation_determinant_lower_bound=float(lower_step),
        analytic_determinant_lower_bound=float(lower_step**steps),twist_determinant=1.,
        maximum_translation_norm_bound_m=maximum_delta,minimum_transition_width_m=width_min,
        finite_composition_not_ODE=True,rigid_XZ_sections_inside_flat_cores=True,Y_fixed=True,
        interpolation='Componentwise shape-preserving PCHIP with constant support-end caps',
        full_Y_derivatives_included=True,discrete_mesh_and_packing_checks_pending=True)

def _window(y,yc,ys):
    low=y<yc[0];gap=np.where(low,yc[0]-ys[0],ys[1]-yc[1]);distance=np.where(low,yc[0]-y,y-yc[1])
    t=np.clip(distance/gap,0,1);w=(1-t)**4*(4*t+1)
    dw=np.where((distance>0)&(distance<gap),-20*t*(1-t)**3/gap*np.where(low,-1.,1.),0.)
    return w,dw

def _radial(q,centre_y,core,core_y,width,width_y):
    radius=np.linalg.norm(q,axis=1);unit=np.divide(q,radius[:,None],out=np.zeros_like(q),where=radius[:,None]>0)
    t=np.clip((radius-core)/width,0,1);w=(1-t)**4*(4*t+1)
    active=(radius>core)&(radius<core+width);derivative=np.where(active,-20*t*(1-t)**3,0.)
    grad=derivative[:,None]/width[:,None]*unit
    radius_y=-np.einsum('ij,ij->i',unit,centre_y)
    t_y=(radius_y-core_y)/width-(radius-core)*width_y/(width*width)
    w_y=derivative*t_y
    return w,grad,w_y

def map_bundle(points,spec):
    """Apply small moving-centre translations then an exact radial XZ twist.

    A Y section in the complete moving flat core receives exactly one proper
    rigid transform. Each translation uses a rank-one determinant bound; the
    radial twist has determinant one. Y-chart slopes enter the full Jacobian
    but cannot change its determinant because Y itself remains literal.
    """
    p=np.asarray(points,float)
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite N by 3 metre coordinates required')
    stations,yc,ys,steps,proof=_controls(spec)
    out=p.copy();jac=np.broadcast_to(np.eye(3),(len(p),3,3)).copy()
    ids=np.flatnonzero((p[:,1]>ys[0])&(p[:,1]<ys[1]))
    if not len(ids):return out,jac,proof
    curves,derivatives=_profiles(stations,tuple(ys));y=p[ids,1];a=curves(y);da=derivatives(y)
    centre,centre_y=a[:,:2],da[:,:2];delta,delta_y=a[:,2:4],da[:,2:4]
    angle,angle_y=a[:,4],da[:,4];core,core_y=a[:,5],da[:,5];width,width_y=a[:,6],da[:,6]
    wy,dwy=_window(y,yc,ys);travel=wy[:,None]*delta;travel_y=dwy[:,None]*delta+wy[:,None]*delta_y
    # Conservative complete swept support; this is spatial filtering, never ownership filtering.
    active=np.linalg.norm(p[ids][:,[0,2]]-centre,axis=1)<core+width+np.linalg.norm(travel,axis=1)
    ids=ids[active];aids=np.flatnonzero(active)
    if not len(ids):return out,jac,proof
    centre,centre_y,travel,travel_y,angle,angle_y,core,core_y,width,width_y,wy,dwy=[v[aids]for v in [centre,centre_y,travel,travel_y,angle,angle_y,core,core_y,width,width_y,wy,dwy]]
    q=p[ids].copy();chain=jac[ids].copy();step=travel/steps;step_y=travel_y/steps
    for k in range(steps):
        fraction=k/steps;c=centre+fraction*travel;cy=centre_y+fraction*travel_y
        w,g,gy=_radial(q[:,[0,2]]-c,cy,core,core_y,width,width_y)
        local=np.broadcast_to(np.eye(3),(len(q),3,3)).copy()
        for row,j in enumerate((0,2)):
            local[:,j,0]+=step[:,row]*g[:,0];local[:,j,2]+=step[:,row]*g[:,1]
            local[:,j,1]+=step[:,row]*gy+w*step_y[:,row]
        q[:,[0,2]]+=w[:,None]*step;chain=local@chain
    c=centre+travel;cy=centre_y+travel_y;offset=q[:,[0,2]]-c
    w,g,gy=_radial(offset,cy,core,core_y,width,width_y)
    theta=wy*angle;theta_y=dwy*angle+wy*angle_y;phi=theta*w
    co,si=np.cos(phi),np.sin(phi);rot=np.empty((len(q),2,2));rot[:,0,0]=rot[:,1,1]=co;rot[:,0,1]=-si;rot[:,1,0]=si
    rotated=np.einsum('nij,nj->ni',rot,offset);krot=np.column_stack([-rotated[:,1],rotated[:,0]])
    gp=theta[:,None]*g;py=theta_y*w+theta*gy
    block=rot+krot[:,:,None]*gp[:,None,:]
    column_y=cy-np.einsum('nij,nj->ni',rot,cy)+krot*py[:,None]
    local=np.broadcast_to(np.eye(3),(len(q),3,3)).copy()
    local[:,0,0]=block[:,0,0];local[:,0,2]=block[:,0,1];local[:,2,0]=block[:,1,0];local[:,2,2]=block[:,1,1]
    local[:,0,1]=column_y[:,0];local[:,2,1]=column_y[:,1]
    moved=c+rotated;nonzero=phi!=0
    q[nonzero,0]=moved[nonzero,0];q[nonzero,2]=moved[nonzero,1]
    chain=local@chain;out[ids]=q;jac[ids]=chain
    return out,jac,proof

def compose_bundles(points,specs):
    """Compose the same ambient bundle fields for every anatomical member."""
    if not isinstance(specs,(tuple,list))or not specs:raise ValueError('Nonempty common field sequence required')
    q=np.asarray(points,float);jac=np.broadcast_to(np.eye(3),(len(q),3,3)).copy();proofs=[];lower=1.
    for spec in specs:
        q,j,p=map_bundle(q,spec);jac=j@jac;lower*=p['analytic_determinant_lower_bound'];proofs.append(p)
    return q,jac,dict(analytic_determinant_lower_bound=float(lower),fields=proofs,
        Y_fixed=True,finite_composition_not_ODE=True,discrete_mesh_and_anatomical_checks_pending=True)
