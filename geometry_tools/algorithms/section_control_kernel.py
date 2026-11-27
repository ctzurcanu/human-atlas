"""Planar untangling in actual shared source controls, without point welding.

points + weights @ control_delta defines the only allowed motion. For one
2D control, each oriented triangle area is affine even when several of its
vertices depend on that control: the quadratic cross(delta,delta) vanishes.
This is a local constraint repair, not a global injectivity/anatomy proof.
"""
import time
import numpy as np


def repair_section_controls(points, faces, weights, *, minimum_double_area,
                            maximum_displacement_m, maximum_steps=32, wall_seconds=1.):
    started=time.monotonic()
    source=np.asarray(points,float)
    f=np.asarray(faces,int)
    A=np.asarray(weights,float)
    floor=np.broadcast_to(np.asarray(minimum_double_area,float),(len(f),))
    if source.ndim!=2 or source.shape[1]!=2 or A.ndim!=2 or len(A)!=len(source):
        raise ValueError('Planar points and explicit source-control weights required')
    if (not np.isfinite(np.r_[source.ravel(),A.ravel(),floor.ravel()]).all()
            or np.any(floor<=0) or maximum_displacement_m<=0
            or maximum_steps<1 or not 0<wall_seconds<=10):
        raise ValueError('Finite geometry and positive finite repair budgets required')
    delta=np.zeros((A.shape[1],2))
    q=source.copy()
    def areas(x):
        t=x[f];a=t[:,1]-t[:,0];b=t[:,2]-t[:,0]
        return a[:,0]*b[:,1]-a[:,1]*b[:,0]
    def clip(poly,n,b):
        out=[];v=poly@n-b
        for i,j in zip(range(len(poly)),np.roll(np.arange(len(poly)),-1)):
            if v[i]>=0:out.append(poly[i])
            if (v[i]>=0)!=(v[j]>=0):
                w=v[i]/(v[i]-v[j]);out.append(poly[i]*(1-w)+poly[j]*w)
        return np.asarray(out).reshape(-1,2)
    current=areas(q)
    initial=int(np.count_nonzero(current<floor))
    moves=[];reason=None
    touched=[np.flatnonzero(np.any(A[f,c]!=0,axis=1)) for c in range(A.shape[1])]
    theta=np.arange(32)*2*np.pi/32
    ball=maximum_displacement_m*np.c_[np.cos(theta),np.sin(theta)]
    for _ in range(maximum_steps):
        bad=np.flatnonzero(current<floor)
        if not len(bad):break
        if time.monotonic()-started>=wall_seconds:
            reason='native wall deadline';break
        candidates=np.flatnonzero(np.any(A[np.unique(f[bad])]!=0,axis=0))
        accepted=False
        for c in candidates:
            ids=touched[c];tri=f[ids];t=q[tri]
            gradient=np.zeros((len(ids),2))
            for k in range(3):
                one,two=t[:,(k+1)%3],t[:,(k+2)%3]
                g=np.c_[one[:,1]-two[:,1],two[:,0]-one[:,0]]
                gradient+=A[tri[:,k],c,None]*g
            kernel=ball.copy()
            for n,b in zip(gradient,floor[ids]-current[ids]+gradient@delta[c]):
                kernel=clip(kernel,n,b)
                if len(kernel)<3:break
            if len(kernel)<3:continue
            choices=[]
            for one,two in zip(kernel,np.roll(kernel,-1,axis=0)):
                edge=two-one;den=float(edge@edge)
                w=0. if den==0 else np.clip(-(one@edge)/den,0,1)
                choices.append(one+w*edge)
            proposed=.999*min(choices,key=lambda p:p@p)+.001*kernel.mean(0)
            old=delta[c].copy();delta[c]=proposed
            candidate=source+A@delta
            after=areas(candidate)
            if np.any(after[ids]<floor[ids]*(1-1e-9)) or np.count_nonzero(after<floor)>=len(bad):
                delta[c]=old;continue
            q,current=candidate,after
            moves.append(dict(control=int(c),remaining_bad_faces=int(np.count_nonzero(current<floor))))
            accepted=True;break
        if not accepted:reason='no bounded shared-control kernel';break
    else:reason='step budget'
    return q,delta,dict(passed=not np.any(current<floor),initial_bad_faces=initial,
        remaining_bad_faces=np.flatnonzero(current<floor).tolist(),moves=moves,
        stop_reason=reason,maximum_control_displacement_m=float(np.linalg.norm(delta,axis=1).max(initial=0)),
        maximum_point_displacement_m=float(np.linalg.norm(q-source,axis=1).max(initial=0)),
        declared_maximum_displacement_m=maximum_displacement_m,
        exact_control_replay=bool(np.array_equal(q,source+A@delta)),
        maximum_steps=maximum_steps,native_wall_seconds=wall_seconds,
        elapsed_seconds=time.monotonic()-started)


def repair_section_control_patch(points, faces, weights, *, minimum_double_area,
                                 maximum_displacement_m, initial=None,
                                 maximum_steps=80, wall_seconds=2.):
    """Couple real controls when single-control kernels are empty.

    The caller supplies mathematical section triangles, not new tissue faces.
    Independent section/surface checks are mandatory even if this solve passes.
    """
    from scipy.optimize import minimize
    started=time.monotonic()
    original=np.asarray(points,float)
    f=np.asarray(faces,int)
    A=np.asarray(weights,float)
    floor=np.broadcast_to(np.asarray(minimum_double_area,float),(len(f),))
    if original.ndim!=2 or original.shape[1]!=2 or A.ndim!=2 or len(A)!=len(original) or not 0<A.shape[1]<=64:
        raise ValueError('A bounded patch of at most64 explicit planar controls required')
    if not np.isfinite(np.r_[original.ravel(),A.ravel(),floor.ravel()]).all() or np.any(floor<=0) or maximum_displacement_m<=0 or maximum_steps<1 or not 0<wall_seconds<=10:
        raise ValueError('Finite source and bounded repair budgets required')
    # Origin subtraction prevents large global coordinate offsets from
    # dominating derivatives. Exported coordinates are independently replayed.
    origin=original.mean(0)
    source=original-origin
    bound=maximum_displacement_m
    delta=np.zeros((A.shape[1],2))if initial is None else np.asarray(initial,float)
    if delta.shape!=(A.shape[1],2) or not np.isfinite(delta).all() or np.linalg.norm(delta,axis=1).max(initial=0)>bound:
        raise ValueError('Initial controls exceed the declared source ball')
    du=A[f[:,1]]-A[f[:,0]]
    dv=A[f[:,2]]-A[f[:,0]]
    def geometry(x):
        if time.monotonic()-started>=wall_seconds:
            raise TimeoutError('native coupled section-control deadline')
        q=source+A@x.reshape(-1,2)*bound
        t=q[f];u=t[:,1]-t[:,0];v=t[:,2]-t[:,0]
        return u,v,u[:,0]*v[:,1]-u[:,1]*v[:,0]
    _,_,before=geometry(np.zeros(delta.size))
    scales=np.maximum(abs(before),floor)
    def constraint(x):
        _,_,area=geometry(x)
        return np.r_[(area-floor)/scales,1-np.sum(x.reshape(-1,2)**2,axis=1)]
    def jacobian(x):
        u,v,_=geometry(x)
        jac=np.empty((len(f)+A.shape[1],delta.size))
        jac[:len(f),0::2]=(du*v[:,1,None]-dv*u[:,1,None])*bound/scales[:,None]
        jac[:len(f),1::2]=(dv*u[:,0,None]-du*v[:,0,None])*bound/scales[:,None]
        jac[len(f):]=0
        for i,value in enumerate(x.reshape(-1,2)):
            jac[len(f)+i,2*i:2*i+2]=-2*value
        return jac
    result=None;reason=None;passed=False
    try:
        result=minimize(lambda x:float(x@x/2),(delta/bound).ravel(),jac=lambda x:x,
            constraints=[dict(type='ineq',fun=constraint,jac=jacobian)],
            method='SLSQP',options=dict(maxiter=maximum_steps,ftol=1e-11,disp=False))
        delta=result.x.reshape(-1,2)*bound
        q=original+A@delta
        t=q[f];u=t[:,1]-t[:,0];v=t[:,2]-t[:,0]
        after=u[:,0]*v[:,1]-u[:,1]*v[:,0]
        passed=bool(result.success and np.all(after>=floor*(1-1e-8)) and np.linalg.norm(delta,axis=1).max()<=bound*(1+1e-10))
        reason=str(result.message)
    except TimeoutError as error:
        delta=np.zeros_like(delta);q=original.copy();after=before;reason=str(error)
    return q,delta,dict(passed=passed,method='coupled-real-source-controls',
        remaining_bad_faces=np.flatnonzero(after<floor*(1-1e-8)).tolist(),
        maximum_control_displacement_m=float(np.linalg.norm(delta,axis=1).max(initial=0)),
        exact_control_replay=bool(np.array_equal(q,original+A@delta)),
        maximum_steps=maximum_steps,native_wall_seconds=wall_seconds,
        iterations=int(result.nit)if result is not None else 0,
        stop_reason=reason,elapsed_seconds=time.monotonic()-started)
