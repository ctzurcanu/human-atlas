"""Bounded Gaussian-kernel velocity interpolation along landmark trajectories.

All anatomy follows one smooth time-dependent flow. This is a prescribed-path
initializer, not the minimum-energy geodesic optimizer in Joshi & Miller (2000),
doi:10.1109/83.855431. Numerical integration, discrete mesh quality, anatomy and
calibre require independent verification; theoretical smooth-flow invertibility
does not certify the emitted triangles.
"""
import time
import numpy as np
from scipy.linalg import cho_factor,cho_solve
from scipy.integrate import RK45
from scipy.spatial.distance import cdist
from scipy.interpolate import CubicSpline


def prepare_landmark_velocity(source,target,*,kernel_width_m=.018,kernel_weights=None,regularization=1e-6,time_samples=65):
    a=np.asarray(source,float);b=np.asarray(target,float)
    widths=np.atleast_1d(np.asarray(kernel_width_m,float));weights=np.ones(len(widths))if kernel_weights is None else np.asarray(kernel_weights,float)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or not 4<=len(a)<=600:
        raise ValueError('Between four and 600 explicit corresponding 3D landmarks required')
    if not np.isfinite(a).all()or not np.isfinite(b).all()or widths.ndim!=1 or not 1<=len(widths)<=3 or np.any((widths<.005)|(widths>.05)):
        raise ValueError('Finite landmarks and explicit local kernel width required')
    if weights.shape!=widths.shape or np.any(weights<=0)or not np.isfinite(weights).all():raise ValueError('Positive explicit kernel weights required')
    if not 1e-9<=regularization<=1e-2 or not 9<=time_samples<=129:
        raise ValueError('Explicit finite kernel regularization and temporal sampling required')
    times=np.linspace(0,1,time_samples);coefs=[];residual=[];minimum_separation=np.inf
    velocity=b-a
    for t in times:
        c=a+t*velocity;dist=cdist(c,c,'sqeuclidean')
        np.fill_diagonal(dist,np.inf);minimum_separation=min(minimum_separation,float(np.sqrt(dist.min())))
        np.fill_diagonal(dist,0.);K=sum(w*np.exp(-dist/(2*s**2))for s,w in zip(widths,weights))
        H=K+np.eye(len(a))*regularization
        coef=cho_solve(cho_factor(H,lower=True,check_finite=False),velocity,check_finite=False)
        coefs.append(coef);residual.append(float(np.linalg.norm(K@coef-velocity,axis=1).max()))
    if minimum_separation<1e-7:raise ValueError('Corresponding landmark trajectories collide')
    return dict(source=a,target=b,times=times,coefficients=np.array(coefs),kernel_width_m=widths,kernel_weights=weights),dict(
        method='Prescribed linear landmark trajectories with regularized Gaussian velocity interpolation',
        landmarks=len(a),timeSamples=time_samples,kernelWidths_m=widths.tolist(),kernelWeights=weights.tolist(),regularization=regularization,
        minimumSampledTrajectorySeparation_m=minimum_separation,maximumVelocityInterpolationResidual_m=max(residual),
        optimalGeodesicClaimed=False,anatomicalAcceptance=False)


def advance_landmark_flow(points,bank,*,start_time=0.,end_time=1.,wall_seconds=15.,maximum_steps=256,
                          relative_tolerance=1e-7,absolute_tolerance_m=1e-8,fixed_above=1.08,free_below=1.,
                          temporal_interpolation='linear'):
    p=np.asarray(points,float);a=np.asarray(bank['source'],float);b=np.asarray(bank['target'],float)
    times=np.asarray(bank['times'],float);coefs=np.asarray(bank['coefficients'],float);widths=np.atleast_1d(np.asarray(bank['kernel_width_m'],float))
    weights=np.asarray(bank.get('kernel_weights',np.ones(len(widths))),float)
    if not 0<=start_time<end_time<=1 or not 0<wall_seconds<=20 or not 1<=maximum_steps<=512:
        raise ValueError('Finite explicit continuation and work bounds required')
    if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()or not free_below<fixed_above:
        raise ValueError('Finite physical controls and an explicit fixed-root support required')
    if coefs.shape!=(len(times),len(a),3)or not np.all(np.diff(times)>0):raise ValueError('Valid saved temporal velocity bank required')
    if temporal_interpolation not in ['linear','cubic']:
        raise ValueError('Explicit linear or cubic temporal coefficient interpolation required')
    spline=CubicSpline(times,coefs,axis=0)if temporal_interpolation=='cubic'else None
    began=time.monotonic();evaluations=0;steps=0
    class Deadline(Exception):pass
    def rhs(t,y):
        nonlocal evaluations
        if time.monotonic()-began>=wall_seconds:raise Deadline()
        evaluations+=1;q=y.reshape(-1,3);k=int(np.clip(np.searchsorted(times,t)-1,0,len(times)-2))
        u=(t-times[k])/(times[k+1]-times[k]);alpha=spline(t)if spline is not None else(1-u)*coefs[k]+u*coefs[k+1];c=a+t*(b-a)
        out=np.empty_like(q)
        for first in range(0,len(q),2048):
            x=q[first:first+2048];dist=cdist(x,c,'sqeuclidean');K=sum(w*np.exp(-dist/(2*width**2))for width,w in zip(widths,weights))
            v=K@alpha;s=np.clip((fixed_above-x[:,1])/(fixed_above-free_below),0,1)
            out[first:first+len(x)]=v*(s*s*(3-2*s))[:,None]
        return out.ravel()
    solver=RK45(rhs,start_time,p.ravel(),end_time,rtol=relative_tolerance,atol=absolute_tolerance_m,max_step=1/16)
    timed_out=False
    try:
        while solver.status=='running'and steps<maximum_steps:
            if time.monotonic()-began>=wall_seconds:raise Deadline()
            solver.step();steps+=1
    except Deadline:timed_out=True
    complete=solver.status=='finished'and not timed_out
    q=solver.y.reshape(-1,3).copy();fixed=p[:,1]>=fixed_above
    if not np.array_equal(q[fixed],p[fixed]):raise ValueError('Fixed proximal physical controls moved')
    return q,dict(complete=complete,status='Integrated'if complete else'NativeDeadline'if timed_out else str(solver.status),
        startTime=float(start_time),endTime=float(solver.t),requestedEndTime=float(end_time),wholeFlowComplete=bool(complete and end_time==1.),acceptedSteps=steps,rhsEvaluations=evaluations,
        relativeTolerance=relative_tolerance,absoluteTolerance_m=absolute_tolerance_m,
        temporalCoefficientInterpolation=temporal_interpolation,
        fixedAbove_m=fixed_above,freeBelow_m=free_below,fixedOriginalControls=int(fixed.sum()),
        nativeWallSeconds=wall_seconds,elapsedSeconds=time.monotonic()-began,
        timedOut=timed_out,maximumStepsExhausted=steps>=maximum_steps and not complete,
        anatomicalAcceptance=False,discreteMeshAndCalibreAuditRequired=True)
