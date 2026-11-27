"""One deterministic, metric-bounded displacement field on an explicit patch.

The caller owns the anatomical patch, fixed controls and support direction.
This does not identify a junction, weld surfaces, or accept a candidate.
"""
import time
import numpy as np
from scipy.optimize import minimize


def separate_patch(vertices, faces, first, second, free, direction, *,
                   gap_m, maximum_motion_m, maximum_gradient, wall_seconds):
    start = time.monotonic()
    p = np.asarray(vertices, float); f = np.asarray(faces, int)
    free = np.asarray(sorted(set(map(int, free))), int)
    n = np.asarray(direction, float)
    if not len(free) or not np.isfinite(p).all() or abs(np.linalg.norm(n)-1) > 1e-10:
        raise ValueError('Explicit free controls, finite positions and unit direction required')
    if min(gap_m, maximum_motion_m, maximum_gradient, wall_seconds) <= 0:
        raise ValueError('Positive geometric bounds and deadline required')
    # Micrometres condition this small local solve without altering world data.
    scale = 1e6
    x = (p-p[np.r_[first,second]].mean(0))*scale
    index = np.full(len(p), -1, int); index[free] = np.arange(len(free))
    touched = np.flatnonzero(np.any(np.isin(f, free), axis=1))
    tri = x[f[touched]]
    E = np.swapaxes(tri[:,1:]-tri[:,:1], 1, 2)
    if np.any(np.linalg.norm(np.cross(E[:,:,0], E[:,:,1]), axis=1) <= 1e-15):
        raise ValueError('Nondegenerate source facets required')
    P = np.swapaxes(np.linalg.pinv(E), 1, 2)
    H = np.zeros((len(touched),3,len(free)))
    for k, factor in [(0,-P.sum(2)), (1,P[:,:,0]), (2,P[:,:,1])]:
        for j, raw in enumerate(f[touched,k]):
            if index[raw] >= 0: H[j,:,index[raw]] += factor[j]
    C=[]; rhs=[]
    for a in first:
        for b in second:
            row=np.zeros(len(free))
            if index[a]>=0: row[index[a]]+=1
            if index[b]>=0: row[index[b]]-=1
            C.append(row);rhs.append(gap_m*scale-(x[a]-x[b])@n)
    C=np.asarray(C);rhs=np.asarray(rhs)
    def deadline(_=None):
        if time.monotonic()-start > wall_seconds:
            raise TimeoutError('Scalar patch native deadline')
    def objective(u):
        deadline();return .5*float(u@u)
    def metric(u):
        g=np.einsum('fij,j->fi',H,u)
        return maximum_gradient**2-np.einsum('fi,fi->f',g,g)
    def metric_jac(u):
        g=np.einsum('fij,j->fi',H,u)
        return -2*np.einsum('fi,fij->fj',g,H)
    r=minimize(objective,np.zeros(len(free)),jac=lambda u:u,method='SLSQP',
        bounds=[(-maximum_motion_m*scale,maximum_motion_m*scale)]*len(free),
        constraints=[dict(type='ineq',fun=lambda u:C@u-rhs,jac=lambda u:C),
                     dict(type='ineq',fun=metric,jac=metric_jac)],
        callback=deadline,options=dict(maxiter=100,ftol=1e-12,disp=False))
    deadline()
    q=p.copy();q[free]+=r.x[:,None]*n/scale
    gradient=float(np.linalg.norm(np.einsum('fij,j->fi',H,r.x),axis=1).max())
    gap=float(((q[first]-p[first[0]])@n).min()-((q[second]-p[first[0]])@n).max())
    proof=dict(success=bool(r.success),message=str(r.message),iterations=int(r.nit),
        maximumGradient=gradient,projectedGap_m=gap,
        maximumMotion_m=float(np.abs(r.x).max()/scale),
        elapsedSeconds=time.monotonic()-start,anatomicalAcceptance=False)
    if not r.success or gradient>maximum_gradient*(1+1e-7) or gap<gap_m-1e-12:
        raise ValueError('Inadmissible metric-bounded patch: '+str(proof))
    return q,proof
