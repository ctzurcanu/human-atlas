"""Small source-core proposals with explicit fixed annular controls.

SAT planes are only candidate generators. Actual source self/contact checks,
normal reconstruction and anatomical review remain the caller's responsibility.
"""
import time
import numpy as np
from scipy.optimize import minimize


def core_separation_candidates(vertices, faces, pair, movable, *, gap_m,
                               maximum_motion_m, wall_seconds=3., barriers=()):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);pair=np.asarray(pair,int)
    free=np.asarray(movable,int)
    if (v.ndim!=2 or v.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3
        or pair.shape!=(2,) or np.any(pair<0) or np.any(pair>=len(f))
        or np.any(f<0) or np.any(f>=len(v)) or not len(free) or len(set(free))!=len(free)
        or np.any(free<0) or np.any(free>=len(v)) or not np.isfinite(v).all()
        or not np.isfinite([gap_m,maximum_motion_m,wall_seconds]).all()
        or not 0<gap_m<=maximum_motion_m or not 0<wall_seconds<=5):
        raise ValueError('Finite source triangles, explicit free controls and bounded positive tolerances required')
    barrier_matrix=[];barrier_rhs=[]
    lookup={int(x):i for i,x in enumerate(free)}
    for barrier in barriers:
        face=int(barrier['face']);target=np.asarray(barrier['triangle'],float)
        n=np.asarray(barrier['axis'],float);gap=float(barrier['minimum_gap_m'])
        if (not 0<=face<len(f) or target.shape!=(3,3) or n.shape!=(3,)
            or not np.isfinite(np.r_[target.ravel(),n,gap]).all()
            or gap<0 or abs(np.linalg.norm(n)-1)>1e-10):
            raise ValueError('Explicit fixed-recipient triangle and oriented unit barrier required')
        original=float((v[f[face]]@n).min()-(target@n).max())
        if original<gap-1e-12:raise ValueError('Barrier must preserve an originally separated source interval')
        for x in f[face]:
            for y in target:
                row=np.zeros((len(free),3))
                if int(x)in lookup:row[lookup[int(x)]]+=n
                barrier_matrix.append(row.ravel());barrier_rhs.append((gap-n@(v[x]-y))*1000)
    tick=time.monotonic();a,b=v[f[pair]]
    ea=np.roll(a,-1,axis=0)-a;eb=np.roll(b,-1,axis=0)-b
    normals=np.array([np.cross(ea[0],ea[1]),np.cross(eb[0],eb[1])])
    axes=np.r_[normals,np.cross(ea[:,None],eb[None]).reshape(-1,3),np.cross(normals[0],ea),np.cross(normals[1],eb)]
    axes=axes[np.linalg.norm(axes,axis=1)>1e-25];axes/=np.linalg.norm(axes,axis=1)[:,None]
    unique=[]
    for axis in axes:
        if not any(abs(axis@other)>1-1e-10 for other in unique):unique.append(axis)
    proposals=[]
    for axis in unique:
        for sign in [-1.,1.]:
            if time.monotonic()-tick>=wall_seconds:raise TimeoutError('Source-core plane proposal deadline')
            n=axis*sign;matrix=[];rhs=[]
            for x in f[pair[0]]:
                for y in f[pair[1]]:
                    row=np.zeros((len(free),3))
                    if int(x) in lookup:row[lookup[int(x)]]+=n
                    if int(y) in lookup:row[lookup[int(y)]]-=n
                    matrix.append(row.ravel());rhs.append((gap_m-n@(v[x]-v[y]))*1000)
            matrix.extend(barrier_matrix);rhs.extend(barrier_rhs)
            A=np.array(matrix);rhs=np.array(rhs)
            result=minimize(lambda d:.5*float(d@d),np.zeros(3*len(free)),jac=lambda d:d,
                bounds=[(-maximum_motion_m*1000,maximum_motion_m*1000)]*(3*len(free)),
                constraints=[dict(type='ineq',fun=lambda d:A@d-rhs,jac=lambda d:A)],
                method='SLSQP',options=dict(maxiter=50,ftol=1e-13))
            if not result.success or np.min(A@result.x-rhs)<-1e-7:continue
            delta=result.x.reshape(-1,3)/1000
            if np.linalg.norm(delta,axis=1).max()>maximum_motion_m:continue
            q=v.copy();q[free]+=delta
            minimum=float(np.min(q[f[pair[0]]]@n)-np.max(q[f[pair[1]]]@n))
            if minimum<gap_m-1e-10:continue
            proposals.append((float(np.linalg.norm(delta)),q,dict(axis=n.tolist(),minimumGap_m=minimum,
                movable=free.tolist(),maximumMotion_m=float(np.linalg.norm(delta,axis=1).max()),
                squaredMotion_m2=float(np.sum(delta*delta)))))
    proposals.sort(key=lambda p:p[0])
    return [(q,proof)for _,q,proof in proposals]
