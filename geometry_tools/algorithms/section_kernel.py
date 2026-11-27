"""Bounded planar one-ring untangling with literal fixed anatomical controls.

For one free vertex, every incident oriented-area constraint is affine. Their
half-plane intersection is its feasible kernel. Moving to that kernel cannot
invert an incident triangle; all other triangles remain unchanged. The caller
must still audit the complete 3D construction and anatomical interpretation.
"""
import time
import numpy as np
from scipy.optimize import minimize


def repair_section_kernel(points,faces,fixed,*,minimum_double_area,
                          maximum_displacement_m,maximum_steps=64,wall_seconds=2.):
    start=time.monotonic();original=np.asarray(points,float);p=original.copy();f=np.asarray(faces,int)
    floor=np.broadcast_to(np.asarray(minimum_double_area,float),(len(f),))
    if np.any(floor<=0)or maximum_displacement_m<=0 or maximum_steps<1 or wall_seconds<=0:
        raise ValueError('Positive finite area, displacement, step and wall limits required')
    fixed=set(map(int,fixed));incident=[[]for _ in p]
    for i,tri in enumerate(f):
        for vertex in tri:incident[int(vertex)].append(i)
    def area(q):
        t=q[f];a=t[:,1]-t[:,0];b=t[:,2]-t[:,0]
        return a[:,0]*b[:,1]-a[:,1]*b[:,0]
    def clip(poly,a,b):
        result=[];value=poly@a-b
        for i,j in zip(range(len(poly)),np.roll(np.arange(len(poly)),-1)):
            if value[i]>=0:result.append(poly[i])
            if(value[i]>=0)!=(value[j]>=0):
                w=value[i]/(value[i]-value[j]);result.append(poly[i]*(1-w)+poly[j]*w)
        return np.asarray(result)
    current=area(p);initial=int(np.count_nonzero(current<floor));moves=[];stop=None
    for _ in range(maximum_steps):
        bad=np.flatnonzero(current<floor)
        if not len(bad):break
        if time.monotonic()-start>=wall_seconds:stop='native wall deadline';break
        vertices,count=np.unique(f[bad],return_counts=True);order=vertices[np.argsort(-count,kind='stable')]
        accepted=False
        for vertex in order:
            vertex=int(vertex)
            if vertex in fixed:continue
            theta=np.arange(32)*2*np.pi/32
            kernel=maximum_displacement_m*np.c_[np.cos(theta),np.sin(theta)]
            ids=np.asarray(incident[vertex]);tri=f[ids];a=[];b=[]
            shift=p[vertex]-original[vertex]
            for index,row in zip(ids,tri):
                k=int(np.flatnonzero(row==vertex)[0]);one=p[row[(k+1)%3]];two=p[row[(k+2)%3]]
                gradient=np.array([one[1]-two[1],two[0]-one[0]])
                a.append(gradient);b.append(floor[index]-current[index]+gradient@shift)
            for normal,offset in zip(a,b):
                kernel=clip(kernel,normal,offset)
                if len(kernel)<3:break
            if len(kernel)<3:continue
            # The convex kernel is inside the explicit Euclidean source ball.
            candidates=[]
            for one,two in zip(kernel,np.roll(kernel,-1,axis=0)):
                edge=two-one;den=float(edge@edge)
                w=0. if den==0 else np.clip(-(one@edge)/den,0,1)
                candidates.append(one+w*edge)
            candidate=min(candidates,key=lambda q:q@q)
            candidate=.999*candidate+.001*kernel.mean(0)
            proposed=p[vertex].copy();p[vertex]=original[vertex]+candidate
            updated=area(p)
            if np.any(updated[ids]<floor[ids]*(1-1e-9))or np.count_nonzero(updated<floor)>=len(bad):
                p[vertex]=proposed;continue
            moves.append(dict(vertex=vertex,displacement_m=float(np.linalg.norm(candidate)),
                incident_faces=len(ids),remaining_bad_faces=int(np.count_nonzero(updated<floor))))
            current=updated;accepted=True;break
        if not accepted:stop='no bounded free-vertex kernel';break
    else:stop='step budget'
    bad=np.flatnonzero(current<floor)
    proof=dict(passed=not len(bad),initial_bad_faces=initial,remaining_bad_faces=bad.tolist(),
        moves=moves,maximum_displacement_m=float(np.linalg.norm(p-original,axis=1).max(initial=0)),
        declared_maximum_displacement_m=maximum_displacement_m,
        fixed_controls_literal=bool(np.array_equal(p[sorted(fixed)],original[sorted(fixed)])),
        minimum_double_area_m2=float(current.min()),stop_reason=stop,
        maximum_steps=maximum_steps,native_wall_seconds=wall_seconds,elapsed_seconds=time.monotonic()-start)
    return p,proof


def repair_section_patch(points,faces,fixed,*,minimum_double_area,
                         maximum_displacement_m,initial=None,maximum_steps=80,wall_seconds=3.):
    """Couple a small free patch when individual vertex kernels are empty.

    Signed triangle areas and Euclidean source balls are explicit nonlinear
    constraints. Analytic gradients avoid geometry-wide finite differencing.
    No fixed or out-of-patch control is an optimization variable.
    """
    started=time.monotonic();source=np.asarray(points,float);p=source.copy()if initial is None else np.array(initial,float,copy=True)
    f=np.asarray(faces,int);floor=np.broadcast_to(np.asarray(minimum_double_area,float),(len(f),));fixed=set(map(int,fixed))
    if maximum_displacement_m<=0 or np.any(floor<=0):raise ValueError('Positive finite source and area bounds required')
    if np.linalg.norm(p-source,axis=1).max(initial=0)>maximum_displacement_m:
        raise ValueError('Initial patch already exceeds the source displacement bound')
    def areas(q,tri=f):
        t=q[tri];a=t[:,1]-t[:,0];b=t[:,2]-t[:,0];return a[:,0]*b[:,1]-a[:,1]*b[:,0]
    bad=np.flatnonzero(areas(p)<floor);selected=set(map(int,np.unique(f[bad])))-fixed
    touched=np.flatnonzero(np.isin(f,list(selected)).any(1));selected|=set(map(int,np.unique(f[touched])))-fixed
    variables=np.array(sorted(selected),int)
    if not len(variables)or len(variables)>64:
        return source,dict(passed=False,stop_reason='No bounded patch of at most 64 free controls',moves=[])
    ids=np.flatnonzero(np.isin(f,variables).any(1));tri=f[ids];scale=np.maximum(abs(areas(source)[ids]),floor[ids]);bound=maximum_displacement_m
    lookup={int(v):i for i,v in enumerate(variables)};entries=[(j,k,lookup[int(v)])for j,row in enumerate(tri)for k,v in enumerate(row)if int(v)in lookup]
    def unpack(x):
        if time.monotonic()-started>=wall_seconds:raise TimeoutError('native patch wall deadline')
        q=p.copy();q[variables]=source[variables]+x.reshape(-1,2)*bound;return q
    def constraint(x):
        q=unpack(x);u=x.reshape(-1,2)
        return np.r_[(areas(q,tri)-floor[ids])/scale,1-np.einsum('ij,ij->i',u,u)]
    def jacobian(x):
        q=unpack(x);t=q[tri];result=np.zeros((len(ids)+len(variables),len(x)))
        for j,k,v in entries:
            one,two=t[j,(k+1)%3],t[j,(k+2)%3]
            result[j,2*v:2*v+2]=np.array([one[1]-two[1],two[0]-one[0]])*bound/scale[j]
        for i,u in enumerate(x.reshape(-1,2)):result[len(ids)+i,2*i:2*i+2]=-2*u
        return result
    x0=((p[variables]-source[variables])/bound).ravel();reason=None;solver=None
    try:
        solver=minimize(lambda x:float(x@x/2),x0,jac=lambda x:x,method='SLSQP',
            constraints=[dict(type='ineq',fun=constraint,jac=jacobian)],
            options=dict(maxiter=maximum_steps,ftol=1e-11,disp=False))
        q=unpack(solver.x);after=areas(q)
        distance=float(np.linalg.norm(q-source,axis=1).max())
        passed=bool(solver.success and np.all(after>=floor*(1-1e-8))and distance<=bound*(1+1e-10))
        reason=str(solver.message)
    except TimeoutError as error:
        q=source.copy();after=areas(q);distance=0.;passed=False;reason=str(error)
    if not passed:q=source.copy()
    return q,dict(passed=passed,method='coupled-oriented-area-patch',stop_reason=reason,
        free_controls=variables.tolist(),patch_faces=len(ids),initial_bad_faces=len(bad),
        remaining_bad_faces=np.flatnonzero(after<floor*(1-1e-8)).tolist(),
        maximum_displacement_m=distance,declared_maximum_displacement_m=bound,
        fixed_controls_literal=bool(np.array_equal(q[sorted(fixed)],source[sorted(fixed)])),
        minimum_double_area_m2=float(after.min()),maximum_steps=maximum_steps,
        iterations=int(solver.nit)if solver is not None else 0,native_wall_seconds=wall_seconds,
        elapsed_seconds=time.monotonic()-started)
