"""Bounded vertex clearance using exact linear corner-field constraints.

With two triangle corners fixed, every outward-normal inequality is linear
in the remaining corner. Add separating planes from actual collisions, solve
a minimum-maximum coordinate move, then test neighboring Float32 points on
the literal geometry with the global crossing guard. No tissue is inferred.
"""
import itertools,time
import numpy as np
from scipy.optimize import linprog
from geometry_tools.algorithms.local_crossing_guard import CrossingBudgetExceeded


def repair_vertex_clearance(desired_positions,stored_positions,faces,corner_words,*,origin,
        global_guard,maximum_displacement_m,source_clusters=None,locked=(),
        wall_seconds=12.,maximum_lp_solves=100,maximum_pair_constraints=16):
    started=time.monotonic();v=np.asarray(stored_positions,float).copy();desired=np.asarray(desired_positions,float)
    f=np.asarray(faces,int);w=np.asarray(corner_words,float);origin=np.asarray(origin,float)
    if desired.shape!=v.shape or v.ndim!=2 or v.shape[1]!=3 or w.shape!=(len(f),3,3)or origin.shape!=(3,)or not np.isfinite(np.r_[v.ravel(),desired.ravel(),w.ravel(),origin]).all():raise ValueError('Finite corresponding stored/source controls, triangle fields and origin required')
    if not 0<maximum_displacement_m or not np.isfinite(maximum_displacement_m)or not 0<wall_seconds or not np.isfinite(wall_seconds)or not isinstance(maximum_lp_solves,int)or maximum_lp_solves<1 or not isinstance(maximum_pair_constraints,int)or not 1<=maximum_pair_constraints<=100:raise ValueError('Explicit finite displacement, time and solve budgets required')
    if not np.array_equal(global_guard.positions,v)or not np.array_equal(global_guard.faces,f):raise ValueError('Global guard must bind the actual source')
    relative=(v-origin).astype(np.float32)
    if not np.array_equal(relative.astype(float)+origin,v):raise ValueError('Stored points must use the declared Float32 grid')
    clusters={int(k):np.asarray(members,int)for k,members in(source_clusters or {}).items()};fixed=set(map(int,locked))
    if np.linalg.norm(v-desired,axis=1).max()>maximum_displacement_m:raise ValueError('Initial source displacement exceeds the bound')
    for k,members in clusters.items():
        if k not in members or members.min()<0 or members.max()>=len(v)or np.linalg.norm(desired[members]-v[k],axis=1).max()>maximum_displacement_m:raise ValueError('Source alias cluster exceeds its displacement bound')
    order=np.argsort(f.ravel(),kind='stable');counts=np.bincount(f.ravel(),minlength=len(v));starts=np.r_[0,np.cumsum(counts)];incident=order//3
    used=np.unique(f);occupancy={}
    for vertex in used:occupancy.setdefault(tuple(v[vertex]),set()).add(int(vertex))
    initial=len(global_guard.remaining);moves=[];solves=0;stop=None

    def separating(pair,vertex):
        a,b=map(int,pair)
        if vertex not in f[a]:a,b=b,a
        if vertex not in f[a]or vertex in f[b]:return []
        other=v[f[b]];normal=np.cross(other[1]-other[0],other[2]-other[0]);normal/=np.linalg.norm(normal)
        local=v[f[a]];d=(local-other[0])@normal;j=int(np.flatnonzero(f[a]==vertex)[0]);remain=np.delete(d,j);tol=128*np.finfo(float).eps*max(1.,float(np.abs(local).max()))
        if np.all(remain>=-tol)and np.any(remain>tol):signs=[1.]
        elif np.all(remain<=tol)and np.any(remain<-tol):signs=[-1.]
        elif np.all(np.abs(remain)<=tol):signs=[1.,-1.]
        else:return []
        margin=max(1e-12,.55*float(np.abs(normal)@np.abs(np.spacing(relative[vertex])).astype(float)))
        return [(sign*normal,(-sign*d[j]+margin)*1e9)for sign in signs]

    while global_guard.remaining:
        if time.monotonic()-started>=wall_seconds:stop='native wall deadline';break
        accepted=False
        for pair in sorted(global_guard.remaining):
            for vertex in np.unique(f[list(pair)]):
                vertex=int(vertex)
                if vertex in fixed:continue
                ids=np.unique(incident[starts[vertex]:starts[vertex+1]]);p=v[f[ids]];normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);corner=np.argmax(f[ids]==vertex,axis=1)
                edge=p[np.arange(len(ids)),(corner+2)%3]-p[np.arange(len(ids)),(corner+1)%3]
                coeff=np.cross(w[ids],edge[:,None]);length=np.linalg.norm(coeff,axis=2)
                if np.any(length==0):continue
                A=(coeff/length[:,:,None]).reshape(-1,3);distance=(np.einsum('fki,fi->fk',w[ids],normal)/length).reshape(-1)
                if np.any(distance<=0):raise ValueError('Current outward field is invalid')
                B=-distance*1e9+.001
                members=clusters.get(vertex,np.array([vertex]));lo=np.max(desired[members]-v[vertex]-maximum_displacement_m,axis=0)*1e9;hi=np.min(desired[members]-v[vertex]+maximum_displacement_m,axis=0)*1e9
                for first in separating(pair,vertex):
                    constraints=[first];seen={tuple(map(int,pair))}
                    for _ in range(maximum_pair_constraints):
                        if solves>=maximum_lp_solves:stop='LP solve budget';break
                        remaining=wall_seconds-(time.monotonic()-started)
                        if remaining<=0:stop='native wall deadline';break
                        planes=np.r_[A,np.array([a for a,b in constraints])];rhs=np.r_[B,np.array([b for a,b in constraints])]
                        ub=np.c_[-planes,np.zeros(len(planes))];limit=-rhs
                        for axis in range(3):
                            for sign in [-1,1]:
                                row=np.zeros(4);row[axis]=sign;row[3]=-1;ub=np.r_[ub,[row]];limit=np.r_[limit,0.]
                        result=linprog([0,0,0,1],A_ub=ub,b_ub=limit,bounds=[*zip(lo,hi),(0,2*maximum_displacement_m*1e9)],method='highs',options={'time_limit':min(.1,remaining),'primal_feasibility_tolerance':1e-9,'dual_feasibility_tolerance':1e-9});solves+=1
                        if not result.success:
                            if result.status==1:stop='LP native deadline'
                            break
                        point=(v[vertex]+result.x[:3]*1e-9-origin).astype(np.float32)
                        choices=[[np.nextafter(x,np.float32(-np.inf)),x,np.nextafter(x,np.float32(np.inf))]for x in point]
                        candidates=sorted(itertools.product(*choices),key=lambda q:float(np.linalg.norm(np.asarray(q,float)-relative[vertex])))
                        new_pairs=set()
                        for q in candidates:
                            if time.monotonic()-started>=wall_seconds:stop='native wall deadline';break
                            q=np.asarray(q,np.float32);value=q.astype(float)+origin
                            if np.array_equal(value,v[vertex])or occupancy.get(tuple(value),set())-{vertex}or np.linalg.norm(desired[members]-value,axis=1).max()>maximum_displacement_m:continue
                            patch=p.copy();patch[f[ids]==vertex]=value;n=np.cross(patch[:,1]-patch[:,0],patch[:,2]-patch[:,0])
                            if np.any(np.einsum('fki,fi->fk',w[ids],n)<=0):continue
                            try:token=global_guard.check_move([vertex],value[None])
                            except CrossingBudgetExceeded:stop='global pair-test budget';break
                            if token is None:new_pairs.update(global_guard.last_new_pairs);continue
                            if len(token.remaining)>=len(global_guard.remaining):continue
                            global_guard.commit(token);occupancy[tuple(v[vertex])].remove(vertex);occupancy.setdefault(tuple(value),set()).add(vertex);v[vertex]=value;relative[vertex]=q
                            moves.append(dict(vertex=vertex,remaining_global_defect_pairs=len(global_guard.remaining),maximum_cluster_source_displacement_m=float(np.linalg.norm(desired[members]-value,axis=1).max())));accepted=True;break
                        if accepted or stop:break
                        added=False
                        for new_pair in sorted(new_pairs-seen):
                            options=separating(new_pair,vertex)
                            if len(options)==1:constraints.append(options[0]);seen.add(new_pair);added=True
                        if not added:break
                    if accepted or stop:break
                if accepted or stop:break
            if accepted or stop:break
        if stop:break
        if not accepted:stop='no improving bounded linear clearance proposal';break
    return dict(positions=v,evidence=dict(passed=not global_guard.remaining and stop is None,initial_global_defect_pairs=initial,remaining_global_defect_pair_count=len(global_guard.remaining),moves=moves,LP_solves=solves,maximum_LP_solves=maximum_lp_solves,native_wall_seconds=wall_seconds,budget_or_stop_reason=stop,faces_words_uvs_unchanged=True,global_crossing_guard=global_guard.evidence(),elapsed_seconds=time.monotonic()-started))
