"""Bounded Newton inversion of a positive, fixed-boundary regular cage."""
import numpy as np


def inverse_cage_points(cage,points,displacement,*,maximum_iterations=24,
                        maximum_backtracks=12,residual_m=2e-11):
    target=np.asarray(points,float);d=np.asarray(displacement,float)
    if target.ndim!=2 or target.shape[1]!=3 or not np.isfinite(target).all():
        raise ValueError('Finite target points are required')
    if not isinstance(maximum_iterations,int)or maximum_iterations<1 or not isinstance(maximum_backtracks,int)or maximum_backtracks<1 or not np.isfinite(residual_m)or residual_m<=0:
        raise ValueError('Explicit finite inverse iteration and residual budgets required')
    proof=cage.check(d)
    if not proof['passed']:raise ValueError('A positive fixed-boundary cage is required')
    jacobians=cage.jacobians(d)
    def mapped(p):
        where=cage.locate(p)
        return p+np.einsum('nk,nki->ni',where['weights'],d[where['node_ids']]),where
    estimate=target-(mapped(target)[0]-target);iterations=0
    for iteration in range(maximum_iterations):
        result,locations=mapped(estimate);residual=result-target;length=np.linalg.norm(residual,axis=1)
        active=np.flatnonzero(length>residual_m)
        if not len(active):break
        j=np.broadcast_to(np.eye(3),(len(active),3,3)).copy();inside=locations['inside'][active]
        j[inside]=jacobians[locations['tetrahedron'][active[inside]]]
        step=np.linalg.solve(j,residual[active,:,None])[...,0];pending=np.arange(len(active));alpha=np.ones(len(active))
        for _ in range(maximum_backtracks):
            trial=estimate[active[pending]]-alpha[pending,None]*step[pending]
            error=np.linalg.norm(mapped(trial)[0]-target[active[pending]],axis=1)
            good=(error<length[active[pending]])|(error<=residual_m)
            estimate[active[pending[good]]]=trial[good];pending=pending[~good]
            if not len(pending):break
            alpha[pending]*=.5
        if len(pending):raise ValueError('Inverse cage line search exhausted its declared budget')
        iterations=iteration+1
    replay,locations=mapped(estimate);error=np.linalg.norm(replay-target,axis=1)
    if np.any(error>residual_m):raise ValueError('Inverse cage residual exceeds the declared tolerance')
    j=np.broadcast_to(np.eye(3),(len(estimate),3,3)).copy();take=locations['inside']
    j[take]=jacobians[locations['tetrahedron'][take]]
    return dict(source_points=estimate,locations=locations,jacobians=j,iterations=iterations,
        maximum_residual_m=float(error.max(initial=0)),positive_cage_gate=proof)
