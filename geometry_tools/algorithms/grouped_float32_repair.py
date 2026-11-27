"""Bounded rigid grid translations of explicitly declared control groups.

This is a numerical candidate generator. Faces, normal words and UVs stay
fixed; every accepted move reduces outward failures or an explicitly replayed
crossing worklist without increasing either. Full geometry still needs review.
"""
from __future__ import annotations
import itertools
import time
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.algorithms.local_crossing_guard import CrossingBudgetExceeded


def repair_grouped_rounding(desired_positions, faces, corner_words, groups, *,
        origin, maximum_displacement_m, stored_positions=None,
        max_evaluations=30000, maximum_rounds=100, wall_seconds=15., locked=(),
        defect_face_pairs=(), translations_m=None,global_guard=None,
        maximum_group_size=8,source_clusters=None,maximum_displacement_by_vertex=None):
    started=time.monotonic()
    desired=np.asarray(desired_positions,float);f=np.asarray(faces,int)
    words=np.asarray(corner_words,float);origin=np.asarray(origin,float)
    if desired.ndim!=2 or desired.shape[1]!=3 or not np.isfinite(desired).all():
        raise ValueError('Finite source controls required')
    if f.ndim!=2 or f.shape[1]!=3 or not len(f) or words.shape!=(len(f),3,3)or not np.isfinite(words).all():
        raise ValueError('Triangle and corresponding corner fields required')
    if origin.shape!=(3,)or not np.isfinite(origin).all():raise ValueError('Finite three-coordinate origin required')
    if not np.isfinite(maximum_displacement_m)or maximum_displacement_m<=0 or not np.isfinite(wall_seconds)or wall_seconds<=0:
        raise ValueError('Positive finite displacement and wall limits required')
    if max_evaluations<0 or maximum_rounds<1:raise ValueError('Finite iteration budgets required')
    relative=(desired-origin).astype(np.float32)if stored_positions is None else(np.asarray(stored_positions,float)-origin).astype(np.float32)
    world=relative.astype(float)+origin
    if global_guard is not None and(not np.array_equal(global_guard.positions,world)or not np.array_equal(global_guard.faces,f)):raise ValueError('Global crossing guard must bind the actual initial coordinates and faces')
    if stored_positions is not None and not np.array_equal(world,np.asarray(stored_positions,float)):
        raise ValueError('Initial stored controls must already use the declared Float32 grid')
    if np.max(np.linalg.norm(world-desired,axis=1))>maximum_displacement_m:
        raise ValueError('Initial encoding exceeds the displacement bound')
    control_bounds=np.full(len(world),maximum_displacement_m)if maximum_displacement_by_vertex is None else np.asarray(maximum_displacement_by_vertex,float)
    if control_bounds.shape!=(len(world),)or not np.isfinite(control_bounds).all()or np.any(control_bounds<0)or np.any(control_bounds>maximum_displacement_m):raise ValueError('Finite per-control bounds within the global bound required')
    if not isinstance(maximum_group_size,int)or not 1<=maximum_group_size<=32:raise ValueError('Explicit group-size bound from one to 32 required')
    clusters={int(k):np.asarray(sorted(set(map(int,members))),int)for k,members in(source_clusters or {}).items()}
    claimed=set()
    for k,members in clusters.items():
        if k not in members or members.min()<0 or members.max()>=len(world)or claimed.intersection(map(int,members))or np.linalg.norm(desired[members]-world[k],axis=1).max()>maximum_displacement_m:raise ValueError('Disjoint original source alias clusters within the displacement bound required')
        claimed.update(map(int,members))
    declared=[];fixed=set(map(int,locked));used=np.unique(f);used_set=set(map(int,used))
    for group in groups:
        group=tuple(sorted(set(map(int,group))))
        if not group or len(group)>maximum_group_size or min(group)<0 or max(group)>=len(world)or not set(group)<=used_set:
            raise ValueError('Explicit groups within the declared control-count bound required')
        if not fixed.intersection(group):declared.append(group)
    declared=list(dict.fromkeys(declared))
    if len(declared)>2000:raise ValueError('Split more than 2000 control groups into saved stages')
    occupancy={}
    for i in used:occupancy.setdefault(tuple(world[i]),set()).add(int(i))
    tree=cKDTree(desired[used]);padding=64*np.finfo(float).eps*max(1.,float(np.abs(world).max()))
    order=np.argsort(f.ravel(),kind='stable');counts=np.bincount(f.ravel(),minlength=len(world));offsets=np.r_[0,np.cumsum(counts)];incident=order//3
    word_length=np.linalg.norm(words,axis=2)
    if np.any(word_length==0):raise ValueError('Nonzero corner normals required')
    quantum=float(np.max(np.abs(np.spacing(np.max(np.abs(relative[used]),axis=0)))))
    if translations_m is None:
        ratio=maximum_displacement_m/quantum
        if not np.isfinite(ratio)or ratio>7:raise ValueError('At most seven global grid quanta are supported')
        radius=int(np.ceil(2*ratio))+1
        delta=np.asarray(list(itertools.product(range(-radius,radius+1),repeat=3)),float)*quantum
    else:
        delta=np.asarray(translations_m,float)
        if delta.ndim!=2 or delta.shape[1]!=3 or not 0<len(delta)<=500 or not np.isfinite(delta).all():raise ValueError('One to 500 explicit finite group translations required')
    delta=delta[np.argsort(np.linalg.norm(delta,axis=1),kind='stable')]
    def margins(ids, proposals=None, group=None):
        p=world[f[ids]]
        if proposals is not None:
            p=np.broadcast_to(p,(len(proposals),)+p.shape).copy()
            for j,vertex in enumerate(group):p[:,f[ids]==vertex,:]=proposals[:,j,None,:]
        n=np.cross(p[...,1,:]-p[...,0,:],p[...,2,:]-p[...,0,:])
        length=np.linalg.norm(n,axis=-1)
        if proposals is None:return np.einsum('fki,fi->fk',words[ids],n)/np.maximum(length[:,None]*word_length[ids],1e-300)
        return np.einsum('fki,cfi->cfk',words[ids],n)/np.maximum(length[:,:,None]*word_length[ids][None],1e-300)
    margin=margins(np.arange(len(f)));bad=np.any(margin<=0,axis=1);initial=int(bad.sum())
    pairs=np.asarray(defect_face_pairs,int).reshape(-1,2)
    if len(pairs)>2000 or len(pairs)and(pairs.min()<0 or pairs.max()>=len(f)or np.any(pairs[:,0]==pairs[:,1])):raise ValueError('At most 2000 explicit distinct source-face pairs required')
    from geometry_tools.distance import triangle_pairs
    def crossings(ids,proposals=None,group=None):
        p=world[f[pairs[ids]]].copy()
        if proposals is None:
            if not len(ids):return np.zeros(0,bool)
            return np.isin(triangle_pairs(p[:,0],p[:,1])['relation'],['crossing','coplanar-overlap'])
        if not len(ids):return np.zeros((len(proposals),0),bool)
        p=np.broadcast_to(p,(len(proposals),)+p.shape).copy()
        for j,vertex in enumerate(group):p[:,f[pairs[ids]]==vertex,:]=proposals[:,j,None,:]
        relation=triangle_pairs(p[:,:,0].reshape(-1,3,3),p[:,:,1].reshape(-1,3,3))['relation']
        return np.isin(relation,['crossing','coplanar-overlap']).reshape(len(proposals),len(ids))
    known_bad=crossings(np.arange(len(pairs)));initial_known_bad=int(known_bad.sum())
    group_faces={g:np.unique(np.concatenate([incident[offsets[v]:offsets[v+1]]for v in g]))for g in declared}
    group_pairs={g:np.flatnonzero(np.isin(f[pairs],g).any(axis=(1,2)))for g in declared}
    blocked=set();moves=[];evaluations=0;skips=0;reason=None
    for _ in range(maximum_rounds):
        if not bad.any()and not known_bad.any():break
        accepted=False
        for group in declared:
            if time.monotonic()-started>=wall_seconds:reason='native wall deadline';break
            if group in blocked:skips+=1;continue
            ids=group_faces[group];old_bad=bad[ids]
            pids=group_pairs[group];old_known=known_bad[pids]
            if not old_bad.any()and not old_known.any():blocked.add(group);continue
            old=world[list(group)].copy();candidate_rel=(relative[list(group)].astype(float)[None]+delta[:,None,:]).astype(np.float32)
            candidate=candidate_rel.astype(float)+origin
            displacement=np.linalg.norm(candidate-desired[list(group)][None],axis=2)
            shift=candidate-old[None]
            rigid=np.max(np.abs(shift-shift[:,:1]),axis=(1,2))<=padding
            eligible=np.all(displacement<=control_bounds[list(group)][None],axis=1)&rigid&np.any(candidate!=old[None],axis=(1,2))
            if global_guard is not None:
                eligible&=np.all(np.linalg.norm(candidate-global_guard.source[list(group)][None],axis=2)<=global_guard.maximum,axis=1)
            for j,vertex in enumerate(group):
                if vertex in clusters:eligible&=np.linalg.norm(candidate[:,j,None,:]-desired[clusters[vertex]][None],axis=2).max(1)<=maximum_displacement_m
            candidate,candidate_rel,displacement=candidate[eligible],candidate_rel[eligible],displacement[eligible]
            good=[];members=set(group)
            for i,points in enumerate(candidate):
                if any(occupancy.get(tuple(point),set())-members for point in points):continue
                if len({tuple(point)for point in points})!=len(group):continue
                good.append(i)
            candidate,candidate_rel,displacement=candidate[good],candidate_rel[good],displacement[good]
            remaining=max_evaluations-evaluations
            if remaining<=0:reason='evaluation budget';break
            candidate,candidate_rel,displacement=candidate[:remaining],candidate_rel[:remaining],displacement[:remaining]
            if not len(candidate):blocked.add(group);continue
            current=margins(ids,candidate,group);evaluations+=len(candidate);current_bad=np.any(current<=0,axis=2)
            allowed=~np.any(current_bad&~old_bad[None],axis=1)
            count=current_bad.sum(1);allowed&=count<=int(old_bad.sum())
            known=crossings(pids,candidate,group);known_count=known.sum(1)
            allowed&=~np.any(known&~old_known[None],axis=1)&(known_count<=int(old_known.sum()))
            allowed&=(count<int(old_bad.sum()))|(known_count<int(old_known.sum()))
            possible=np.flatnonzero(allowed)
            if not len(possible):blocked.add(group);continue
            ranked=sorted(possible,key=lambda i:(int(count[i]),int(known_count[i]),float(displacement[i].max()),float(displacement[i].sum())))
            chosen=None;guarded=None
            for proposed in ranked:
                if time.monotonic()-started>=wall_seconds:reason='native wall deadline';break
                if global_guard is not None:
                    try:guarded=global_guard.check_move(group,candidate[proposed])
                    except CrossingBudgetExceeded:reason='global pair-test budget';break
                    if guarded is None:continue
                chosen=proposed;break
            if chosen is None:
                if reason:break
                blocked.add(group);continue
            if global_guard is not None:global_guard.commit(guarded)
            for vertex,point in zip(group,old):
                occupancy[tuple(point)].remove(vertex)
                if not occupancy[tuple(point)]:occupancy.pop(tuple(point))
            for vertex,point,rel in zip(group,candidate[chosen],candidate_rel[chosen]):
                world[vertex]=point;relative[vertex]=rel;occupancy.setdefault(tuple(point),set()).add(vertex)
            margin[ids]=current[chosen];bad[ids]=current_bad[chosen]
            known_bad[pids]=known[chosen]
            affected=set(map(int,np.unique(f[ids])))|set(map(int,np.unique(f[pairs[pids]])))
            for point in np.r_[old,candidate[chosen]]:
                affected.update(map(int,used[tree.query_ball_point(point,maximum_displacement_m+padding)]))
            blocked={g for g in blocked if not affected.intersection(g)}
            moves.append(dict(vertices=list(group),remaining_bad_faces=int(bad.sum()),
                remaining_declared_defect_pairs=int(known_bad.sum()),
                maximum_displacement_from_source_m=float(displacement[chosen].max()),
                maximum_rigid_translation_residual_m=float(np.max(np.abs((candidate[chosen]-old)-(candidate[chosen,:1]-old[:1]))))))
            accepted=True;break
        if reason:break
        if not accepted:reason='no improving bounded rigid group translation';break
    else:reason='round budget'
    if time.monotonic()-started>=wall_seconds:reason='native wall deadline'
    passed=not bad.any()and not known_bad.any()and reason not in ['native wall deadline','global pair-test budget']
    aggregate=max([float(np.linalg.norm(world-desired,axis=1).max())]+[float(np.linalg.norm(desired[members]-world[k],axis=1).max())for k,members in clusters.items()])
    return dict(positions=world,relative_positions=relative,evidence=dict(passed=passed,
        initial_bad_faces=initial,remaining_bad_faces=int(bad.sum()),
        remaining_bad_face_indices=np.flatnonzero(bad).tolist(),moves=moves,
        initial_declared_defect_pairs=initial_known_bad,
        remaining_declared_defect_pairs=pairs[known_bad].tolist(),
        remaining_declared_defect_pair_count=int(known_bad.sum()),
        declared_groups=len(declared),cached_group_skips=skips,evaluations=evaluations,
        max_evaluations=max_evaluations,native_wall_seconds=wall_seconds,
        maximum_displacement_m=float(np.linalg.norm(world-desired,axis=1).max()),
        declared_displacement_bound_m=float(maximum_displacement_m),
        maximum_group_size=maximum_group_size,source_alias_clusters=len(clusters),
        maximum_aggregate_source_displacement_m=aggregate,
        faces_unchanged=True,normal_words_unchanged=True,
        global_crossing_guard=global_guard.evidence()if global_guard is not None else None,
        budget_or_stop_reason=None if passed else reason,
        geometry_and_anatomical_audits_pending=True,elapsed_seconds=time.monotonic()-started))
