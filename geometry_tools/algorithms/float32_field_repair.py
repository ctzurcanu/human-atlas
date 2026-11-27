"""Bounded common-grid rounding repair for an outward corner-normal field.

All cells must use this same position bank. Faces and attributes are immutable.
The result establishes only outward corner fields and distinct moved controls;
topology, intersections, occupied volumes and anatomy require separate audits.
"""
import itertools
import time
import numpy as np


def repair_outward_rounding(positions,faces,corner_words,*,origin=(0.,0.,0.),
        maximum_displacement_m,max_evaluations=5000,maximum_rounds=40,wall_seconds=8.,locked=(),
        selection='global-best',candidate_mode='local-neighbours',
        require_bad_face_reduction=False):
    start=time.monotonic()
    desired=np.asarray(positions,float)
    f=np.asarray(faces,int)
    words=np.asarray(corner_words,float)
    origin=np.asarray(origin,float)
    if selection not in ('global-best','first-improving'):
        raise ValueError('Explicit global-best or first-improving selection required')
    if candidate_mode not in ('local-neighbours','bounded-global-lattice'):
        raise ValueError('Explicit local-neighbours or bounded-global-lattice candidates required')
    if words.shape!=(len(f),3,3)or np.any(np.linalg.norm(words,axis=2)==0):
        raise ValueError('A nonzero outward normal is required at every source corner')
    if not np.isfinite(maximum_displacement_m)or maximum_displacement_m<=0:
        raise ValueError('A finite positive numerical displacement bound is required')
    relative=(desired-origin).astype(np.float32)
    initial_relative=relative.copy()
    world=relative.astype(float)+origin
    baseline=world.copy()
    if np.any(np.linalg.norm(world-desired,axis=1)>maximum_displacement_m):
        raise ValueError('The declared bound is smaller than initial Float32 rounding')
    used=np.unique(f)
    occupancy={}
    for i in used:occupancy.setdefault(tuple(world[i]),set()).add(int(i))
    order=np.argsort(f.ravel(),kind='stable')
    counts=np.bincount(f.ravel(),minlength=len(world))
    offsets=np.r_[0,np.cumsum(counts)]
    incident=order//3
    word_length=np.linalg.norm(words,axis=2)
    protected=set(map(int,locked))
    quantum=float(np.max(np.abs(np.spacing(np.float32(np.max(np.abs(relative[used]),axis=0))))))
    lattice=None
    if candidate_mode=='bounded-global-lattice':
        ratio=maximum_displacement_m/quantum
        if not np.isfinite(ratio)or ratio>7:
            raise ValueError('Bounded global lattice requires at most seven global quanta')
        radius=int(np.ceil(ratio))+1
        lattice=np.asarray(list(itertools.product(range(-radius,radius+1),repeat=3)),float)*quantum
    def margins(ids,vertex=None,value=None):
        p=world[f[ids]].copy()
        if vertex is not None:p[f[ids]==vertex]=value
        n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
        length=np.linalg.norm(n,axis=1)
        return np.einsum('tki,ti->tk',words[ids],n)/np.maximum(length[:,None]*word_length[ids],1e-300)
    margin=margins(np.arange(len(f)))
    bad=np.any(margin<=0,axis=1)
    initial=int(bad.sum())
    evaluations=0;moves=[];reason=None
    for iteration in range(maximum_rounds):
        if not np.any(bad):break
        best=None
        for vertex in np.unique(f[bad]):
            vertex=int(vertex)
            if vertex in protected:continue
            ids=np.unique(incident[offsets[vertex]:offsets[vertex+1]])
            old_bad=bad[ids]
            old_score=(int(old_bad.sum()),float(np.maximum(-margin[ids],0).sum()))
            point=relative[vertex].copy()
            choices=[]
            for axis in range(3):
                down=np.nextafter(point[axis],np.float32(-np.inf))
                up=np.nextafter(point[axis],np.float32(np.inf))
                choices.append([down,point[axis],up])
            candidates=list(itertools.product(*choices))
            if lattice is not None:
                coarse=(initial_relative[vertex].astype(float)+lattice).astype(np.float32)
                displacement=np.linalg.norm(coarse.astype(float)+origin-desired[vertex],axis=1)
                candidates.extend(map(tuple,coarse[displacement<=maximum_displacement_m]))
            for axis in range(3):
                for sign in (-1,1):
                    q=point.copy();q[axis]=np.float32(float(q[axis])+sign*quantum)
                    candidates.append(tuple(q))
            seen=set()
            for value in candidates:
                if evaluations>=max_evaluations:reason='evaluation budget';break
                if time.monotonic()-start>=wall_seconds:reason='native wall deadline';break
                q=np.asarray(value,np.float32)
                key=tuple(map(float,q))
                if key in seen:continue
                seen.add(key)
                candidate=q.astype(float)+origin
                if np.array_equal(candidate,world[vertex]):continue
                movement=float(np.linalg.norm(candidate-desired[vertex]))
                if movement>maximum_displacement_m:continue
                occupied=occupancy.get(tuple(candidate),set())
                if occupied-{vertex}:continue
                evaluations+=1
                current=margins(ids,vertex,candidate)
                current_bad=np.any(current<=0,axis=1)
                if np.any(current_bad&~old_bad):continue
                score=(int(current_bad.sum()),float(np.maximum(-current,0).sum()))
                if require_bad_face_reduction and score[0]>=old_score[0]:continue
                if not score<old_score:continue
                rank=score+(movement,vertex)
                if best is None or rank<best[0]:
                    best=(rank,vertex,q,candidate,ids,current,current_bad)
            if reason:break
            if selection=='first-improving' and best is not None:break
        if best is None:
            reason=reason or 'no improving bounded grid move'
            break
        _,vertex,q,candidate,ids,current,current_bad=best
        old=tuple(world[vertex]);occupancy[old].remove(vertex)
        if not occupancy[old]:occupancy.pop(old)
        occupancy.setdefault(tuple(candidate),set()).add(vertex)
        relative[vertex]=q;world[vertex]=candidate
        margin[ids]=current;bad[ids]=current_bad
        moves.append(dict(vertex=vertex,remaining_bad_faces=int(bad.sum()),
            displacement_from_requested_m=float(np.linalg.norm(candidate-desired[vertex]))))
        if reason:break
    else:
        reason='round budget'
    if time.monotonic()-start>=wall_seconds:reason='native wall deadline'
    passed=not np.any(bad)and reason!='native wall deadline'
    return dict(positions=world,relative_positions=relative,evidence=dict(passed=passed,
        initial_bad_faces=initial,remaining_bad_faces=int(bad.sum()),
        remaining_bad_face_indices=np.flatnonzero(bad).tolist(),moves=moves,
        maximum_displacement_m=float(np.linalg.norm(world-desired,axis=1).max()),
        maximum_change_from_nearest_float32_m=float(np.linalg.norm(world-baseline,axis=1).max()),
        declared_displacement_bound_m=float(maximum_displacement_m),evaluations=evaluations,
        search_selection=selection,
        candidate_mode=candidate_mode,require_bad_face_reduction=bool(require_bad_face_reduction),
        max_evaluations=max_evaluations,native_wall_seconds=wall_seconds,
        budget_or_stop_reason=None if passed else reason,
        faces_unchanged=True,normal_words_unchanged=True,
        geometry_and_anatomical_audits_pending=True,elapsed_seconds=time.monotonic()-start))
