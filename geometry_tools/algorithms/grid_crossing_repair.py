"""Bounded common-Float32-grid proposals for an explicit crossing worklist.

Faces, tissue ownership, UVs and corner fields are immutable. Each accepted
move reduces verified crossings without failing any incident outward field.
The caller must still audit new global crossings, topology and occupied space.
"""
import itertools
import time
import numpy as np
from geometry_tools.distance import triangle_pairs
from geometry_tools.algorithms.local_crossing_guard import CrossingBudgetExceeded


def repair_crossing_rounding(desired_positions, stored_positions, faces, corner_words,
        defect_face_pairs, *, origin, maximum_displacement_m=5e-8,
        source_clusters=None, locked=(), max_evaluations=30000,
        maximum_rounds=100, wall_seconds=15., candidate_mode='local-neighbours',global_guard=None,
        maximum_facet_rotation_degrees=None):
    started=time.monotonic()
    desired=np.asarray(desired_positions,float);world=np.asarray(stored_positions,float).copy()
    f=np.asarray(faces,int);w=np.asarray(corner_words,float);pairs=np.asarray(defect_face_pairs,int).reshape(-1,2);origin=np.asarray(origin,float)
    if desired.shape!=world.shape or world.ndim!=2 or world.shape[1]!=3 or origin.shape!=(3,)or not np.isfinite(np.r_[desired.ravel(),world.ravel(),origin]).all():raise ValueError('Finite corresponding controls and origin required')
    if w.shape!=(len(f),3,3)or not np.isfinite(w).all()or np.any(np.linalg.norm(w,axis=2)==0):raise ValueError('Nonzero finite corresponding corner fields required')
    if not np.isfinite(maximum_displacement_m)or maximum_displacement_m<=0 or not np.isfinite(wall_seconds)or wall_seconds<=0:raise ValueError('Finite positive source and wall limits required')
    if not isinstance(max_evaluations,int)or max_evaluations<0 or not isinstance(maximum_rounds,int)or maximum_rounds<1:raise ValueError('Finite evaluation and round budgets required')
    if candidate_mode not in ['local-neighbours','bounded-combined-lattice']:raise ValueError('Explicit local or bounded combined grid candidates required')
    if maximum_facet_rotation_degrees is not None and(not np.isfinite(maximum_facet_rotation_degrees)or not 0<=maximum_facet_rotation_degrees<90):raise ValueError('An explicit finite facet-rotation bound below 90 degrees is required')
    reference=world.copy();ref=reference[f];reference_normals=np.cross(ref[:,1]-ref[:,0],ref[:,2]-ref[:,0]);reference_lengths=np.linalg.norm(reference_normals,axis=1)
    if len(pairs)>2000 or len(pairs)and(pairs.min()<0 or pairs.max()>=len(f)or np.any(pairs[:,0]==pairs[:,1])):raise ValueError('At most 2000 explicit distinct source-face pairs required')
    relative=(world-origin).astype(np.float32)
    if not np.array_equal(relative.astype(float)+origin,world):raise ValueError('Stored controls must already use the declared Float32 grid')
    if np.linalg.norm(world-desired,axis=1).max(initial=0)>maximum_displacement_m:raise ValueError('Initial encoding exceeds the source bound')
    used=np.unique(f);occupancy={}
    if global_guard is not None and(not np.array_equal(global_guard.positions,world)or not np.array_equal(global_guard.faces,f)):raise ValueError('Global crossing guard must bind the actual initial coordinates and faces')
    for vertex in used:occupancy.setdefault(tuple(world[vertex]),set()).add(int(vertex))
    clusters={int(k):np.array(sorted(set(map(int,members))),int)for k,members in(source_clusters or {}).items()}
    for k,members in clusters.items():
        if k not in members or members.min()<0 or members.max()>=len(world)or np.linalg.norm(desired[members]-world[k],axis=1).max()>maximum_displacement_m:raise ValueError('Initial source alias cluster exceeds the bound')
    order=np.argsort(f.ravel(),kind='stable');counts=np.bincount(f.ravel(),minlength=len(world));starts=np.r_[0,np.cumsum(counts)];incident=order//3
    def field_ok(ids,vertex=None,value=None):
        p=world[f[ids]].copy()
        if vertex is not None:p[f[ids]==vertex]=value
        n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
        if vertex is not None and maximum_facet_rotation_degrees is not None:
            lengths=np.linalg.norm(n,axis=1)
            return bool(np.all(lengths>0)and np.all(np.einsum('fi,fi->f',n,reference_normals[ids])>=np.cos(np.radians(maximum_facet_rotation_degrees))*lengths*reference_lengths[ids]))
        return bool(np.all(np.einsum('fki,fi->fk',w[ids],n)>0))
    if not field_ok(np.arange(len(f))):raise ValueError('Repair crossing worklists only after the stored corner fields pass')
    def crossing(ids,vertex=None,value=None):
        tri=world[f[pairs[ids]]].copy()
        if vertex is not None:tri[f[pairs[ids]]==vertex]=value
        return np.isin(triangle_pairs(tri[:,0],tri[:,1])['relation'],['crossing','coplanar-overlap'])
    bad=crossing(np.arange(len(pairs)));initial=int(bad.sum());blocked=set();fixed=set(map(int,locked));moves=[];evaluations=0;reason=None
    quantum=float(np.max(np.abs(np.spacing(relative[used]))))
    lattice=None
    if candidate_mode=='bounded-combined-lattice':
        ratio=maximum_displacement_m/quantum
        if not np.isfinite(ratio)or ratio>7:raise ValueError('Combined search supports at most seven global grid quanta')
        radius=int(np.ceil(ratio))+1
        lattice=np.array(list(itertools.product(range(-radius,radius+1),repeat=3)),float)*quantum
    for _ in range(maximum_rounds):
        if not bad.any():break
        accepted=False
        for vertex in np.unique(f[pairs[bad]]):
            if time.monotonic()-started>=wall_seconds:reason='native wall deadline';break
            vertex=int(vertex)
            if vertex in fixed or vertex in blocked:continue
            ids=np.unique(incident[starts[vertex]:starts[vertex+1]])
            pids=np.flatnonzero(np.any(f[pairs]==vertex,axis=(1,2)))
            old_count=int(bad[pids].sum());point=relative[vertex];choices=[]
            for axis in range(3):choices.append([np.nextafter(point[axis],np.float32(-np.inf)),point[axis],np.nextafter(point[axis],np.float32(np.inf))])
            candidates=list(itertools.product(*choices))
            if lattice is not None:
                centre=(desired[vertex]-origin).astype(np.float32)
                q=(centre.astype(float)+lattice).astype(np.float32)
                eligible=np.linalg.norm(q.astype(float)+origin-desired[vertex],axis=1)<=maximum_displacement_m
                candidates.extend(map(tuple,q[eligible]))
            for axis in range(3):
                for sign in(-1,1):
                    q=point.copy();q[axis]=np.float32(float(q[axis])+sign*quantum);candidates.append(tuple(q))
            for pid in pids[bad[pids]]:
                for a,b in(pairs[pid],pairs[pid,::-1]):
                    if vertex not in f[a]:continue
                    other=world[f[b]];normal=np.cross(other[1]-other[0],other[2]-other[0]);normal/=np.linalg.norm(normal)
                    distances=(world[f[a]]-other[0])@normal
                    j=int(np.flatnonzero(f[a]==vertex)[0]);remaining=np.delete(distances,j)
                    if remaining[0]*remaining[1]<=0:continue
                    side=np.sign(remaining[0]);projected=world[vertex]-distances[j]*normal
                    q=(projected-origin).astype(np.float32);step=float(np.max(np.abs(np.spacing(q))))
                    for amount in(.01,.25,1.,2.):candidates.append(tuple((projected+side*step*amount*normal-origin).astype(np.float32)))
                    for axis in range(3):
                        for sign in(-1,1):
                            q2=q.copy();q2[axis]=np.nextafter(q2[axis],np.float32(sign*np.inf));candidates.append(tuple(q2))
            best=None;seen=set();members=clusters.get(vertex,np.array([vertex]))
            for proposal in candidates:
                if evaluations>=max_evaluations:reason='evaluation budget';break
                if time.monotonic()-started>=wall_seconds:reason='native wall deadline';break
                q=np.asarray(proposal,np.float32);key=tuple(q)
                if key in seen:continue
                seen.add(key);value=q.astype(float)+origin
                if np.array_equal(value,world[vertex])or occupancy.get(tuple(value),set())-{vertex}:continue
                movement=float(np.linalg.norm(desired[members]-value,axis=1).max())
                if movement>maximum_displacement_m:continue
                evaluations+=1
                if not field_ok(ids,vertex,value):continue
                after=crossing(pids,vertex,value)
                if np.any(after&~bad[pids])or int(after.sum())>=old_count:continue
                try:guarded=global_guard.check_move([vertex],value[None])if global_guard is not None else None
                except CrossingBudgetExceeded:reason='global pair-test budget';break
                if global_guard is not None and guarded is None:continue
                rank=(int(after.sum()),movement)
                if best is None or rank<best[0]:best=(rank,q,value,after,guarded)
            if best is None:
                blocked.add(vertex)
                if reason:break
                continue
            rank,q,value,after,guarded=best
            if global_guard is not None:global_guard.commit(guarded)
            old=tuple(world[vertex]);occupancy[old].remove(vertex)
            if not occupancy[old]:occupancy.pop(old)
            occupancy.setdefault(tuple(value),set()).add(vertex);relative[vertex]=q;world[vertex]=value;bad[pids]=after
            affected=set(map(int,np.unique(f[ids])))|set(map(int,np.unique(f[pairs[pids]])));blocked.difference_update(affected)
            moves.append(dict(vertex=vertex,remaining_declared_defect_pairs=int(bad.sum()),maximum_cluster_source_displacement_m=rank[1]));accepted=True;break
        if reason:break
        if not accepted:reason='no improving bounded grid move';break
    else:reason='round budget'
    passed=not bad.any()and reason not in ['native wall deadline','global pair-test budget']
    maximum=max([float(np.linalg.norm(world-desired,axis=1).max(initial=0))]+[float(np.linalg.norm(desired[members]-world[k],axis=1).max())for k,members in clusters.items()])
    output_words=w.astype(np.int16);transport=None
    if maximum_facet_rotation_degrees is not None:
        from geometry_tools.algorithms.surface_normals import transport_corner_words
        changed=np.flatnonzero(np.any(world[f]!=reference[f],axis=(1,2)))
        if len(changed):
            output_words[changed],transport=transport_corner_words(reference,world,f[changed],w[changed],maximum_facet_rotation_degrees=maximum_facet_rotation_degrees)
    return dict(positions=world,corner_words=output_words,evidence=dict(passed=passed,initial_declared_defect_pairs=initial,remaining_declared_defect_pairs=pairs[bad].tolist(),remaining_declared_defect_pair_count=int(bad.sum()),outward_fields_preserved=True,faces_words_uvs_unchanged=maximum_facet_rotation_degrees is None,normal_transport=transport,moves=moves,evaluations=evaluations,max_evaluations=max_evaluations,candidate_mode=candidate_mode,maximum_aggregate_source_displacement_m=maximum,declared_source_bound_m=maximum_displacement_m,native_wall_seconds=wall_seconds,budget_or_stop_reason=None if passed else reason,global_crossing_guard=global_guard.evidence()if global_guard is not None else None,global_crossings_topology_occupancy_and_anatomical_review_pending=True,elapsed_seconds=time.monotonic()-started))
