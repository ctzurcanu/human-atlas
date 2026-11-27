"""Bounded existing-edge contractions of a shared material boundary complex.

This only proposes local repairs. Each material must satisfy the manifold
edge-link condition, surviving corner fields cannot gain failures, and either
outward failures or an explicit verified crossing worklist must decrease.
Endpoint membership is preserved unless explicitly relaxed while testing
every affected material link. Full geometry audits remain required.
"""
from __future__ import annotations
from collections import Counter
import time
import numpy as np
from geometry_tools.algorithms.local_crossing_guard import CrossingBudgetExceeded


def collapse_short_edges(positions, surfaces, cells, corner_words, *,
        desired_positions=None, maximum_displacement_m=5e-8,
        maximum_steps=100, wall_seconds=15., candidate_face_rings=0,
        initial_clusters=None, defect_face_pairs=(),
        preserve_material_membership=True,locked=(),global_guard=None):
    started = time.monotonic()
    v = np.asarray(positions, float)
    desired = v if desired_positions is None else np.asarray(desired_positions, float)
    if v.ndim != 2 or v.shape[1] != 3 or desired.shape != v.shape or not np.isfinite(v).all() or not np.isfinite(desired).all():
        raise ValueError('Finite corresponding source and stored positions required')
    if not np.isfinite(maximum_displacement_m) or maximum_displacement_m <= 0 or not np.isfinite(wall_seconds) or wall_seconds <= 0:
        raise ValueError('Positive finite displacement and wall limits required')
    if not isinstance(maximum_steps, int) or maximum_steps < 0:
        raise ValueError('Finite nonnegative step budget required')
    if not isinstance(candidate_face_rings, int) or not 0 <= candidate_face_rings <= 2:
        raise ValueError('At most two explicitly bounded candidate face rings')
    if len(cells) > 4096:
        raise ValueError('At most 4096 explicitly owned materials per stage')
    if np.linalg.norm(v-desired, axis=1).max(initial=0) > maximum_displacement_m:
        raise ValueError('Input encoding already exceeds the source bound')
    names = list(surfaces)
    counts = [len(surfaces[s]) for s in names]
    offsets = np.r_[0, np.cumsum(counts)]
    f = np.concatenate([np.asarray(surfaces[s], np.int64) for s in names]).copy()
    fixed=set(map(int,locked))
    if global_guard is not None and(not np.array_equal(global_guard.positions,v)or not np.array_equal(global_guard.faces,f)):raise ValueError('Global crossing guard must bind the actual source coordinates and faces')
    w = np.concatenate([np.asarray(corner_words[s], float) for s in names])
    if f.shape != (len(f), 3) or w.shape != (len(f), 3, 3) or not np.isfinite(w).all() or f.min() < 0 or f.max() >= len(v):
        raise ValueError('Corresponding triangular faces and corner fields required')
    masks = {s: 0 for s in names}
    for i, uses in enumerate(cells.values()):
        for use in uses:
            masks[use['surface']] |= 1 << i
    if any(not mask for mask in masks.values()):
        raise ValueError('Every physical chart needs an explicit material owner')
    owner_dtype=np.uint64 if len(cells)<=63 else object
    owner = np.concatenate([np.full(n, masks[s], owner_dtype) for s, n in zip(names, counts)])
    def failures(triangles, words):
        p = v[triangles]
        n = np.cross(p[:, 1]-p[:, 0], p[:, 2]-p[:, 0])
        return np.any(np.einsum('fki,fi->fk', words, n) <= 0, axis=1)
    bad = failures(f, w)
    initial = int(bad.sum())
    active = np.ones(len(f), bool)
    known_pairs = np.asarray(defect_face_pairs, dtype=np.int64).reshape(-1, 2)
    if len(known_pairs) > 2000:
        raise ValueError('Split more than 2000 diagnostic face pairs into stages')
    if len(known_pairs) and (known_pairs.min() < 0 or known_pairs.max() >= len(f) or np.any(known_pairs[:, 0] == known_pairs[:, 1])):
        raise ValueError('Explicit distinct source-face pairs required')
    from geometry_tools.distance import triangle_pairs
    def known_failures(faces, live):
        result = np.zeros(len(known_pairs), bool)
        ids = np.flatnonzero(np.all(live[known_pairs], axis=1))
        if len(ids):
            pairs = known_pairs[ids]
            relation = triangle_pairs(v[faces[pairs[:, 0]]], v[faces[pairs[:, 1]]])['relation']
            result[ids] = np.isin(relation, ['crossing', 'coplanar-overlap'])
        return result
    known_bad = known_failures(f, active)
    initial_known_bad = int(known_bad.sum())
    order = np.argsort(f.ravel(), kind='stable')
    degree = np.bincount(f.ravel(), minlength=len(v))
    starts = np.r_[0, np.cumsum(degree)]
    incident = order//3
    cache = {}
    def star(vertex):
        if vertex not in cache:
            cache[vertex] = set(map(int, incident[starts[vertex]:starts[vertex+1]]))
        return cache[vertex]
    def closed_vertex_link(triangles, vertex):
        links = [tuple(int(q) for q in tri if q != vertex) for tri in triangles]
        adjacency = {}
        for a, b in links:
            adjacency.setdefault(a, []).append(b)
            adjacency.setdefault(b, []).append(a)
        if not adjacency or any(len(q) != 2 for q in adjacency.values()):
            return False
        reached, pending = set(), [next(iter(adjacency))]
        while pending:
            q = pending.pop()
            if q not in reached:
                reached.add(q)
                pending.extend(adjacency[q])
        return len(reached) == len(adjacency)
    clusters = {}
    assigned = set()
    used = set(map(int, np.unique(f)))
    for kept, members in (initial_clusters or {}).items():
        kept, members = int(kept), set(map(int, members))
        if kept not in members or kept not in used or min(members) < 0 or max(members) >= len(v) or assigned & members or (members-{kept}) & used:
            raise ValueError('Initial source clusters must be disjoint and have exactly one live control')
        if np.linalg.norm(desired[list(members)]-v[kept], axis=1).max() > maximum_displacement_m:
            raise ValueError('Initial source cluster exceeds the displacement bound')
        clusters[kept] = members
        assigned.update(members)
    moves = []
    rejected = Counter()
    stop = None
    for _ in range(maximum_steps):
        failed = np.flatnonzero(bad & active)
        if not len(failed) and not known_bad.any():
            break
        selection = set(map(int, failed)) | set(map(int, known_pairs[known_bad].ravel()))
        for _ in range(candidate_face_rings):
            selection.update(index for vertex in np.unique(f[sorted(selection)]) for index in star(int(vertex)))
        selected = f[sorted(selection)]
        edges = np.unique(np.sort(np.concatenate([selected[:, [0, 1]], selected[:, [1, 2]], selected[:, [2, 0]]]), axis=1), axis=0)
        distance = np.linalg.norm(v[edges[:, 0]]-v[edges[:, 1]], axis=1)
        nearby = distance <= 2*maximum_displacement_m
        edges, distance = edges[nearby], distance[nearby]
        candidates = edges[np.argsort(distance, kind='stable')]
        accepted = False
        for a, b in candidates:
            if time.monotonic()-started >= wall_seconds:
                stop = 'native wall deadline'
                break
            a, b = int(a), int(b)
            sa, sb = star(a), star(b)
            ia, ib = np.array(sorted(sa), int), np.array(sorted(sb), int)
            if not len(ia) or not len(ib):
                continue
            ma = int(np.bitwise_or.reduce(owner[ia]))
            mb = int(np.bitwise_or.reduce(owner[ib]))
            if preserve_material_membership and ma != mb:
                rejected['different material endpoint membership'] += 1
                continue
            patch = np.array(sorted(sa | sb), int)
            before = f[patch]
            # Test each individual closed material, including reciprocal charts.
            link_ok = True
            active_mask=ma|mb
            while active_mask:
                bit_value=active_mask & -active_mask
                active_mask-=bit_value
                mask=np.uint64(bit_value)if owner_dtype is np.uint64 else bit_value
                local = before[(owner[patch] & mask) != 0]
                at_a = local[np.any(local == a, axis=1)]
                at_b = local[np.any(local == b, axis=1)]
                if not len(at_a) or not len(at_b):
                    existing, vertex = (at_a, a) if len(at_a) else (at_b, b)
                    if not closed_vertex_link(existing, vertex):link_ok = False
                    if not link_ok:break
                    continue
                shared = local[np.any(local == a, axis=1) & np.any(local == b, axis=1)]
                opposite = set(map(int, shared.ravel()))-{a, b}
                common = (set(map(int, at_a.ravel()))-{a}) & (set(map(int, at_b.ravel()))-{b})
                if len(shared) != 2 or len(opposite) != 2 or common != opposite or not closed_vertex_link(at_a, a) or not closed_vertex_link(at_b, b):
                    link_ok = False
                    break
            if not link_ok:
                rejected['material manifold link'] += 1
                continue
            for keep, remove in [(a, b), (b, a)]:
                if remove in fixed:
                    rejected['locked source control']+=1
                    continue
                members = clusters.get(keep, {keep}) | clusters.get(remove, {remove})
                bound = float(np.linalg.norm(desired[list(members)]-v[keep], axis=1).max())
                if bound > maximum_displacement_m:
                    rejected['aggregate source bound'] += 1
                    continue
                after = before.copy()
                after[after == remove] = keep
                retain = np.all(np.diff(np.sort(after, axis=1), axis=1) != 0, axis=1)
                remaining = after[retain]
                if len(np.unique(np.sort(remaining, axis=1), axis=0)) != len(remaining):
                    rejected['duplicate physical face'] += 1
                    continue
                after_bad = failures(remaining, w[patch[retain]])
                before_count, after_count = int(bad[patch].sum()), int(after_bad.sum())
                if np.any(after_bad & ~bad[patch[retain]]) or after_count > before_count:
                    rejected['outward field or no strict improvement'] += 1
                    continue
                proposed_known = known_bad.copy()
                if len(known_pairs):
                    # Replay only the declared triangles, rather than copying
                    # the entire atlas surface for each numerical proposal.
                    proposed_known=np.zeros(len(known_pairs),bool)
                    live=np.all(active[known_pairs],axis=1)&~np.isin(known_pairs,patch[~retain]).any(1)
                    selected=np.flatnonzero(live);local=f[known_pairs[selected]].copy();local[local==remove]=keep
                    if len(selected):proposed_known[selected]=np.isin(triangle_pairs(v[local[:,0]],v[local[:,1]])['relation'],['crossing','coplanar-overlap'])
                if np.any(proposed_known & ~known_bad) or int(proposed_known.sum()) > int(known_bad.sum()) or (after_count == before_count and int(proposed_known.sum()) == int(known_bad.sum())):
                    rejected['no outward or declared crossing improvement'] += 1
                    continue
                guarded=None
                if global_guard is not None:
                    try:guarded=global_guard.check_faces(patch,after,retain)
                    except CrossingBudgetExceeded:stop='global pair-test budget';break
                    if guarded is None:
                        rejected['new global crossing']+=1
                        continue
                    global_guard.commit(guarded)
                # Only existing endpoints and unchanged retained attributes are used.
                touched = set(map(int, before.ravel())) | set(map(int, remaining.ravel()))
                for vertex in touched:
                    star(vertex).difference_update(map(int, patch))
                removed_faces = patch[~retain]
                active[removed_faces] = False
                bad[removed_faces] = False
                f[patch[retain]] = remaining
                bad[patch[retain]] = after_bad
                known_bad = proposed_known
                for index, tri in zip(patch[retain], remaining):
                    for vertex in tri:
                        star(int(vertex)).add(int(index))
                clusters[keep] = members
                clusters.pop(remove, None)
                moves.append(dict(keep=keep, remove=remove,
                    deleted_physical_faces=removed_faces.tolist(),
                    maximum_cluster_source_displacement_m=bound,
                    remaining_bad_faces=int((bad & active).sum())))
                accepted = True
                break
            if accepted:
                break
            if stop:break
        if stop:
            break
        if not accepted:
            stop = 'no improving bounded manifold edge contraction'
            break
    else:
        stop = 'step budget'
    out, source_faces = {}, {}
    for i, sid in enumerate(names):
        ids = np.arange(offsets[i], offsets[i+1])
        ids = ids[active[ids]]
        out[sid] = f[ids]
        source_faces[sid] = ids-offsets[i]
    return dict(surfaces=out, source_faces=source_faces,
        evidence=dict(passed=not np.any(bad & active) and not known_bad.any() and stop != 'native wall deadline',
            initial_bad_faces=initial, remaining_bad_faces=int((bad & active).sum()),
            remaining_global_source_faces=np.flatnonzero(bad & active).tolist(),
            initial_declared_defect_pairs=initial_known_bad,
            remaining_declared_defect_pairs=known_pairs[known_bad].tolist(),
            material_endpoint_membership_preserved=preserve_material_membership,
            contractions=moves, rejected=dict(rejected),
            source_clusters=[dict(kept=k, members=sorted(members)) for k, members in sorted(clusters.items())],
            removed_physical_faces=int((~active).sum()),
            surviving_corner_attributes_literal=True, original_position_array_literal=True,
            source_cluster_bound_m=maximum_displacement_m,
            maximum_accepted_cluster_displacement_m=max((q['maximum_cluster_source_displacement_m'] for q in moves), default=0.),
            maximum_steps=maximum_steps, native_wall_seconds=wall_seconds,
            candidate_face_rings=candidate_face_rings,
            global_crossing_guard=global_guard.evidence()if global_guard is not None else None,
            budget_or_stop_reason=stop, elapsed_seconds=time.monotonic()-started,
            full_geometry_occupancy_and_anatomical_review_pending=True))
