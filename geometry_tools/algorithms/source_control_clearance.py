"""Bounded local clearance with literal sections and movable original cores.

Every source physical point has explicit ownership: a rigidly translated
recorded section, or its own source core control. Original mesh edges supply
regularization and projected edge-length guards. No geometry is invented.
"""
import time

import numpy as np
from scipy.spatial import cKDTree

from geometry_tools.open3d_backend import library
from .section_course_socp import solve_weighted_section_socp
from .profile_graph_transport import finite_diagnostic
from .feature_normals import indexed_feature_normals,metric_oriented_rows,inward_sheet_clearance


def source_control_basis(vertices, faces, profiles, physical_ids, *, micro_coherence_m=0.,core_patch_diameter_m=0.):
    v = np.asarray(vertices, float)
    _, representative, inv = np.unique(physical_ids, return_index=True, return_inverse=True)
    owner = np.full(len(representative), -1, int)
    centres = []
    for k, raw in enumerate(profiles):
        physical = inv[raw]
        if len(np.unique(physical)) != len(raw) or np.any(owner[physical] >= 0):
            raise ValueError('Recorded sections need disjoint actual source controls')
        owner[physical] = k
        centres.append(v[raw].mean(0))
    for raw in representative[owner < 0]:
        owner[inv[raw]] = len(centres)
        centres.append(v[raw])
    if micro_coherence_m:
        # A common displacement preserves the actual micron-scale source
        # offsets. This does not merge source points, faces or identities.
        parent = np.arange(len(centres))
        def root(k):
            while parent[k] != k:
                parent[k] = parent[parent[k]]
                k = parent[k]
            return k
        for a, b in cKDTree(v[representative]).query_pairs(micro_coherence_m):
            ra, rb = root(owner[a]), root(owner[b])
            if ra != rb:
                parent[rb] = ra
        labels = np.array([root(k) for k in range(len(centres))])
        unique, compact = np.unique(labels, return_inverse=True)
        centres = [np.mean(np.asarray(centres)[compact == k], axis=0) for k in range(len(unique))]
        owner = compact[owner]
    source_edges = np.unique(np.sort(faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1), axis=0)
    if core_patch_diameter_m:
        if not 0<core_patch_diameter_m<=.003:
            raise ValueError('Explicit source-core patch diameter within (0,3] mm required')
        # Group displacement variables only. Every original coordinate and
        # face survives; a patch is joined exclusively through source edges.
        # Actual section groups are protected from this agglomeration.
        physical=v[representative]
        pinned=set(int(owner[inv[raw]])for section in profiles for raw in section)
        parent=np.arange(len(centres));members=[np.flatnonzero(owner==k)for k in range(len(centres))]
        def root(k):
            while parent[k]!=k:parent[k]=parent[parent[k]];k=int(parent[k])
            return k
        edges=np.unique(np.sort(inv[source_edges],axis=1),axis=0)
        length=np.linalg.norm(physical[edges[:,1]]-physical[edges[:,0]],axis=1)
        for edge in edges[np.argsort(length,kind='stable')]:
            a,b=map(root,owner[edge])
            if a==b or a in pinned or b in pinned:continue
            ids=np.r_[members[a],members[b]];p=physical[ids]
            if np.max(np.linalg.norm(p[:,None]-p[None,:],axis=2))>core_patch_diameter_m:continue
            parent[b]=a;members[a]=ids
        labels=np.array([root(k)for k in range(len(centres))]);_,compact=np.unique(labels,return_inverse=True)
        owner=compact[owner]
        centres=[physical[owner==k].mean(0)for k in range(int(owner.max())+1)]
    W = np.eye(len(centres))[owner[inv]]
    links = np.unique(np.sort(owner[inv[source_edges]], axis=1), axis=0)
    links = links[links[:, 0] != links[:, 1]]
    return W, np.asarray(centres), links, source_edges


def clear_source_controls(vertices, faces, profiles, physical_ids, bone, skin, *,
                          maximum_motion_m=.003, clearance_m=.0002,
                          maximum_iterations=8, wall_seconds=2., guard_face_pairs=(),
                          protected_raw=(),preserve_dihedral_sign=True,
                          recipient_surfaces=(), rigid_translation=False,initial_positions=None,solver_seconds=.4,
                          harmonic_cores=False,core_patch_diameter_m=0.,recipient_clearance_m=None,
                          elastic_recipient_proposal=False,recipient_plane_guards=()):
    if not 0 < maximum_motion_m <= .01 or not 0 < wall_seconds <= 3:
        raise ValueError('Explicit motion and finite native wall bounds required')
    if not 0<solver_seconds<=1 or solver_seconds>wall_seconds:
        raise ValueError('Explicit optimizer deadline at most one second and within native wall limit required')
    recipient_clearance_m=clearance_m if recipient_clearance_m is None else float(recipient_clearance_m)
    if not 0<recipient_clearance_m<=.001:raise ValueError('Explicit positive recipient clearance within 1 mm required')
    start = time.monotonic()
    v, f = np.asarray(vertices, float), np.asarray(faces, int)
    if harmonic_cores and core_patch_diameter_m:raise ValueError('Choose harmonic cores or source-edge patches, not both')
    W, centres, _links, edges = source_control_basis(v, f, profiles, physical_ids, micro_coherence_m=.00002,
                                                   core_patch_diameter_m=core_patch_diameter_m)
    if harmonic_cores:
        from .profile_core_transport import profile_translation_weights
        W=profile_translation_weights(v,f,profiles,physical_ids=physical_ids)
        centres=np.array([v[s].mean(0)for s in profiles])
    if rigid_translation:
        W=np.ones((len(v),1));centres=v.mean(0,keepdims=True)
    count = len(centres)
    initial_delta=np.zeros((count,3))
    if initial_positions is not None:
        initial=np.asarray(initial_positions,float)
        if initial.shape!=v.shape or not np.isfinite(initial).all():
            raise ValueError('Finite corresponding source warm-start positions required')
        initial_delta=np.linalg.lstsq(W,initial-v,rcond=None)[0]
        if np.max(abs(v+W@initial_delta-initial))>1e-10:
            raise ValueError('Warm start must retain every rigid original source control group')
        if np.linalg.norm(initial_delta,axis=1).max()>maximum_motion_m+1e-10:
            raise ValueError('Warm start exceeds the declared source motion bound')
    _, representatives, identities = np.unique(physical_ids, return_index=True, return_inverse=True)
    physical_edges = np.unique(np.sort(identities[edges], axis=1), axis=0)
    physical_edges = physical_edges[physical_edges[:, 0] != physical_edges[:, 1]]
    pa, pb = physical_edges.T
    raw_a, raw_b = representatives[pa], representatives[pb]
    edge_controls = W[raw_a]-W[raw_b]
    rest_edges = v[raw_a]-v[raw_b]
    H = edge_controls.T@edge_controls+.01*np.eye(count)
    protected_raw=np.asarray(protected_raw,int)
    if protected_raw.size and ((protected_raw<0).any()or(protected_raw>=len(v)).any()):
        raise ValueError('Protected source controls must be actual raw vertex IDs')
    fixed=np.flatnonzero(np.any(W[protected_raw]!=0,axis=0))if protected_raw.size else np.array([],int)
    if len(fixed)and np.max(abs(initial_delta[fixed]))>1e-12:
        raise ValueError('Warm start moved a protected original source control')
    free=np.setdiff1d(np.arange(count),fixed)
    K=np.zeros_like(H)
    K[np.ix_(free,free)]=np.linalg.inv(H[np.ix_(free,free)])
    K = (K+K.T)/2
    points = [v, v[f].mean(1)]
    weights = [W, W[f].mean(1)]
    for a, b in [(0, 1), (1, 2), (2, 0)]:
        points.append((v[f[:, a]]+v[f[:, b]])/2)
        weights.append((W[f[:, a]]+W[f[:, b]])/2)
    points, weights = np.concatenate(points), np.concatenate(weights)
    plane_guards=[]
    for guard in recipient_plane_guards:
        raw=np.asarray(guard['sourceRaw'],int);normal=np.asarray(guard['normal'],float);point=np.asarray(guard['point_m'],float)
        margin=float(guard.get('clearance_m',recipient_clearance_m))
        if (not len(raw) or raw.min()<0 or raw.max()>=len(v) or normal.shape!=(3,) or point.shape!=(3,) or
            not np.isfinite(normal).all() or not np.isfinite(point).all() or
            abs(np.linalg.norm(normal)-1)>1e-12 or not 0<margin<=.001):
            raise ValueError('A source-face guard needs actual raw controls and a finite unit plane')
        controls=np.any(W[raw]!=0,axis=0)
        support=np.any(abs(weights[:,controls])>1e-12,axis=1)
        plane_guards.append((guard['partId'],raw,point,normal,margin,support))
    # Open recipient sheets have no volume. Keep the initial control-centre
    # side as a numerical route choice, never assert closed-solid containment.
    # Full triangle and anatomical review remain mandatory after this proposal.
    recipient=[]
    for surface in recipient_surfaces:
        index,closed=surface[:2]
        policy=surface[2] if len(surface)>2 else 'initial-control-side'
        label=surface[3] if len(surface)>3 else str(len(recipient))
        if policy not in ['initial-control-side','bone-exterior','oriented-exterior','reference-positive','unsigned-sheet','coherent-section-exterior','directed-exterior']:
            raise ValueError('Explicit recipient-sheet side policy required')
        if policy in ['coherent-section-exterior','directed-exterior']and not closed:
            raise ValueError('A coherent section exterior needs an explicitly closed oriented recipient')
        escape=None if len(surface)<5 else surface[4]
        if policy=='directed-exterior':
            escape=np.asarray(escape,float)
            if escape.shape!=(3,)or not np.isfinite(escape).all()or abs(np.linalg.norm(escape)-1)>1e-12:
                raise ValueError('An explicit finite unit escape direction is required')
        side_centres=centres+initial_delta
        hit,_=index.closest_points(side_centres)
        normal,_=indexed_feature_normals(index,hit)
        side=np.where(np.einsum('ij,ij->i',side_centres-hit['points'],normal)>=0,1.,-1.)
        recipient.append((index,bool(closed),side[np.argmax(weights,axis=1)],policy,label,escape))

    def recipient_constraints(q):
        normals=[];gaps=[];rows=[];witnesses=[]
        for surface_id,(index,closed,side,policy,label,escape) in enumerate(recipient):
            hit,_=index.closest_points(q)
            if closed:
                tensor=library().core.Tensor(((q-index.origin)/index.scale).astype(np.float32))
                inside=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
                normal=q-hit['points'];normal[inside]*=-1
                lengths=np.linalg.norm(normal,axis=1)
                if np.any(lengths<1e-12):raise ValueError('Ambiguous recipient surface witness')
                normal/=lengths[:,None]
                gap=hit['distance_m']*np.where(inside,-1.,1.)
                near_or_inside=inside|(hit['distance_m']<=recipient_clearance_m+.00001)
                # If a rigid source patch straddles the recipient, its clear
                # vertices must take the same escape as its penetrating ones.
                # Otherwise their old-side halfspaces forbid the translation.
                # Restrict propagation to actual shared source control support.
                contact_controls=np.any(abs(weights[near_or_inside])>1e-12,axis=0)
                contact_rows=np.any(abs(weights[:,contact_controls])>1e-12,axis=1)
                if policy=='coherent-section-exterior'and np.any(near_or_inside&(np.count_nonzero(abs(weights)>1e-12,axis=1)==1)):
                    # A narrow rod can cut opposite vertices of the same rigid
                    # section. Independent closest-point gradients then demand
                    # contradictory translations. Select one exterior tangent
                    # plane at that section's actual course centre instead.
                    # This is a conservative local proposal; full surface and
                    # containment replays remain mandatory after movement.
                    # Harmonic branch-core weights can span distant sections;
                    # their weighted centroid is not a physical vessel axis.
                    # Use a shared plane only for an actual rigid section.
                    coherent=near_or_inside&(np.count_nonzero(abs(weights)>1e-12,axis=1)==1)
                    course=weights[coherent]@centres+(q-points)[coherent]
                    ch,_=index.closest_points(course)
                    tensor=library().core.Tensor(((course-index.origin)/index.scale).astype(np.float32))
                    ci=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
                    direction=course-ch['points'];direction[ci]*=-1
                    length=np.linalg.norm(direction,axis=1)
                    fallback,_=indexed_feature_normals(index,ch)
                    direction=np.divide(direction,length[:,None],out=fallback.copy(),where=length[:,None]>1e-12)
                    normal[coherent]=direction
                    gap[coherent]=np.einsum('ij,ij->i',q[coherent]-ch['points'],direction)
                elif policy=='directed-exterior'and contact_rows.any():
                    # Trace from an independently checked interior probe to
                    # the actual closed recipient's first exit. A common
                    # declared direction avoids opposite escapes along one
                    # source course. This local plane is not a packing proof.
                    if harmonic_cores:raise ValueError('Directed patch escape requires local source controls, not global harmonic supports')
                    chosen=np.flatnonzero(contact_rows)
                    feature,_=indexed_feature_normals(index,hit)
                    feature=feature[chosen];surface_point=hit['points'][chosen]
                    probe=surface_point-feature*.000001
                    tensor=library().core.Tensor(((probe-index.origin)/index.scale).astype(np.float32))
                    valid=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=7).numpy()>.5
                    if not valid.all():
                        probe[~valid]=surface_point[~valid]+feature[~valid]*.000001
                        tensor=library().core.Tensor(((probe-index.origin)/index.scale).astype(np.float32))
                        if not (index.scene.compute_occupancy(tensor,nthreads=2,nsamples=7).numpy()>.5).all():
                            raise ValueError('Directed recipient exit has an unresolved interior probe')
                    direction=np.broadcast_to(escape,(len(chosen),3)).copy()
                    rays=np.c_[(probe-index.origin)/index.scale,direction]
                    travel=index.scene.cast_rays(library().core.Tensor(rays.astype(np.float32)),nthreads=2)['t_hit'].numpy()*index.scale
                    if not np.isfinite(travel).all()or np.any(travel<0)or np.any(travel>2*maximum_motion_m+.002):
                        raise ValueError('Directed recipient exit exceeds its explicit finite ray bound')
                    exit_point=probe+direction*travel[:,None]
                    normal[chosen]=direction
                    gap[chosen]=np.einsum('ij,ij->i',q[chosen]-exit_point,direction)
            else:
                normal,boundary=indexed_feature_normals(index,hit)
                if policy=='bone-exterior':
                    bone_hit,_=bone.closest_points(hit['points'])
                    radial=hit['points']-bone_hit['points']
                    polarity=np.einsum('ij,ij->i',normal,radial)
                    if np.any(abs(polarity)<1e-10):
                        raise ValueError('Bone surface does not determine the declared sheet exterior')
                    normal*=np.where(polarity>=0,1.,-1.)[:,None]
                elif policy=='initial-control-side':
                    normal*=side[:,None]
                unsigned=np.ones(len(q),bool)if policy=='unsigned-sheet'else boundary
                normal,gap=metric_oriented_rows(q,hit,normal,unsigned_boundary=unsigned)
            for guarded_label,raw,point,plane_normal,margin,support in plane_guards:
                if guarded_label==label:
                    normal[support]=plane_normal
                    gap[support]=(q[support]-point)@plane_normal
            # A finite sheet must not create an infinite obstacle plane.
            active=hit['distance_m']<=maximum_motion_m+recipient_clearance_m
            if closed:active|=inside
            normals.append(normal[active]);gaps.append(gap[active]);rows.append(weights[active])
            witnesses.append(np.column_stack([np.full(active.sum(),surface_id),np.flatnonzero(active),hit['face_ids'][active]]))
        return (np.concatenate(normals),np.concatenate(gaps),np.concatenate(rows),np.concatenate(witnesses).astype(int)) if recipient else (np.empty((0,3)),np.empty(0),np.empty((0,count)),np.empty((0,3),int))
    d = v[edges[:, 1]]-v[edges[:, 0]]
    length = np.linalg.norm(d, axis=1)
    ew = W[edges[:, 1]]-W[edges[:, 0]]
    live = (length > 1e-12)&np.any(ew != 0, axis=1)
    ew, en, eb = ew[live], d[live]/length[live, None], -.35*length[live]
    # Every vertex of each actual crossing source triangle must lie beyond
    # the declared common plane. This is stricter than sampling only its
    # centroid and edge midpoints. Final original-surface replay is separate.
    for label,raw,point,normal,margin,support in plane_guards:
        ew=np.r_[ew,W[raw]];en=np.r_[en,np.broadcast_to(normal,(len(raw),3))]
        eb=np.r_[eb,margin-(v[raw]-point)@normal]
    # Preserve pre-existing tiny clearances between distinct source controls.
    # This is a nonpenetration inequality, never a tolerance weld or an
    # anatomical attachment. Independent source coordinates remain independent.
    near_pairs = cKDTree(v[representatives]).query_pairs(.00002, output_type='ndarray')
    if len(near_pairs):
        near_raw = representatives[near_pairs]
        md = v[near_raw[:, 1]]-v[near_raw[:, 0]]
        ml = np.linalg.norm(md, axis=1)
        mw = W[near_raw[:, 1]]-W[near_raw[:, 0]]
        keep = (ml > 1e-10)&np.any(mw != 0, axis=1)
        mw, mn, mb = mw[keep], md[keep]/ml[keep, None], -.2*ml[keep]
    else:
        mw, mn, mb = np.empty((0, count)), np.empty((0, 3)), np.empty(0)
    delta = initial_delta.copy()
    runs = []

    def shape_target(displacement):
        current_edges = rest_edges+edge_controls@displacement
        covariance = np.zeros((len(representatives), 3, 3))
        outer = np.einsum('ni,nj->nij', current_edges, rest_edges)
        np.add.at(covariance, pa, outer)
        np.add.at(covariance, pb, outer)
        u, _, vt = np.linalg.svd(covariance)
        sign = np.ones((len(u), 3))
        sign[:, -1] = np.linalg.det(u@vt)
        rotations = (u*sign[:, None, :])@vt
        desired = np.einsum('nij,nj->ni', (rotations[pa]+rotations[pb])/2, rest_edges)
        return K@(edge_controls.T@(desired-rest_edges))
    initial_cross = np.cross(v[f[:, 1]]-v[f[:, 0]], v[f[:, 2]]-v[f[:, 0]])
    initial_area = np.linalg.norm(initial_cross, axis=1)
    if np.any(initial_area <= 1e-18):
        raise ValueError('A degenerate source triangle cannot define a local orientation guard')
    face_normal = initial_cross/initial_area[:, None]

    def face_constraints(displacement):
        t = (v+W@displacement)[f]
        normal = face_normal
        gradients = np.stack([np.cross(normal, t[:, 2]-t[:, 1]),
                              np.cross(t[:, 2]-t[:, 0], normal),
                              np.cross(normal, t[:, 1]-t[:, 0])], axis=1)
        # Normalize to a length residual for the solver's metre tolerance.
        support = np.linalg.norm(gradients, axis=2).sum(1)
        A = np.einsum('fki,fkn->fni', gradients, W[f])/support[:, None, None]
        signed = np.einsum('ij,ij->i', normal, np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]))
        offset = (.22*initial_area-signed)/support+np.einsum('fni,ni->f', A, displacement)
        return A, offset, signed/initial_area

    # Four actual neighbouring source controls constrain local fold-through.
    # These are numerical orientation guards, not constructed tissue cells.
    _, rep, physical = np.unique(physical_ids, return_index=True, return_inverse=True)
    adjacent = {}
    for tri in physical[f]:
        for a, b, other in [(tri[0], tri[1], tri[2]), (tri[1], tri[2], tri[0]), (tri[2], tri[0], tri[1])]:
            adjacent.setdefault(tuple(sorted((a, b))), []).append(other)
    cells = np.array([[rep[a], rep[b], rep[other[0]], rep[other[1]]]
                      for (a, b), other in adjacent.items() if len(other) == 2 and preserve_dihedral_sign], int).reshape(-1,4)
    extra_cells = []
    for first, second in guard_face_pairs:
        a, b = physical[f[first]], physical[f[second]]
        candidates = [list(a)+[x] for x in b if x not in a]
        candidates += [list(b)+[x] for x in a if x not in b]
        for i in range(3):
            for j in range(3):
                candidates.append([a[i], a[(i+1)%3], b[j], b[(j+1)%3]])
        extra_cells.extend([rep[np.array(sorted(x), int)] for x in candidates if len(set(x)) == 4])
    if extra_cells:
        cells = np.unique(np.sort(np.r_[cells, extra_cells], axis=1), axis=0)
    source_t = v[cells]
    source_e = source_t[:, 1:]-source_t[:, :1]
    volume = np.einsum('ij,ij->i', source_e[:, 0], np.cross(source_e[:, 1], source_e[:, 2]))
    denominator = np.prod(np.linalg.norm(source_e, axis=2), axis=1)
    eligible = np.abs(volume) > .02*denominator
    cells, volume = cells[eligible], volume[eligible]
    signs = np.sign(volume)
    volume = np.abs(volume)

    def cell_constraints(displacement):
        t = (v+W@displacement)[cells]
        e = t[:, 1:]-t[:, :1]
        gb, gc, gd = np.cross(e[:, 1], e[:, 2]), np.cross(e[:, 2], e[:, 0]), np.cross(e[:, 0], e[:, 1])
        gradients = np.stack([-gb-gc-gd, gb, gc, gd], axis=1)*signs[:, None, None]
        support = np.linalg.norm(gradients, axis=2).sum(1)
        if np.any(support < 1e-20):
            raise ValueError('Collapsed source orientation guard')
        A = np.einsum('fki,fkn->fni', gradients, W[cells])/support[:, None, None]
        signed = np.einsum('ij,ij->i', e[:, 0], np.cross(e[:, 1], e[:, 2]))*signs
        offset = (.12*volume-signed)/support+np.einsum('fni,ni->f', A, displacement)
        return A, offset, signed/volume

    def constraints(q):
        near, _ = bone.closest_points(q)
        tensor = library().core.Tensor(((q-bone.origin)/bone.scale).astype(np.float32))
        inside = bone.scene.compute_occupancy(tensor, nthreads=2, nsamples=3).numpy() > .5
        toward = q-near['points']
        toward[inside] *= -1
        lengths = np.linalg.norm(toward, axis=1)
        if np.any(lengths < 1e-12):
            raise ValueError('Ambiguous bone witness at an exact surface point')
        bn = toward/lengths[:, None]
        bone_gap = near['distance_m']*np.where(inside, -1., 1.)
        sn, _ = skin.closest_points(q)
        outward,boundary=indexed_feature_normals(skin,sn)
        inward,skin_gap,_=inward_sheet_clearance(q,sn,outward,boundary)
        return bn, bone_gap, inward, skin_gap

    complete = False
    for iteration in range(maximum_iterations):
        if time.monotonic()-start > wall_seconds:
            return v.copy(), dict(complete=False, status='NativeDeadline', runs=runs,
                                  elapsedSeconds=time.monotonic()-start)
        bn, bg, sn, sg = constraints(points+weights@delta)
        rn,rg,rw,recipient_witnesses=recipient_constraints(points+weights@delta)
        _, _, area_ratio = face_constraints(delta)
        _, _, cell_ratio = cell_constraints(delta)
        if iteration >= 3 and min(bg.min(), sg.min()) >= clearance_m-1e-8 and rg.min(initial=1.)>=recipient_clearance_m-1e-8 and area_ratio.min() >= .2-1e-6 and cell_ratio.min(initial=1) >= .1-1e-6:
            complete = True
            break
        # Rebuild local halfspaces at the current pose. Obsolete tangent planes
        # of nonconvex bones/skin are not permanent anatomical constraints.
        shifts = weights@delta
        target_gap = clearance_m+.00001
        rhs = np.r_[eb, mb, target_gap-bg+np.einsum('ij,ij->i', bn, shifts),
                    target_gap-sg+np.einsum('ij,ij->i', sn, shifts),
                    recipient_clearance_m+.00001-rg+np.einsum('ij,ij->i',rn,rw@delta)]
        area_rows, area_offsets, _ = face_constraints(delta)
        cell_rows, cell_offsets, _ = cell_constraints(delta)
        proposal, proof = solve_weighted_section_socp(
            K, np.r_[ew, mw, weights, weights,rw], np.r_[en, mn, bn, sn,rn], rhs,
            maximum_motion_m=maximum_motion_m, residual_m=1e-8,
            time_limit_seconds=min(solver_seconds, max(.01, wall_seconds-(time.monotonic()-start))),
            maximum_iterations=60, precision_matrix=H if not len(fixed)else None,
            vector_rows=np.r_[area_rows, cell_rows], vector_offsets=np.r_[area_offsets, cell_offsets],
            target_displacement=shape_target(delta),
            elastic_row_groups=None if not elastic_recipient_proposal else np.r_[
                np.full(len(ew)+len(mw)+2*len(points),-1,int),recipient_witnesses[:,0]])
        # Retain useful aggregate diagnostics; thousands of trivial per-row
        # lower bounds are reproducible from the source constraints.
        proof.pop('constraintRowNecessaryRadiusLowerBoundMeters', None)
        proof.pop('freeStations', None)
        base=len(ew)+len(mw)+2*len(points)
        for witness in proof.get('leadingHalfspaceDualRows',[]):
            k=witness['row']-base
            if 0<=k<len(recipient_witnesses):
                surface_id,sample_id,face_id=recipient_witnesses[k]
                witness.update(recipientSurface=recipient[surface_id][4],recipientFace=int(face_id),
                    sample=int(sample_id),point_m=(points+weights@delta)[sample_id].tolist(),
                    normal=rn[k].tolist(),signedGap_m=float(rg[k]),
                    sourceControls=np.flatnonzero(rw[k]).tolist(),sheetPolicy=recipient[surface_id][3])
        runs.append(finite_diagnostic(proof))
        if proposal is None:
            break
        delta = proposal
    q = v+W@delta
    if protected_raw.size and not np.array_equal(q[protected_raw],v[protected_raw]):
        raise ValueError('A protected original source mouth moved')
    _, bg, _, sg = constraints(points+weights@delta)
    _,rg,_,_=recipient_constraints(points+weights@delta)
    actual_recipient_gap=1.
    final_samples=points+weights@delta
    for index,closed,side,policy,label,escape in recipient:
        hit,_=index.closest_points(final_samples);gap=hit['distance_m']
        if closed:
            tensor=library().core.Tensor(((final_samples-index.origin)/index.scale).astype(np.float32))
            inside=index.scene.compute_occupancy(tensor,nthreads=2,nsamples=3).numpy()>.5
            gap=gap*np.where(inside,-1.,1.)
        actual_recipient_gap=min(actual_recipient_gap,float(gap.min()))
    plane_gap=min((float(np.min((q[raw]-point)@normal-margin))for label,raw,point,normal,margin,support in plane_guards),default=1.)
    worst = int(np.argmin(sg))
    witness = points[worst]+weights[worst]@delta
    skin_hit, _ = skin.closest_points(witness[None])
    all_skin_hits,_=skin.closest_points(points+weights@delta)
    _,skin_rim=indexed_feature_normals(skin,all_skin_hits)
    _, _, area_ratio = face_constraints(delta)
    _, _, cell_ratio = cell_constraints(delta)
    complete = bool(min(bg.min(), sg.min()) >= clearance_m-1e-8 and rg.min(initial=1.)>=recipient_clearance_m-1e-8 and actual_recipient_gap>=recipient_clearance_m-1e-8 and plane_gap>=-1e-8 and area_ratio.min() >= .2-1e-6 and cell_ratio.min(initial=1) >= .1-1e-6)
    elapsed=time.monotonic()-start
    timed_out=elapsed>wall_seconds
    complete=complete and not timed_out
    gauge = max((float(np.abs(np.linalg.norm(q[s][:, None]-q[s][None, :], axis=2)-
                                    np.linalg.norm(v[s][:, None]-v[s][None, :], axis=2)).max())
                 for s in profiles), default=0.)
    return q, dict(complete=complete, status='Proposed' if complete else 'ClearanceUnresolved',
                   runs=runs, controls=count, rigidSections=len(profiles),
                   minimumSampledBoneGap_m=float(bg.min()), minimumSampledSkinGap_m=float(sg.min()),
                   worstSkinSample=dict(index=worst, point_m=witness.tolist(),
                                        face=int(skin_hit['face_ids'][0]),
                                        nearestPoint_m=skin_hit['points'][0].tolist()),
                   maximumControlMotion_m=float(np.linalg.norm(delta, axis=1).max()),
                   sectionGaugeError_m=gauge, sourceFacesRetained=True,
                   minimumProjectedSourceFaceAreaRatio=float(area_ratio.min()),
                   sourceEdgeOrientationGuards=len(cells), minimumSourceEdgeOrientationRatio=float(cell_ratio.min(initial=1)),
                   sourceDihedralSignPreserved=preserve_dihedral_sign,
                   sourceMicroclearanceGuards=len(mb), sourceControlsWelded=False,
                   sourceMicroCoherenceRadius_m=.00002,
                   fixedSourceControls=len(fixed),protectedSourceRawControls=len(protected_raw),
                   timedOut=timed_out,nativeWallSeconds=wall_seconds,
                   optimizerWallSeconds=solver_seconds,
                   recipientSurfaceCount=len(recipient),minimumSampledRecipientGap_m=float(rg.min(initial=1.)),
                   recipientClearance_m=recipient_clearance_m,
                   recipientSheetSidePolicies=[r[3]for r in recipient],rigidTranslation=rigid_translation,
                   warmStartUsed=initial_positions is not None,
                   harmonicSourceCoreTransport=harmonic_cores,
                   sourceCorePatchDiameterBound_m=core_patch_diameter_m,
                   coherentPlanesLimitedToContactControlSupport=True,
                   closestFeatureNormals='Angle-weighted vertex, incident-face edge sum, original face interior',
                   finiteOpenBoundaryFeaturesUnsigned=True,
                   skinOpenRimWitnesses=int(skin_rim.sum()),
                   actualTriangleAuditRequired=True, openSkinContainmentNotEstablished=True,
                   elasticRecipientProposal=elastic_recipient_proposal,
                   finalClearanceChecksUnrelaxed=True,
                   sourceFacePlaneGuards=len(plane_guards),minimumSourceFacePlaneSurplus_m=plane_gap,
                   minimumActualRecipientGap_m=actual_recipient_gap,
                   elapsedSeconds=elapsed)
