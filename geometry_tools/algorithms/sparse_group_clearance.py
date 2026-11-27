"""Sparse shared-control clearance; no dense inverse or independent port motion.

Rigid groups translate together. Open-skin sidedness proposes constraints only;
actual emitted triangles must be checked separately. All lengths are metres.
"""
import time

import clarabel
import numpy as np
from scipy import sparse

from geometry_tools.open3d_backend import library
from .feature_normals import indexed_feature_normals, inward_sheet_clearance


def sparse_group_clearance(vertices, faces, groups, bone, skin, *,
                           fixed_above=1.08, maximum_motion_m=.005,
                           clearance_m=.0002, wall_seconds=10., maximum_iterations=4,
                           coherent_bone_sections=False,active_raw=None):
    started = time.monotonic()
    v, f = np.asarray(vertices, float), np.asarray(faces, int)
    if not 0 < wall_seconds <= 12 or not 0 < maximum_motion_m <= .020 or not 1 <= maximum_iterations <= 6:
        raise ValueError('Explicit bounded sparse clearance budgets required')
    owner = np.full(len(v), -1, int)
    for i, group in enumerate(groups):
        if np.any(owner[group] >= 0):
            raise ValueError('Rigid groups overlap')
        owner[group] = i
    count = len(groups)
    for i in np.flatnonzero(owner < 0):
        owner[i] = count; count += 1
    T = sparse.csr_matrix((np.ones(len(v)), (np.arange(len(v)), owner)), shape=(len(v), count))
    e = np.unique(np.sort(f[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1), axis=0)
    vectors = v[e[:, 1]]-v[e[:, 0]]; lengths = np.linalg.norm(vectors, axis=1)
    good = lengths > 1e-12
    e, vectors, lengths = e[good], vectors[good], lengths[good]
    C = T[e[:, 1]]-T[e[:, 0]]
    H = C.T@C+sparse.eye(count)*.03
    P = sparse.triu(sparse.kron(H, sparse.eye(3)), format='csc')
    samples = [sparse.eye(len(v), format='csr')]
    for weights in [[1/3, 1/3, 1/3], [.5, .5, 0], [0, .5, .5], [.5, 0, .5]]:
        samples.append(sparse.csr_matrix((np.tile(weights, len(f)),
                      (np.repeat(np.arange(len(f)), 3), f.ravel())), shape=(len(f), len(v))))
    S = sparse.vstack(samples, format='csr'); W = S@T; points = S@v
    delta = np.zeros((count, 3)); records = []
    fixed = np.unique(owner[v[:, 1] >= fixed_above])
    active_samples=np.ones(S.shape[0],bool)
    if active_raw is not None:
        selected=np.asarray(active_raw,int)
        if selected.ndim!=1 or not len(selected)or np.any(selected<0)or np.any(selected>=len(v)):
            raise ValueError('Explicit nonempty original active controls required')
        active_groups=np.unique(owner[selected])
        fixed=np.union1d(fixed,np.setdiff1d(np.arange(count),active_groups))
        active_samples=np.asarray(W[:,active_groups].sum(1)).ravel()>0

    def directional(matrix, normals):
        a = matrix.tocoo()
        return sparse.coo_matrix(((a.data[:, None]*normals[a.row]).ravel(),
                (np.repeat(a.row, 3), (a.col[:, None]*3+np.arange(3)).ravel())),
                shape=(matrix.shape[0], 3*count)).tocsr()

    def gaps(p):
        near, _ = bone.closest_points(p)
        tensor = library().core.Tensor(((p-bone.origin)/bone.scale).astype(np.float32))
        inside = bone.scene.compute_occupancy(tensor, nthreads=2, nsamples=3).numpy() > .5
        d = p-near['points']; d[inside] *= -1
        length = np.linalg.norm(d, axis=1)
        if (length < 1e-12).any() and not coherent_bone_sections:
            raise ValueError('Ambiguous exact bone witness')
        bn = np.divide(d,length[:,None],out=np.zeros_like(d),where=length[:,None]>1e-12)
        bg = near['distance_m']*np.where(inside, -1., 1.)
        if coherent_bone_sections:
            # A rigid ring crossing a cortex cannot move in two opposed
            # directions. All its supported samples use one exterior tangent
            # proposal from the ring centre. Exact triangle replay is still
            # required; this local tangent is not a clearance certificate.
            current=v+T@delta
            centres=np.column_stack([np.bincount(owner,weights=current[:,j],minlength=count)
                                     /np.bincount(owner,minlength=count)for j in range(3)])
            ch,_=bone.closest_points(centres)
            ct=library().core.Tensor(((centres-bone.origin)/bone.scale).astype(np.float32))
            ci=bone.scene.compute_occupancy(ct,nthreads=2,nsamples=3).numpy()>.5
            cn=centres-ch['points'];cn[ci]*=-1
            clen=np.linalg.norm(cn,axis=1)
            fallback,_=indexed_feature_normals(bone,ch)
            cn=np.divide(cn,clen[:,None],out=fallback.copy(),where=clen[:,None]>1e-12)
            live=np.unique(W[bg<clearance_m+.00001].indices)
            affected=np.asarray(W[:,live].sum(1)).ravel()>0 if len(live)else np.zeros(len(p),bool)
            normal=W@cn;norm=np.linalg.norm(normal,axis=1)
            affected&=norm>1e-12
            normal[affected]/=norm[affected,None]
            bn[affected]=normal[affected]
            bg[affected]=np.einsum('ij,ij->i',p[affected]-(W@ch['points'])[affected],normal[affected])
        hit, _ = skin.closest_points(p)
        outward, rim = indexed_feature_normals(skin,hit)
        sn, sg, _ = inward_sheet_clearance(p,hit,outward,rim)
        return bn, bg, sn, sg

    I = sparse.eye(3*count, format='csr')
    bound = np.full((count, 3), maximum_motion_m/np.sqrt(3))
    bound[fixed] = 0.
    edge_A = directional(C, vectors/lengths[:, None])
    for iteration in range(maximum_iterations):
        if time.monotonic()-started >= wall_seconds:
            return v.copy(), dict(complete=False, status='NativeDeadline', runs=records)
        q = v+T@delta; shift = W@delta
        bn, bg, sn, sg = gaps(points+shift)
        if min(bg[active_samples].min(), sg[active_samples].min()) >= clearance_m-1e-8:
            break
        # A fixed 75% projection prevents legitimate tight turns. Keep a
        # positive local orientation only; actual triangles decide collisions.
        rows = [edge_A, I, -I]; rhs = [-.95*lengths, -bound.ravel(), -bound.ravel()]
        bands = [('edge', len(lengths), None), ('lower-motion', 3*count, None), ('upper-motion', 3*count, None)]
        for normal, gap in [(bn, bg), (sn, sg)]:
            active = (gap < maximum_motion_m+clearance_m)&active_samples
            rows.append(directional(W[active], normal[active]))
            rhs.append(clearance_m+.00001-gap[active]+
                       np.einsum('ij,ij->i', normal[active], shift[active]))
            bands.append(('bone' if normal is bn else 'skin', int(active.sum()), np.flatnonzero(active)))
        # Linearized face orientation guards use only original source faces.
        tri = q[f]; n = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        size = np.linalg.norm(n, axis=1)
        if (size <= 1e-20).any():
            raise ValueError('Collapsed source triangle')
        n /= size[:, None]
        gradients = np.stack([np.cross(n, tri[:, 2]-tri[:, 1]),
                              np.cross(tri[:, 2]-tri[:, 0], n),
                              np.cross(n, tri[:, 1]-tri[:, 0])], axis=1)
        support = np.linalg.norm(gradients, axis=2).sum(1)
        fa = sum((directional(T[f[:, k]], gradients[:, k]/support[:, None]) for k in range(3)))
        rows.append(fa); rhs.append(-.95*size/support+fa@delta.ravel())
        bands.append(('face-orientation', len(f), None))
        A = sparse.vstack(rows, format='csc'); b = np.concatenate(rhs)/.001
        settings = clarabel.DefaultSettings(); settings.verbose = False
        settings.max_iter = 60; settings.time_limit = min(2., max(.01, wall_seconds-(time.monotonic()-started)))
        settings.tol_gap_abs = settings.tol_gap_rel = settings.tol_feas = 1e-8
        solution = clarabel.DefaultSolver(P, np.zeros(3*count), -A, -b,
                   [clarabel.NonnegativeConeT(len(b))], settings).solve()
        record = dict(iteration=iteration, status=str(solution.status),
                      solverSeconds=float(solution.solve_time), sparseNonzeros=int(A.nnz),
                      controls=count, minimumBoneGap_m=float(bg.min()), minimumSkinGap_m=float(sg.min()))
        records.append(record)
        if str(solution.status) != 'Solved':
            dual = np.array(solution.z); leading = np.argsort(-np.abs(dual))[:12]
            witnesses = []
            for j in leading:
                offset = 0
                for tag, length, indices in bands:
                    if offset <= j < offset+length:
                        item = dict(kind=tag, row=int(j-offset), dual=float(dual[j]))
                        if indices is not None:
                            index = int(indices[j-offset]); item.update(sample=index, point_m=(points+shift)[index].tolist())
                        witnesses.append(item); break
                    offset += length
            record['leadingDualWitnesses'] = witnesses
            return v.copy(), dict(complete=False, status=str(solution.status), runs=records,
                                  maximumMotionBound_m=maximum_motion_m,
                                  elapsedSeconds=time.monotonic()-started)
        trial = np.array(solution.x).reshape(count, 3)*.001
        residual = float(np.max(np.concatenate(rhs)-A@trial.ravel(), initial=0))
        if residual > 1e-8 or not np.isfinite(trial).all():
            return v.copy(), dict(complete=False, status='ResidualFailed', residual_m=residual, runs=records)
        # Locked proximal controls must remain literal, including solver noise.
        trial[fixed] = 0.; delta = trial
    q = v+T@delta; _, bg, _, sg = gaps(S@q)
    complete = bool(min(bg[active_samples].min(), sg[active_samples].min()) >= clearance_m-1e-8 and time.monotonic()-started < wall_seconds)
    return q if complete else v.copy(), dict(complete=complete, status='Proposed' if complete else 'ClearanceUnresolved',
              minimumSampledBoneGap_m=float(bg[active_samples].min()), minimumSampledSkinGap_m=float(sg[active_samples].min()),
              maximumControlMotion_m=float(np.linalg.norm(delta, axis=1).max()), runs=records,
              nativeWallSeconds=wall_seconds, elapsedSeconds=time.monotonic()-started,
              openSkinContainmentNotEstablished=True, actualTriangleAuditRequired=True,
              coherentBoneSectionProposal=coherent_bone_sections,
              explicitActiveControlScope=active_raw is not None,
              checkedSamples=int(active_samples.sum()),uncheckedFixedSamples=int((~active_samples).sum()))
