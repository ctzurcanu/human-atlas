"""Transport source sections along actual vascular links, pinning interfaces."""
import numpy as np
from .profile_chain_frames import branch_anchored_frames


def finite_diagnostic(value):
    """Retain infeasible-bound evidence without invalid JSON float literals."""
    if isinstance(value, dict):
        return {k: finite_diagnostic(v) for k, v in value.items()}
    if isinstance(value, list):
        return [finite_diagnostic(v) for v in value]
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return '+infinity' if value > 0 else '-infinity' if value < 0 else 'NaN'
    return value


def profile_chains(count, links, protected=()):
    adjacency = [set() for _ in range(count)]
    edges = set()
    for a, b in np.asarray(links, int).reshape(-1, 2):
        if min(a, b) < 0 or max(a, b) >= count or a == b:
            raise ValueError('Invalid explicit profile link')
        adjacency[a].add(int(b)); adjacency[b].add(int(a)); edges.add(tuple(sorted((int(a), int(b)))))
    ends = {i for i, n in enumerate(adjacency) if len(n) != 2} | set(map(int, protected))
    seen = set(); chains = []
    for a in sorted(ends):
        for b in sorted(adjacency[a]):
            if tuple(sorted((a, b))) in seen:
                continue
            chain = [a, b]; seen.add(tuple(sorted((a, b))))
            while chain[-1] not in ends:
                choices = adjacency[chain[-1]]-{chain[-2]}
                if len(choices) != 1:
                    raise ValueError('Ambiguous source chain topology')
                c = next(iter(choices)); edge = tuple(sorted((chain[-1], c)))
                if edge in seen:
                    raise ValueError('Repeated source edge in an anchored chain')
                seen.add(edge); chain.append(c)
            chains.append(chain)
    return chains, sorted(edges-seen)


def transport_linked_profiles(source, target, rotations, links, protected=()):
    p, q, R = np.asarray(source, float), np.asarray(target, float), np.asarray(rotations, float)
    if p.shape != q.shape or R.shape != (len(p), 3, 3):
        raise ValueError('One source/target point and proper pose per profile required')
    chains, closed = profile_chains(len(p), links, protected)
    out = R.copy(); records = []
    for chain in chains:
        ids = np.array(chain, int)
        if len(ids) < 3:
            continue
        posed = branch_anchored_frames(p[ids], q[ids], R[ids[0]], R[ids[-1]])
        out[ids] = posed
        records.append(dict(profiles=chain, endpointRotationsLiteral=True))
    fixed = np.array(sorted(set(map(int, protected))), int)
    if len(fixed) and not np.array_equal(out[fixed], R[fixed]):
        raise ValueError('Protected interface frame changed')
    return out, dict(chains=records, unresolvedClosedEdges=closed,
                     allSourceCentresRetained=True, allTargetCentresRetained=True,
                     protectedInterfaceFramesLiteral=True)


def round_linked_courses(source, target, links, support, index, protected=(), *,
                         reference_rotations=None, minimum_advance_fraction=.25,
                         smoothing_weight=2., maximum_motion_m=.003, iterations=2):
    """Smooth each declared source chain with node AND span halfspaces.

    Conic solves are bounded. Failed chains remain literal diagnostics, with
    failure evidence, and cannot authorize installation. Halfspaces are local
    native proposals; the resulting actual triangles need separate review.
    """
    from .section_course_socp import solve_weighted_section_socp
    from geometry_tools.open3d_backend import library
    from .tube_sections import arclength
    p, q = np.asarray(source, float), np.asarray(target, float).copy()
    radius = np.asarray(support, float)
    reference = np.asarray(reference_rotations, float) if reference_rotations is not None else None
    chains, closed = profile_chains(len(p), links, protected)
    evidence = []
    for chain in chains:
        ids = np.array(chain, int); n = len(ids)
        if n < 3:
            continue
        original = q[ids].copy()
        h = np.diff(arclength(p[ids])); typical = np.median(h)
        D = np.zeros((n-2, n))
        for i in range(n-2):
            D[i, i:i+3] = [typical/h[i], -typical*(1/h[i]+1/h[i+1]), typical/h[i+1]]
        H = np.eye(n)+smoothing_weight*(D.T@D)
        free, ends = np.arange(1, n-1), np.array([0, n-1])
        base = original.copy()
        base[free] = np.linalg.solve(H[np.ix_(free, free)], original[free]-H[np.ix_(free, ends)]@original[ends])
        initial_motion = float(np.linalg.norm(base-original, axis=1).max())
        if initial_motion > maximum_motion_m:
            evidence.append(dict(profiles=chain, status='SmoothingTrustBoundExceeded', motion_m=initial_motion))
            continue
        K = np.zeros((n, n)); K[np.ix_(free, free)] = np.linalg.inv(H[np.ix_(free, free)])
        edge_weights = np.zeros((n-1, n))
        edge_weights[np.arange(n-1), np.arange(n-1)] = -1.
        edge_weights[np.arange(n-1), np.arange(1, n)] = 1.
        source_edges = np.diff(p[ids], axis=0)
        if reference is None:
            direction = np.diff(original, axis=0)
        else:
            direction = np.einsum('nij,nj->ni', reference[ids[:-1]], source_edges)
        advance_normals = direction/np.linalg.norm(direction, axis=1)[:, None]
        advance_offsets = minimum_advance_fraction*np.linalg.norm(source_edges, axis=1)-np.einsum(
            'ij,ij->i', advance_normals, edge_weights@base)
        W = [np.eye(n)]
        sup = [radius[ids]]
        for fraction in [.25, .5, .75]:
            w = np.zeros((n-1, n))
            w[np.arange(n-1), np.arange(n-1)] = 1-fraction
            w[np.arange(n-1), np.arange(1, n)] = fraction
            W.append(w); sup.append(np.maximum(radius[ids[:-1]], radius[ids[1:]]))
        W = np.concatenate(W); sup = np.concatenate(sup)
        current = base.copy()
        weights, normals, offsets, runs = [edge_weights], [advance_normals], [advance_offsets], []
        success = True
        for _ in range(iterations):
            points = W@current
            closest, _ = index.closest_points(points)
            local = ((points-index.origin)/index.scale).astype(np.float32)
            inside = index.scene.compute_occupancy(library().core.Tensor(local), nthreads=2, nsamples=3).numpy() > .5
            direction = points-closest['points']; direction[inside] *= -1
            length = np.linalg.norm(direction, axis=1)
            if np.any(length < 1e-12):
                success = False; runs.append(dict(status='AmbiguousRecipientNormal')); break
            normal = direction/length[:, None]
            signed = closest['distance_m']*np.where(inside, -1., 1.)
            rhs = sup+.00025-signed+np.einsum('ij,ij->i', normal, W@(current-base))
            weights.append(W); normals.append(normal); offsets.append(rhs)
            delta, solve = solve_weighted_section_socp(K, np.concatenate(weights), np.concatenate(normals),
                np.concatenate(offsets), maximum_motion_m=maximum_motion_m,
                residual_m=1e-8, time_limit_seconds=.25, maximum_iterations=60)
            runs.append(finite_diagnostic(solve))
            if delta is None:
                success = False; break
            current = base+delta
        if success:
            current[ends] = original[ends]
            q[ids] = current
        evidence.append(dict(profiles=chain, status='Proposed' if success else 'ConstraintProposalFailed',
                             endpointsLiteral=True, solveRuns=runs,
                             maximumMotion_m=float(np.linalg.norm(current-original, axis=1).max())))
    fixed = np.array(sorted(set(map(int, protected))), int)
    if len(fixed) and not np.array_equal(q[fixed], np.asarray(target)[fixed]):
        raise ValueError('Protected source branch endpoint moved')
    return q, dict(chains=evidence, unresolvedClosedEdges=closed,
                   sampledNodeAndSpanHalfspacesOnly=True, wholeSurfaceClearanceCertified=False)
