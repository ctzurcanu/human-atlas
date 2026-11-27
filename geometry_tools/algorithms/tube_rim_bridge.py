"""Join two explicitly paired open sheet boundaries without a lumen disk.

Dynamic programming chooses a minimum-cost monotone annular triangulation,
enumerating cyclic phases. Input controls remain literal and may be nonplanar.
This supplies combinatorics only; whole geometry and patency require checking.
"""
from collections import Counter
import numpy as np

from .annular_pipe_envelopes import boundary_loops


def _edge_chain(faces):
    counts = Counter((int(a), int(b)) for t in faces for a, b in zip(t, np.roll(t, -1)))
    return Counter({e: n - counts[e[::-1]] for e, n in counts.items() if n > counts[e[::-1]]})


def bridge_boundaries(positions, first_boundary, second_boundary, *, align_projected_phase=True):
    """Faces whose boundary opposes both original oriented sheet boundaries."""
    p = np.asarray(positions, float)
    a = np.asarray(first_boundary, int)[::-1]; original_b = np.asarray(second_boundary, int)
    if min(len(a), len(original_b)) < 3 or len(set(a)) != len(a) or len(set(original_b)) != len(original_b):
        raise ValueError('Two genuine simple boundary cycles are required')
    if set(a) & set(original_b):
        raise ValueError('Already touching boundaries require an explicit shared-mouth operation')
    n, m = len(a), len(original_b); best = None
    ca,cb=p[a].mean(0),p[original_b].mean(0);axis=cb-ca
    if np.linalg.norm(axis)==0:raise ValueError('Distinct section centres are required for a rim loft')
    axis/=np.linalg.norm(axis);u=np.cross(axis,np.eye(3)[np.argmin(abs(axis))]);u/=np.linalg.norm(u);w=np.cross(axis,u)
    def angle(ids,center):
        q=p[ids]-center
        return np.arctan2(q@w,q@u)
    aa=angle(a,ca)
    for phase in range(m):
        b = np.roll(original_b, phase)
        bb=angle(b,cb)
        def phase_ok(i,j):
            return not align_projected_phase or abs(np.angle(np.exp(1j*(aa[i%n]-bb[j%m]))))<np.pi/2
        if not phase_ok(0,0):continue
        cost = np.full((n + 1, m + 1), np.inf); cost[0, 0] = 0.
        previous = np.full((n + 1, m + 1), -1, int)
        for i in range(n + 1):
            for j in range(m + 1):
                if not np.isfinite(cost[i, j]):
                    continue
                options = []
                if i < n:
                    options.append((i + 1, j, 0, [a[i % n], a[(i + 1) % n], b[j % m]]))
                if j < m:
                    options.append((i, j + 1, 1, [a[i % n], b[(j + 1) % m], b[j % m]]))
                for ii, jj, step, tri in options:
                    # A whole circuit of one ring while stalled on the other
                    # reuses a cross-seam edge four times. Keep relative phase
                    # within an open half turn; cyclic starts are all enumerated.
                    if abs(ii / n - jj / m) >= .5:
                        continue
                    if not phase_ok(ii,jj):continue
                    q = p[tri]; area2 = np.linalg.norm(np.cross(q[1] - q[0], q[2] - q[0]))
                    if not area2:
                        continue
                    lengths = np.sum((q - np.roll(q, -1, axis=0)) ** 2, axis=1)
                    # Area plus a scale-consistent penalty for skinny triangles.
                    value = cost[i, j] + area2 + .05 * lengths.max() ** 2 / area2
                    if value < cost[ii, jj]:
                        cost[ii, jj] = value; previous[ii, jj] = step
        if not np.isfinite(cost[n, m]):
            continue
        faces = []; i, j = n, m
        while i or j:
            step = previous[i, j]
            if step == 0:
                i -= 1; faces.append([a[i % n], a[(i + 1) % n], b[j % m]])
            elif step == 1:
                j -= 1; faces.append([a[i % n], b[(j + 1) % m], b[j % m]])
            else:
                raise ValueError('Incomplete annular triangulation')
        candidate = (float(cost[n, m]), phase, np.asarray(faces, int))
        if best is None or candidate[:2] < best[:2]:
            best = candidate
    if best is None:
        raise ValueError('No nondegenerate monotone annulus exists')
    _, phase, f = best
    edges=np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]])
    unique,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
    if np.any(counts>2) or len(np.unique(f))-len(unique)+len(f)!=0:
        raise ValueError('The proposed strip is not a manifold annulus')
    expected = Counter((int(b), int(a)) for ring in [first_boundary, second_boundary]
                       for a, b in zip(ring, np.roll(ring, -1)))
    if _edge_chain(f) != expected or len(boundary_loops(f)) != 2:
        raise ValueError('The proposed strip does not oppose both source boundaries exactly')
    return f, dict(first_controls=n, second_controls=m, triangles=len(f), cyclic_phase=phase,
                   source_controls_unchanged=True, exact_opposed_boundary_chain=True,
                   nonpinching_relative_phase_band=True,
                   projected_phase_alignment=bool(align_projected_phase),
                   lumen_disk_created=False, self_intersections_not_yet_checked=True)


def end_annulus(wall_faces, inner_faces, inner_loop):
    """Propose the actual one-ring end annulus from explicit lining topology.

    The two boundary cycles and genus-zero annulus are mandatory. The caller
    still reviews whether this band is the correct anatomical joining rim.
    """
    wall = np.asarray(wall_faces, int); inside = {tuple(sorted(t)) for t in inner_faces}
    candidates = np.array([i for i, t in enumerate(wall) if tuple(sorted(t)) not in inside and
                           any(int(k) in set(map(int, inner_loop)) for k in t)], int)
    band = wall[candidates]; loops = boundary_loops(band)
    if len(loops) != 2 or sum(set(q) == set(inner_loop) for q in loops) != 1:
        raise ValueError('The proposed end band is not bounded by the genuine inner lip and one outer lip')
    outer = next(q for q in loops if set(q) != set(inner_loop))
    # The outer sheet has the opposite boundary to its removed end annulus.
    return candidates, outer[::-1]
