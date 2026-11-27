"""Conform a boundary chain after exact coordinate serialization.

This moves no point and fills no positive-area hole. It cancels oppositely
oriented identical triangles and resolves zero-area three-edge cracks by
subdividing the existing long boundary edge at its existing collinear point.
True openings and remaining defects are returned for the caller's strict gate.
"""
from collections import Counter
import numpy as np
from geometry_tools.algorithms.finite_vessel_wall import oriented_boundary_loops


def _oriented(triangle):
    a = tuple(map(int, triangle))
    return min(a, a[1:] + a[:1], a[2:] + a[:2])


def conform_encoded_boundary(vertices, faces, parents=None):
    v = np.asarray(vertices, float)
    original = np.asarray(faces, int)
    ancestry = np.arange(len(original)) if parents is None else np.asarray(parents, int)
    if ancestry.shape != (len(original),):
        raise ValueError('Every source triangle needs an ancestry entry')
    if len(original) == 0 or original.min() < 0 or original.max() >= len(v):
        raise ValueError('A nonempty indexed boundary is required')
    # Identifications are literal coordinate equality, never a distance weld.
    _, first, inverse = np.unique(v, axis=0, return_index=True, return_inverse=True)
    mapped = first[inverse][original]
    triangles = v[mapped]
    n = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    zero = np.all(n == 0., axis=1)
    buckets = {}
    for i in np.flatnonzero(~zero):
        buckets.setdefault(_oriented(mapped[i]), []).append(int(i))
    emitted, parent_ids, handled = [], [], set()
    opposed = 0
    for key, ids in buckets.items():
        if key in handled: continue
        reverse = _oriented((key[0], key[2], key[1]))
        other = buckets.get(reverse, [])
        handled.update((key, reverse))
        count = len(ids)-len(other)
        opposed += 2*min(len(ids), len(other))
        if abs(count) > 1:
            raise ValueError('A same-oriented duplicate boundary has nonunit multiplicity')
        if count:
            selected = ids[0] if count > 0 else other[0]
            emitted.append(mapped[selected]); parent_ids.append(ancestry[selected])
    f = np.asarray(emitted, int).reshape(-1, 3)
    parent_ids = np.asarray(parent_ids, int)
    loops = oriented_boundary_loops(f)
    refinements, unmodified = [], []
    for loop in loops:
        if len(loop) != 3:
            unmodified.append(list(map(int, loop))); continue
        points = v[loop]
        if np.any(np.cross(points[1]-points[0], points[2]-points[0]) != 0.):
            unmodified.append(list(map(int, loop))); continue
        pairs = [(int(loop[i]), int(loop[(i+1)%3])) for i in range(3)]
        a, c = max(pairs, key=lambda e: float(np.linalg.norm(v[e[1]]-v[e[0]])))
        b = next(int(i) for i in loop if i not in (a, c))
        direction = v[c]-v[a]
        fraction = float((v[b]-v[a])@direction/(direction@direction))
        if not 0. < fraction < 1.:
            raise ValueError('The existing boundary point is not strictly on its long edge')
        incident = np.flatnonzero(np.any(f == a, axis=1) & np.any(f == c, axis=1))
        if len(incident) != 1:
            raise ValueError('The long boundary edge must have exactly one incident triangle')
        face_id = int(incident[0]); old = f[face_id]
        k = next(k for k in range(3) if {int(old[k]), int(old[(k+1)%3])} == {a,c})
        x, y, d = map(int, np.roll(old, -k))
        replacement = np.array([[x,b,d], [b,y,d]], int)
        triangle = v[replacement]
        normal = np.cross(triangle[:,1]-triangle[:,0], triangle[:,2]-triangle[:,0])
        old_normal = np.cross(v[old[1]]-v[old[0]], v[old[2]]-v[old[0]])
        if np.any(normal@old_normal <= 0.):
            raise ValueError('Collinear subdivision changed the incident face orientation')
        parent = int(parent_ids[face_id])
        f[face_id] = replacement[0]; f = np.r_[f, replacement[1:]]
        parent_ids = np.r_[parent_ids, parent]
        refinements.append(dict(longEdge=[a,c], existingMiddleControl=b, sourceParent=parent,
            fraction=fraction, addedTrianglesFromSubdivision=1,
            newlyCoveredAreaMetersSquared=0., insertedSourceControls=0))
    return f, parent_ids, dict(exactCoordinateIdentificationsOnly=True, coordinatesMoved=False,
        zeroAreaTrianglesRemoved=int(zero.sum()), opposedCoincidentTrianglesCancelled=opposed,
        collinearBoundaryRefinements=refinements, trueOrUnresolvedBoundaryLoops=unmodified,
        anatomicalOpeningCapsAdded=0, parentTriangleChartsRetained=True)
