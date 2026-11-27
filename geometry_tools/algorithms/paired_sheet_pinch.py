"""Restore one collapsed point between explicitly corresponding source sheets.

This is not a wall-thickness generator or proximity weld. The caller supplies
the two complete opposite fans and their already corresponding neighbors.
"""
import numpy as np


def restore_paired_pinch(vertices, faces, outer, inner, neighbors, *, maximum_motion_m=3e-5):
    v = np.asarray(vertices, float)
    f = np.asarray(faces, int)
    pairs = np.asarray(neighbors, int)
    assert pairs.ndim == 2 and pairs.shape[1] == 2 and len(pairs) >= 3
    assert len(np.unique(pairs)) == pairs.size
    assert outer != inner and np.array_equal(v[outer], v[inner]), 'A literal collapsed source pair is required'
    assert outer not in pairs and inner not in pairs
    a = f[np.any(f == outer, axis=1)]
    b = f[np.any(f == inner, axis=1)]
    assert len(a) == len(b) == len(pairs), 'Only complete interior disk fans are supported'
    assert set(a.ravel()) == {outer, *pairs[:, 0]} and set(b.ravel()) == {inner, *pairs[:, 1]}
    correspondence = {outer: inner, **dict(pairs.tolist())}
    def cyclic(tri):
        return min(tuple(np.roll(tri, k)) for k in range(3))
    assert {cyclic([correspondence[int(x)] for x in t][::-1]) for t in a} == {cyclic(t) for t in b}, 'Opposite complete source fans must correspond'
    for fan, centre in [(a, outer), (b, inner)]:
        edges = [tuple(sorted(set(t)-{centre})) for t in fan]
        assert all(len(e) == 2 for e in edges)
        ids, degree = np.unique(edges, return_counts=True)
        assert len(ids) == len(pairs) and np.all(degree == 2)
        visited = {int(ids[0])}
        while True:
            reached = visited | {x for e in edges if any(y in visited for y in e) for x in e}
            if reached == visited:
                break
            visited = reached
        assert len(visited) == len(ids), 'The fan link must be one simple cycle'
    delta = v[pairs[:, 1]] - v[pairs[:, 0]]
    spacing = np.linalg.norm(delta, axis=1)
    assert np.all(spacing > 1e-9) and np.all(spacing <= maximum_motion_m)
    distances = np.linalg.norm(v[pairs[:, 0]]-v[outer], axis=1)
    assert np.all(distances > maximum_motion_m)
    weights = 1 / distances
    offset = np.average(delta, axis=0, weights=weights)
    normal = np.cross(v[a[:, 1]]-v[a[:, 0]], v[a[:, 2]]-v[a[:, 0]]).sum(0)
    normal /= np.linalg.norm(normal)
    assert offset @ normal < -1e-9, 'The restored inner point must stay inside its explicit outer fan'
    assert 1e-9 < np.linalg.norm(offset) <= maximum_motion_m
    proposed = v.copy()
    proposed[inner] += offset
    assert np.count_nonzero(np.any(proposed != v, axis=1)) == 1
    return proposed, dict(outer=int(outer), inner=int(inner), neighbors=pairs.tolist(),
                          offset_m=offset.tolist(), maximumMotion_m=float(np.linalg.norm(offset)),
                          neighboringSpacing_m=spacing.tolist(), weights=weights.tolist(),
                          oppositeCompleteFans=True, onlyOneExistingControlMoved=True,
                          facesAdded=0, facesRemoved=0, anatomicalThicknessEstablished=False)
