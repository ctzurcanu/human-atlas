"""Two-way PL surface bound through a shared injective planar projection.

Caller supplies two complete triangle disks with a literal common boundary.
The numerical overlay checks projected coverage and non-overlap; it is not an
exact-arithmetic or anatomical acceptance certificate.
"""
from collections import Counter

import numpy as np

from .tube_continuity import continuity_audit


def _cross(a, b):
    return a[0]*b[1]-a[1]*b[0]


def _area(p):
    return abs(sum(_cross(a, b) for a, b in zip(p, np.roll(p, -1, axis=0))))*.5 if len(p)>2 else 0.


def _intersection(first, second):
    polygon = list(first)
    if _cross(second[1]-second[0], second[2]-second[0]) < 0:
        second = second[::-1]
    for a, b in zip(second, np.roll(second, -1, axis=0)):
        output = []
        if not polygon:
            break
        previous = polygon[-1]
        previous_side = _cross(b-a, previous-a)
        for current in polygon:
            side = _cross(b-a, current-a)
            if (side >= 0) != (previous_side >= 0):
                output.append(previous+(current-previous)*previous_side/(previous_side-side))
            if side >= 0:
                output.append(current)
            previous, previous_side = current, side
        polygon = output
    return np.asarray(polygon).reshape(-1, 2)


def projected_patch_bound(vertices, old_faces, new_faces, *, maximum_distance_m):
    points = np.asarray(vertices, float)
    old, new = [np.asarray(f, int) for f in [old_faces, new_faces]]
    boundary = lambda f: {e for e, n in Counter(tuple(sorted(map(int, (a, b))))
        for t in f for a, b in zip(t, np.roll(t, -1))).items() if n == 1}
    if boundary(old) != boundary(new):
        raise ValueError('A literal common source boundary is required')
    for faces in [old, new]:
        if not continuity_audit(points, faces, expected_components=1, expected_mouths=1,
                                 expected_euler=1)['topological_continuity_pass']:
            raise ValueError('Each source patch must be one oriented triangle disk')
    ids = np.unique(np.r_[old.ravel(), new.ravel()])
    origin = points[ids].mean(0)
    scale = float(np.ptp(points[ids], axis=0).max())
    _, _, axes = np.linalg.svd(points[ids]-origin, full_matrices=False)
    local = (points-origin)@axes.T/scale
    patches = [local[f] for f in [old, new]]
    fields, areas = [], []
    for patch in patches:
        signed = np.array([_cross(t[1, :2]-t[0, :2], t[2, :2]-t[0, :2]) for t in patch])
        if np.min(abs(signed)) < 1e-12 or not (np.all(signed>0) or np.all(signed<0)):
            raise ValueError('Every source triangle must have the same nonsingular projected orientation')
        fields.append([np.linalg.solve(np.c_[t[:, :2], np.ones(3)], t[:, 2]) for t in patch])
        areas.append(abs(signed)*.5)
        for i in range(len(patch)):
            for j in range(i):
                if _area(_intersection(patch[i, :, :2], patch[j, :, :2])) > 1e-11:
                    raise ValueError('Projected source patch has overlapping triangle interiors')
    coverage = [np.zeros(len(p)) for p in patches]
    gap = 0.
    for i, a in enumerate(patches[0]):
        for j, b in enumerate(patches[1]):
            overlap = _intersection(a[:, :2], b[:, :2])
            area = _area(overlap)
            coverage[0][i] += area
            coverage[1][j] += area
            if len(overlap):
                difference = np.c_[overlap, np.ones(len(overlap))]@(fields[0][i]-fields[1][j])
                gap = max(gap, float(np.max(abs(difference))))
    for expected, actual in zip(areas, coverage):
        if not np.allclose(actual, expected, atol=1e-11, rtol=1e-10):
            raise ValueError('Projected triangle overlays do not cover both entire source patches')
    gradient = max(np.linalg.norm(field[:2]) for patch in fields for field in patch)
    margin = max(2e-12, 64*np.finfo(float).eps*scale*(1+gradient))
    bound = gap*scale+margin
    if bound > maximum_distance_m:
        raise ValueError('The complete two-way projected surface bound exceeds the declared distance')
    return dict(complete=True,passed=True,maximumTwoWaySurfaceDistanceUpperBound_m=bound,
                maximumOverlayHeightDifference_m=gap*scale,numericalMargin_m=margin,
                projectedPatchesInjective=True,completeOverlayCoverageChecked=True,
                literalBoundaryPreserved=True,exactArithmetic=False,anatomicalAcceptance=False)
