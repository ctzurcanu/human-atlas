"""Material witnesses for proper intersections of literal source triangles.

The same plane/interval geometry as the existing crossing check supplies a
barycentric point. Parallel/coplanar contacts remain explicitly unsupported;
they must never be silently discarded from a final independent crossing gate.
"""
import numpy as np


def barycentric(point, triangle):
    t = np.asarray(triangle, float); p = np.asarray(point, float)
    edges = np.array([t[1] - t[0], t[2] - t[0]])
    uv = np.linalg.solve(edges @ edges.T, edges @ (p - t[0]))
    weights = np.r_[1. - uv.sum(), uv]
    error = float(np.linalg.norm(weights @ t - p))
    bound = 256 * np.finfo(float).eps * max(1., np.abs(t).max())
    if error > bound or weights.min() < -1e-8:
        raise ValueError('The contact witness is not on the actual source triangle')
    return weights


def intersection_witness(first, second):
    a = np.asarray(first, float); b = np.asarray(second, float)
    origin = a[0].copy(); aa = a - origin; bb = b - origin
    na = np.cross(aa[1] - aa[0], aa[2] - aa[0])
    nb = np.cross(bb[1] - bb[0], bb[2] - bb[0])
    na /= np.linalg.norm(na); nb /= np.linalg.norm(nb)
    line = np.cross(na, nb); length = np.linalg.norm(line)
    if length <= 1e-7:
        raise ValueError('Coplanar contact requires a polygonal witness')
    line /= length
    da = (aa - bb[0]) @ nb; db = (bb - aa[0]) @ na
    def interval(t, signed):
        points = []
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            if abs(signed[i]) <= 1e-14:
                points.append(t[i])
            if signed[i] * signed[j] < 0.:
                points.append(t[i] + (t[j] - t[i]) * signed[i] / (signed[i] - signed[j]))
        if len(points) < 2:
            raise ValueError('A proper source intersection is required')
        points = np.asarray(points); projection = points @ line
        return points[int(projection.argmin())], points[int(projection.argmax())]
    ca = np.asarray(interval(aa, da)); cb = np.asarray(interval(bb, db))
    ia, ib = ca @ line, cb @ line
    low, high = max(ia[0], ib[0]), min(ia[1], ib[1])
    if high <= low:
        raise ValueError('Source triangle intersection intervals do not overlap')
    p = ca[0] + line * ((low + high) * .5 - ia[0]) + origin
    return p, barycentric(p, a), barycentric(p, b)
