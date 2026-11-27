"""Float64 ray probes against literal triangles; unhit directions stay unknown.

This measures geometry only. It does not establish clinical wall thickness,
tissue occupancy, correspondence, or a clear route for a finite-radius tube.
"""
import numpy as np


def literal_surface_rays(triangles, origins, directions, *, minimum_incidence=1e-8):
    t = np.asarray(triangles, float)
    o = np.asarray(origins, float)
    d = np.asarray(directions, float)
    if t.ndim != 3 or t.shape[1:] != (3, 3) or o.ndim != 2 or o.shape[1] != 3 or d.shape != o.shape:
        raise ValueError('Literal triangles and corresponding finite 3D rays required')
    if not np.isfinite(t).all() or not np.isfinite(o).all() or not np.isfinite(d).all():
        raise ValueError('Finite geometry required')
    length = np.linalg.norm(d, axis=1)
    if np.any(length == 0) or not 0 < minimum_incidence < 1:
        raise ValueError('Nonzero directions and explicit conditioning bound required')
    d = d/length[:, None]
    center = (t.min(axis=(0, 1))+t.max(axis=(0, 1)))*.5
    v = t-center
    e1, e2 = v[:, 1]-v[:, 0], v[:, 2]-v[:, 0]
    normal = np.cross(e1, e2); area = np.linalg.norm(normal, axis=1)
    if np.any(area == 0):
        raise ValueError('Nondegenerate source triangles required')
    hit = np.zeros(len(o), bool); distance = np.full(len(o), np.inf)
    faces = np.full(len(o), -1, int); counts = np.zeros(len(o), int)
    incidence = np.zeros(len(o)); bary = np.full((len(o), 3), np.nan)
    eps = 128*np.finfo(float).eps
    for i, (origin, direction) in enumerate(zip(o-center, d)):
        p = np.cross(direction, e2); det = np.einsum('ij,ij->i', e1, p)
        good = abs(det)/area >= minimum_incidence
        inv = np.divide(1., det, out=np.zeros_like(det), where=good)
        delta = origin-v[:, 0]
        u = np.einsum('ij,ij->i', delta, p)*inv
        q = np.cross(delta, e1)
        w = (q@direction)*inv
        ray_distance = np.einsum('ij,ij->i', e2, q)*inv
        valid = good & (u >= -eps) & (w >= -eps) & (u+w <= 1+eps) & (ray_distance >= 0)
        counts[i] = int(valid.sum())
        if valid.any():
            f = int(np.argmin(np.where(valid, ray_distance, np.inf)))
            hit[i] = True; distance[i] = ray_distance[f]; faces[i] = f
            incidence[i] = abs(det[f])/area[f]
            bary[i] = [1-u[f]-w[f], u[f], w[f]]
    return dict(hit=hit, distance_m=distance, face_ids=faces, triangle_hit_counts=counts,
                incidence_cosine=incidence, barycentric=bary), dict(complete=True,
        method='Float64 literal Moller-Trumbore ray/triangle probes, localized coordinates',
        triangleCount=len(t), rayCount=len(o), sourceGeometryChanged=False,
        minimumIncidenceCosine=minimum_incidence, barycentricSlack=eps,
        unhitDirectionsRemainUnknown=True, countsAreTrianglesNotDistinctTissueWalls=True,
        clinicalThicknessEstablished=False, exactArithmeticClaim=False)
