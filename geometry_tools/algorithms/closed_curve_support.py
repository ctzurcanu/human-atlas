"""Closest witnesses on an explicit closed 3D curve, without filling its span."""
import numpy as np


def curve_support_planes(curve, margin_m, *, sides=8, axial_margin_m=.004):
    """Convex selection footprint as affine halfspaces about an actual rim.

    This is a spatial support over an existing tissue domain. It does not
    triangulate or fill the opening. The returned planes can also be used in
    spatial queries; mesh construction needs its own precision policy.
    """
    curve=np.asarray(curve,float)
    if curve.ndim!=2 or curve.shape[1]!=3 or len(curve)<3 or not np.isfinite(curve).all():
        raise ValueError('An explicit finite 3D rim is required')
    if not margin_m>0 or not axial_margin_m>0 or not isinstance(sides,int) or not 8<=sides<=32:
        raise ValueError('Positive metre margins and 8–32 support sides are required')
    center=curve.mean(0);_,singular,basis=np.linalg.svd(curve-center,full_matrices=False)
    if singular[1]<=np.finfo(float).eps*max(singular[0],1.):
        raise ValueError('A rim must span a surface rather than a line')
    angle=np.arange(sides)*2*np.pi/sides
    normals=np.r_[np.cos(angle)[:,None]*basis[0]+np.sin(angle)[:,None]*basis[1],basis[2:3],-basis[2:3]]
    offsets=((curve-center)@normals.T).max(0)+np.r_[np.full(sides,margin_m),axial_margin_m,axial_margin_m]
    return center,normals,offsets


def closest_closed_curve(points, curve, *, batch_size=512):
    points, curve = np.asarray(points, float), np.asarray(curve, float)
    if points.ndim != 2 or points.shape[1] != 3 or curve.ndim != 2 or curve.shape[1] != 3 or len(curve) < 3:
        raise ValueError('Finite 3D query points and at least three ordered curve controls are required')
    if not np.isfinite(points).all() or not np.isfinite(curve).all() or batch_size <= 0:
        raise ValueError('Finite coordinates and a positive batch size are required')
    delta = np.roll(curve, -1, axis=0)-curve
    length2 = np.einsum('ij,ij->i', delta, delta)
    if np.any(length2 == 0):
        raise ValueError('A closed curve must not repeat its last or a neighbouring control')
    result = dict(distance_m=np.empty(len(points)), points=np.empty_like(points),
                  edge_ids=np.empty(len(points), int), edge_fraction=np.empty(len(points)))
    for start in range(0, len(points), batch_size):
        q = points[start:start+batch_size]
        t = np.clip(np.einsum('nki,ki->nk', q[:, None]-curve, delta)/length2, 0., 1.)
        foot = curve[None]+t[..., None]*delta
        d2 = np.einsum('nki,nki->nk', q[:, None]-foot, q[:, None]-foot)
        edge = np.argmin(d2, axis=1); row = np.arange(len(q)); sl = slice(start, start+len(q))
        result['distance_m'][sl] = np.sqrt(d2[row, edge])
        result['points'][sl] = foot[row, edge]
        result['edge_ids'][sl] = edge
        result['edge_fraction'][sl] = t[row, edge]
    return result
