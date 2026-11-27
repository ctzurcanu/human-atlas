"""Bound piecewise-linear source deformation without changing source profiles.

Each returned row is one tangent derivative of the displacement field.
Bounding both derivative-vector norms by g/sqrt(2) bounds the full 3x2
gradient by g. With g<1 the tangent projection cannot invert. The actual
source triangles and normals still require independent replay.
"""
import numpy as np


def facet_gradient_rows(mesh, vertex_weights):
    w = np.asarray(vertex_weights, float)
    if w.ndim != 2 or len(w) != len(mesh.vertices) or not np.isfinite(w).all():
        raise ValueError('Finite source weights required for every vertex')
    tri = mesh.triangles
    e = tri[:, 1] - tri[:, 0]
    f = tri[:, 2] - tri[:, 0]
    length = np.linalg.norm(e, axis=1)
    cross = np.cross(e, f)
    area2 = np.linalg.norm(cross, axis=1)
    if np.any(length == 0) or np.any(area2 == 0):
        raise ValueError('Degenerate source facet cannot define a gradient')
    tangent = e / length[:, None]
    along = np.einsum('ij,ij->i', tangent, f)
    height = area2 / length
    dw1 = w[mesh.faces[:, 1]] - w[mesh.faces[:, 0]]
    dw2 = w[mesh.faces[:, 2]] - w[mesh.faces[:, 0]]
    gx = dw1 / length[:, None]
    gy = (dw2 - along[:, None] * gx) / height[:, None]
    return np.stack([gx, gy], axis=1)


def compact_gradient_rows(rows):
    """Retain strongest literal opposite-pair row; keep general rows intact."""
    rows = np.asarray(rows, float).reshape(-1, rows.shape[-1])
    pairs = {}
    general = []
    for row in rows:
        ids = np.flatnonzero(row)
        if not len(ids):
            continue
        if len(ids) == 2 and row[ids[0]] == -row[ids[1]]:
            key = tuple(ids)
            if key not in pairs or abs(row[ids[0]]) > abs(pairs[key][ids[0]]):
                pairs[key] = row.copy()
        else:
            general.append(row)
    out = list(pairs.values()) + general
    return np.unique(np.asarray(out).reshape(-1, rows.shape[-1]), axis=0)


def gradient_evidence(mesh, weights, displacements):
    grad = facet_gradient_rows(mesh, weights) @ np.asarray(displacements)
    return dict(maximumGradientFrobenius=float(np.linalg.norm(grad, axis=(1, 2)).max(initial=0)),
                maximumGradientColumnNorm=float(np.linalg.norm(grad, axis=2).max(initial=0)))
