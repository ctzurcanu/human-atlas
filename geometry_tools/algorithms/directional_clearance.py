"""Choose a consistent exit side of a closed obstacle for a fitting proposal.

Nearest-surface normals on opposite sides of an intrusion can ask a finite
tissue sheet to split around the obstacle. A declared exit direction instead
provides tangent planes on the same visible side. These planes are proposal
constraints only: full occupied-volume and surface checks remain mandatory.
"""
from __future__ import annotations

import numpy as np


def directional_exit_planes(index, points, direction, *, minimum_facing_cosine=.05):
    """Ray-cast from outside the target bounds, back towards each query.

    ``index`` is an Open3D SurfaceIndex of an explicitly closed obstacle.
    All inputs and returned points use its original world coordinate frame.
    Misses and grazing facets are excluded explicitly; no plane is invented.
    """
    from geometry_tools.open3d_backend import library

    points = np.asarray(points, dtype=np.float64)
    direction = np.asarray(direction, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError('Finite N by 3 query points required')
    length = np.linalg.norm(direction)
    if direction.shape != (3,) or not np.isfinite(length) or length == 0:
        raise ValueError('An explicit finite nonzero exit direction is required')
    if not 0 < minimum_facing_cosine <= 1:
        raise ValueError('A positive facing threshold at most one is required')
    direction = direction / length
    if not len(points):
        return dict(valid=np.zeros(0, bool), points=points.copy(),
                    normals=points.copy(), gap_m=np.zeros(0), face_ids=np.zeros(0, int))

    # A bounding-box support plane guarantees every origin is outside, even
    # for queries that themselves lie beyond or behind the obstacle.
    local_points = (points - index.origin) / index.scale
    low = (index.mesh.vertices.min(0) - index.origin) / index.scale
    high = (index.mesh.vertices.max(0) - index.origin) / index.scale
    support = np.where(direction >= 0, high, low) @ direction
    offset = np.maximum(support - local_points @ direction, 0) + 1.
    origins = local_points + offset[:, None] * direction
    rays = np.c_[origins, np.broadcast_to(-direction, points.shape)].astype(np.float32)
    hit = index.scene.cast_rays(library().core.Tensor(rays), nthreads=2)
    distance = hit['t_hit'].numpy().astype(np.float64)
    face_ids = hit['primitive_ids'].numpy().astype(np.int64)
    normals = hit['primitive_normals'].numpy().astype(np.float64)
    normal_length = np.linalg.norm(normals, axis=1)
    normals /= np.maximum(normal_length[:, None], 1e-30)
    valid = np.isfinite(distance) & (normal_length > 0)
    valid &= (normals @ direction) >= minimum_facing_cosine
    world_hit = np.full_like(points, np.nan)
    world_hit[valid] = (origins[valid] - distance[valid, None] * direction) * index.scale + index.origin
    gap = np.full(len(points), np.nan)
    gap[valid] = np.einsum('ij,ij->i', points[valid] - world_hit[valid], normals[valid])
    return dict(valid=valid, points=world_hit, normals=normals,
                gap_m=gap, face_ids=face_ids)
