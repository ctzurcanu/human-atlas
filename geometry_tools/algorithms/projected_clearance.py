"""Whole-triangle directional clearance against a fixed surface projection.

This is a directional packing constraint, not containment or anatomical
registration. A subsequent complete contact check is required, especially for
the zero-area projected faces. Distances and coordinates are in metres.
"""
import numpy as np
import shapely
from shapely.strtree import STRtree


def directional_clearance(moving, fixed, direction=(0., 0., 1.), *, max_pairs=250000):
    direction = np.asarray(direction, float)
    direction /= np.linalg.norm(direction)
    axis = np.eye(3)[np.argmin(abs(direction))]
    first = np.cross(direction, axis)
    first /= np.linalg.norm(first)
    basis = np.column_stack([first, np.cross(direction, first), direction])
    projected = [np.asarray(m.triangles) @ basis for m in (moving, fixed)]
    triangles, ids, planes, polygons, skipped = [], [], [], [], []
    for q in projected:
        n = np.cross(q[:, 1]-q[:, 0], q[:, 2]-q[:, 0])
        eligible = abs(n[:, 2]) > 1e-20
        rows = np.flatnonzero(eligible)
        ids.append(rows)
        triangles.append(q[rows])
        skipped.append(int((~eligible).sum()))
        nn = n[rows]
        planes.append(np.column_stack([-nn[:, :2]/nn[:, 2, None],
            np.einsum('ij,ij->i', nn, q[rows, 0])/nn[:, 2]]))
        polygons.append(shapely.polygons(q[rows, :, :2]))
    pairs = STRtree(polygons[1]).query(polygons[0], predicate='intersects')
    if pairs.shape[1] > max_pairs:
        raise ValueError('Directional clearance projection pair budget exceeded')
    if not pairs.shape[1]:
        return dict(projection_overlaps=False, candidate_pairs=0,
                    zero_area_projection_faces=skipped)
    intersections = shapely.intersection(polygons[0][pairs[0]], polygons[1][pairs[1]])
    xy, owners = shapely.get_coordinates(intersections, return_index=True)
    if not len(xy):
        raise ValueError('Projected overlap has no finite witnesses')
    delta = planes[0][pairs[0, owners]] - planes[1][pairs[1, owners]]
    excess = np.einsum('ij,ij->i', np.column_stack([xy, np.ones(len(xy))]), delta)
    worst = int(excess.argmax())
    pair = pairs[:, owners[worst]]
    witnesses = []
    for k in range(2):
        z = np.r_[xy[worst], 1.] @ planes[k][pair[k]]
        witnesses.append((np.r_[xy[worst], z] @ basis.T).tolist())
    return dict(projection_overlaps=True, candidate_pairs=int(pairs.shape[1]),
        intersection_vertices=len(xy), maximum_forward_excess_m=float(excess[worst]),
        direction=direction.tolist(), face_ids=[int(ids[k][pair[k]]) for k in range(2)],
        moving_witness_m=witnesses[0], fixed_witness_m=witnesses[1],
        zero_area_projection_faces=skipped,
        scope='All intersecting nondegenerate triangle projections; separate 3D contact check required')
