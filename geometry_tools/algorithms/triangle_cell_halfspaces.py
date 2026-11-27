"""Separating halfspaces for bounded, fixed-Y rigid-cell translation.

A moving triangle stays disjoint from a fixed triangle when all its corners
remain on the certified separating side. This is a conservative convex
constraint; intersecting baseline pairs are returned as protected sites.
"""
import numpy as np
from geometry_tools.distance import triangle_pairs


def moving_triangle_cell_halfspaces(first, first_cells, second, second_cells,
                                   *, margin_m=.00005):
    """Joint source-side constraints, allowing both actual triangles to move.

    All nine corner differences on one certified source separating normal must
    retain a positive gap. Fixed-Y translations make these inequalities exact.
    This is a local proposal; actual all-neighbour crossings still need audit.
    """
    first, second = np.asarray(first, float), np.asarray(second, float)
    ca, cb = np.asarray(first_cells, int), np.asarray(second_cells, int)
    if (first.shape != second.shape or first.ndim != 3 or first.shape[1:] != (3, 3)
            or ca.shape != first.shape[:2] or cb.shape != ca.shape
            or not np.isfinite(np.r_[first.ravel(), second.ravel()]).all()
            or margin_m < 0):
        raise ValueError('Corresponding finite triangles and source corner cells required')
    answer = triangle_pairs(first, second)
    rows, unresolved, certificates = [], [], []
    for i, distance in enumerate(answer['distance_m']):
        if distance <= 1e-10 or answer['degenerate'][i]:
            unresolved.append(i)
            continue
        n = (answer['first'][i]-answer['second'][i])/distance
        gap = float(np.min(first[i]@n)-np.max(second[i]@n))
        if gap <= 1e-10:
            unresolved.append(i)
            continue
        margin = min(margin_m, max(0., gap-1e-8))
        for j in range(3):
            for k in range(3):
                coefficients = {}
                for cell, sign in [(int(ca[i, j]), 1.), (int(cb[i, k]), -1.)]:
                    coefficients[cell] = coefficients.get(cell, np.zeros(2))+sign*n[[0, 2]]
                rows.append(dict(coefficients={str(key): value.tolist() for key, value in coefficients.items()},
                                 minimum_m=margin-float((first[i, j]-second[i, k])@n)))
        certificates.append(dict(pair=i, supportGap_m=gap, margin_m=margin, normal=n.tolist()))
    return rows, unresolved, certificates


def triangle_cell_halfspaces(first, first_cells, second, *, margin_m=.00075,
                             maximum_motion_m=.008):
    first,second=np.asarray(first,float),np.asarray(second,float);cells=np.asarray(first_cells,int)
    if first.shape!=second.shape or first.ndim!=3 or first.shape[1:]!=(3,3)or cells.shape!=first.shape[:2]or not np.isfinite(np.r_[first.ravel(),second.ravel()]).all():
        raise ValueError('Corresponding actual triangle pairs and corner source cells required')
    if margin_m<0 or not 0<maximum_motion_m<=.010:raise ValueError('Explicit source clearance and motion limits required')
    answer=triangle_pairs(first,second);rows=[];protected=set();evidence=[]
    for i,distance in enumerate(answer['distance_m']):
        if distance>maximum_motion_m+1e-7:continue
        if distance<=1e-10 or answer['degenerate'][i]:
            protected.update(map(int,cells[i]));evidence.append(dict(pair=i,protected=True,distance_m=float(distance)));continue
        n=(answer['first'][i]-answer['second'][i])/distance;gap=float(np.min(first[i]@n)-np.max(second[i]@n))
        if gap<=1e-10:
            protected.update(map(int,cells[i]));evidence.append(dict(pair=i,protected=True,distance_m=float(distance),supportGap_m=gap));continue
        required=min(margin_m,max(0.,gap-1e-8));normal=np.round(n[[0,2]],12);rounding=maximum_motion_m*float(np.linalg.norm(normal-n[[0,2]]));bound=float(np.max(second[i]@n))+required
        for cell in np.unique(cells[i]):
            raw=np.flatnonzero(cells[i]==cell);minimum=float(np.max(bound-first[i,raw]@n))+rounding
            if np.linalg.norm(normal)>1e-14:rows.append(dict(cell=int(cell),normalXZ=normal.tolist(),minimum_m=minimum))
        evidence.append(dict(pair=i,protected=False,distance_m=float(distance),supportGap_m=gap,margin_m=required,normal=n.tolist()))
    return rows,protected,evidence
