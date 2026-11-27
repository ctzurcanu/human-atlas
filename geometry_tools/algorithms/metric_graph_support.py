"""Distance-limited support on an existing source dependency graph.

No spatial neighbors or anatomical relationships are invented. Metric length
is measured between current group centres along the caller's existing edges.
The support radius is a proposal scope, never a deformation certificate.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import dijkstra


def metric_support(graph, centres, seeds, radius_m):
    centres=np.asarray(centres,float);seeds=np.unique(np.asarray(seeds,int))
    if centres.ndim!=2 or centres.shape[1]!=3 or graph.shape!=(len(centres),len(centres)):
        raise ValueError('Corresponding source graph and current centres required')
    if not np.isfinite(centres).all() or not np.isfinite(radius_m) or radius_m<=0:
        raise ValueError('Finite current centres and positive metric radius required')
    if not len(seeds) or seeds.min()<0 or seeds.max()>=len(centres):
        raise ValueError('Explicit existing source seed groups required')
    edges=graph.tocoo();keep=edges.row!=edges.col;rows=edges.row[keep];cols=edges.col[keep]
    length=np.linalg.norm(centres[rows]-centres[cols],axis=1)
    # A 1 nm lower bound gives coincident source dependency centres finite
    # positive graph weight; it does not merge points or create new edges.
    length=np.maximum(length,1e-9)
    metric=sparse.coo_matrix((length,(rows,cols)),shape=graph.shape).tocsr()
    distance=dijkstra(metric,directed=False,indices=seeds,min_only=True,limit=radius_m)
    selected=np.flatnonzero(distance<=radius_m)
    return selected,dict(complete=True,radius_m=float(radius_m),seedCount=len(seeds),
                         selectedCount=len(selected),existingGraphEdgesOnly=True,
                         spatialEdgesInferred=False,minimumEdgeWeight_m=1e-9,
                         distanceIsNotAnatomicalAcceptance=True)
