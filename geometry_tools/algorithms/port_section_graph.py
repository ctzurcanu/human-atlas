"""Place complete source tube sections between declared source mouth poses.

Positive source-length graph weights prevent pulling individual vessel walls
into a collapsed surface-harmonic collar. Mouth positions and proper frames
are explicit; branch chains transport their section frames coherently.
"""
import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

from .profile_graph_transport import transport_linked_profiles


def shortest_rotation(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if min(np.linalg.norm(a), np.linalg.norm(b)) < 1e-12:
        raise ValueError('A collapsed station span cannot determine a section frame')
    a, b = a/np.linalg.norm(a), b/np.linalg.norm(b)
    c = float(np.clip(a@b, -1, 1))
    if c < -1+1e-10:
        axis = np.eye(3)[np.argmin(abs(a))]
        axis -= a*(axis@a); axis /= np.linalg.norm(axis)
        return 2*np.outer(axis, axis)-np.eye(3)
    v = np.cross(a, b)
    K = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
    return np.eye(3)+K+K@K/(1+c)


def pose_section_graph(graph, fixed_positions, fixed_rotations):
    p = graph['centres']
    adjacency = graph['graph']
    a, b = adjacency.nonzero()
    length = np.linalg.norm(p[a]-p[b], axis=1)
    if np.any(length <= 1e-12):
        raise ValueError('Source section graph has a collapsed span')
    weights = coo_matrix((1/length, (a, b)), shape=adjacency.shape).tocsr()
    L = diags(np.asarray(weights.sum(1)).ravel())-weights
    q = p.copy()
    fixed = set(fixed_positions)
    for k, point in fixed_positions.items():
        q[k] = point
    _, labels = connected_components(weights, directed=False)
    for ci in np.unique(labels):
        component = np.flatnonzero(labels == ci)
        if not fixed.intersection(component):
            fixed.update(map(int, component))
    protected = np.array(sorted(fixed), int)
    free = np.setdiff1d(np.arange(len(p)), protected)
    if len(free):
        q[free] = spsolve(L[free][:, free].tocsc(), -L[free][:, protected]@q[protected])
    R = np.repeat(np.eye(3)[None], len(p), axis=0)
    for node in range(len(p)):
        if node in fixed_rotations:
            R[node] = fixed_rotations[node]
            continue
        neighbours = adjacency.indices[adjacency.indptr[node]:adjacency.indptr[node+1]]
        if len(neighbours) < 1:
            continue
        dx, dy = p[neighbours]-p[node], q[neighbours]-q[node]
        if len(neighbours) <= 2:
            R[node] = shortest_rotation(dx[0] if len(dx) == 1 else dx[-1]-dx[0],
                                        dy[0] if len(dy) == 1 else dy[-1]-dy[0])
        else:
            u, _, vt = np.linalg.svd(dy.T@dx)
            sign = np.ones(3); sign[-1] = np.linalg.det(u@vt)
            R[node] = (u*sign[None])@vt
    links = np.array([(i, j) for i, j in zip(a, b) if i < j], int).reshape(-1, 2)
    R, proof = transport_linked_profiles(p, q, R, links, protected=protected)
    for k, point in fixed_positions.items():
        if not np.array_equal(q[k], point):
            raise ValueError('A declared source mouth centre moved')
    return q, R, dict(sourceNodes=len(p), fixedNodes=len(protected),
                      frames=proof, anatomicalAcceptance=False,
                      positionsAreSourceLengthWeighted=True)
