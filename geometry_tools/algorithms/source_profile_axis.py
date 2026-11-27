"""Read-only tube-axis guides from explicit source annuli and mesh incidence.

Contract each actual profile and each remaining source core component. These
are numerical course guides, not newly inferred anatomical branch identities.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra

from .source_tube_topology import recover_regular_source_annuli
from .tube_continuity import boundary_loops


def complete_source_sections(graph,physical_ids):
    """Include actual open mouths in disjoint rigid source control closures.

    A missing quad partner at an open source port does not make its six
    perimeter controls independent. Partly shared numerical annuli and a real
    boundary ring co-move as one closure; no source point is welded or removed.
    """
    physical_ids=np.asarray(physical_ids,int);groups=[]
    for raw in graph['profiles']+graph['boundary_raw']:
        members=set(map(int,physical_ids[raw]))
        if len(members)!=len(raw):raise ValueError('A source ring repeats a physical control')
        touching=[i for i,g in enumerate(groups)if g.intersection(members)]
        if not touching:groups.append(members);continue
        first=touching[0]
        for i in reversed(touching[1:]):members.update(groups.pop(i))
        groups[first].update(members)
    representative={int(k):int(np.flatnonzero(physical_ids==k)[0])for k in np.unique(physical_ids)}
    return [np.array([representative[k]for k in sorted(g)],int)for g in groups]


def source_profile_graph(vertices, faces, *, triangle_pairing, ring_size=6, extend_annuli=False):
    v, f = np.asarray(vertices, float), np.asarray(faces, int)
    rec = recover_regular_source_annuli(v, f, triangle_pairing, ring_size=ring_size)
    p, inv = rec['physical_points'], rec['raw_to_physical']
    profiles = [np.array(r, int) for r in rec['rings']]
    extended = None
    if extend_annuli:
        if ring_size != 6:
            raise ValueError('Source propagation currently requires six-control seeds')
        from .source_annular_propagation import propagate_source_annuli
        extended = propagate_source_annuli(v, f, [rec['representative_raw_ids'][r] for r in profiles])
        profiles = [rec['raw_to_physical'][r] for r in extended['profiles']]
    owner = np.full(len(p), -1, int)
    centres = []
    for i, ids in enumerate(profiles):
        if np.any(owner[ids] >= 0):
            raise ValueError('Source profiles overlap in physical controls')
        owner[ids] = i
        centres.append(p[ids].mean(0))
    ef = inv[f]
    edges = np.unique(np.sort(ef[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1), axis=0)
    free = np.flatnonzero(owner < 0)
    if len(free):
        local = np.full(len(p), -1, int)
        local[free] = np.arange(len(free))
        live = edges[np.all(owner[edges] < 0, axis=1)]
        graph = coo_matrix((np.ones(2*len(live)),
                            (np.r_[local[live[:, 0]], local[live[:, 1]]],
                             np.r_[local[live[:, 1]], local[live[:, 0]]])),
                           shape=(len(free), len(free))).tocsr()
        count, label = connected_components(graph, directed=False)
        for i in range(count):
            ids = free[label == i]
            owner[ids] = len(centres)
            centres.append(p[ids].mean(0))
    centres = np.asarray(centres)
    links = np.unique(np.sort(owner[edges], axis=1), axis=0)
    links = links[links[:, 0] != links[:, 1]]
    boundary_nodes, boundary_raw = [], []
    for loop in boundary_loops(ef):
        point = p[loop].mean(0)
        incident = np.unique(owner[loop])
        if len(incident) == 1 and np.linalg.norm(centres[incident[0]]-point) < 1e-12:
            node = int(incident[0])
        else:
            node = len(centres)
            centres = np.r_[centres, point[None]]
            links = np.r_[links, np.c_[np.full(len(incident), node), incident]]
        boundary_nodes.append(node)
        boundary_raw.append(rec['representative_raw_ids'][loop])
    length = np.linalg.norm(centres[links[:, 0]]-centres[links[:, 1]], axis=1)
    # A zero centroid separation is not a new anatomical attachment. Its
    # existing source incidence receives only a numerical traversal cost.
    length = np.maximum(length, 1e-12)
    graph = coo_matrix((np.r_[length, length],
                        (np.r_[links[:, 0], links[:, 1]], np.r_[links[:, 1], links[:, 0]])),
                       shape=(len(centres), len(centres))).tocsr()
    return dict(centres=centres, graph=graph, raw_node=owner[inv],
                profiles=[rec['representative_raw_ids'][r] for r in profiles],
                boundary_nodes=boundary_nodes, boundary_raw=boundary_raw,
                annularPropagation=extended,
                sourceIncidenceOnly=True)


def primary_axis(guide, face_vertex_ids):
    nodes = np.unique(guide['raw_node'][np.unique(face_vertex_ids)])
    graph, p = guide['graph'], guide['centres']
    raw = np.unique(face_vertex_ids)
    ends = np.array([node for node, ids in zip(guide['boundary_nodes'], guide['boundary_raw'])
                     if np.isin(ids, raw).all()], int)
    if len(ends) < 2:
        raise ValueError('Primary course needs at least two source terminal nodes')
    root = int(ends[np.argmax(p[ends, 1])])
    distance, previous = dijkstra(graph, directed=False, indices=root, return_predecessors=True)
    terminal = ends[np.isfinite(distance[ends])]
    stop = int(terminal[np.argmin(p[terminal, 1])])
    route = [stop]
    while route[-1] != root:
        parent = int(previous[route[-1]])
        if parent < 0 or len(route) > len(nodes)+len(ends):
            raise ValueError('Source course did not replay its actual connected component')
        route.append(parent)
    axis = p[route[::-1]]
    if len(axis) < 3 or np.any(np.linalg.norm(np.diff(axis, axis=0), axis=1) <= 1e-10):
        raise ValueError('Degenerate source course guide')
    return axis, dict(sourceNodes=route[::-1], proximalNode=root, distalNode=stop,
                      courseLength_m=float(distance[stop]), anatomicalAcceptance=False)
