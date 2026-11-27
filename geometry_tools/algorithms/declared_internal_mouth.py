"""Close a declared same-carrier source seam, without caps or face removal."""
import numpy as np
from .tube_continuity import boundary_loops, continuity_audit
from .surface_normals import transport_vertex_words
from ..model import Mesh


def weld_internal_mouth(mesh, first_ids, second_ids, *, maximum_gap_m=.00003, maximum_relative_area_change=.025):
    first_ids, second_ids = [np.asarray(x, dtype=int) for x in [first_ids, second_ids]]
    if first_ids.shape != second_ids.shape or first_ids.ndim != 1 or len(first_ids) < 3:
        raise ValueError('Two complete corresponding source rims are required')
    physical, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
    rings = boundary_loops(inverse[mesh.faces])
    actual = []
    for ids in [first_ids, second_ids]:
        matched = [x for x in rings if set(x) == set(inverse[ids])]
        if len(matched) != 1 or len(set(inverse[ids])) != len(ids):
            raise ValueError('Declared rim must be one complete external source boundary')
        actual.append(physical[matched[0]])
    if set(inverse[first_ids]).intersection(inverse[second_ids]):
        raise ValueError('The two source rims must initially be distinct')
    area_vector = lambda p: np.cross(p-p.mean(0), np.roll(p, -1, axis=0)-p.mean(0)).sum(0)/2
    n0, n1 = [area_vector(x) for x in actual]
    cosine = float(n0@n1/(np.linalg.norm(n0)*np.linalg.norm(n1)))
    if cosine >= -.99:
        raise ValueError('Actual source boundary windings must be opposed')
    p, q = mesh.vertices[first_ids], mesh.vertices[second_ids]
    gap = float(np.linalg.norm(p-q, axis=1).max())
    if gap > maximum_gap_m:
        raise ValueError('Source seam exceeds its explicit matching tolerance')
    midpoint = (p+q)/2
    new_area = np.linalg.norm(area_vector(midpoint))
    changes = [float(abs(new_area/np.linalg.norm(n)-1)) for n in [n0, n1]]
    if max(changes) > maximum_relative_area_change:
        raise ValueError('Common section area exceeds the declared source tolerance')
    points = mesh.vertices.copy()
    moved = set()
    for a, b, common in zip(p, q, midpoint):
        ids = np.flatnonzero(np.all(mesh.vertices == a, axis=1) | np.all(mesh.vertices == b, axis=1))
        points[ids] = common
        moved.update(ids.tolist())
    words, normal_proof = transport_vertex_words(mesh.vertices, points, mesh.faces, mesh.attributes['normals'], maximum_facet_rotation_degrees=10)
    # One effective incident-face normal at each newly shared seam position;
    # source normalSigns and UV discontinuities remain literal.
    triangles = points[mesh.faces]
    facet = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
    for common in midpoint:
        ids = np.flatnonzero(np.all(points == common, axis=1))
        face_ids = np.flatnonzero(np.any(np.isin(mesh.faces, ids), axis=1))
        direction = facet[face_ids].sum(0)
        direction /= np.linalg.norm(direction)
        if not np.all(facet[face_ids]@direction > 0):
            raise ValueError('Joined seam has no outward area-weighted normal cone')
        signs = mesh.attributes['normalSigns'][ids]
        words[ids] = np.rint(direction[None]*32767*signs[:, None]).astype('<i2')
    candidate = Mesh(points, mesh.faces, mesh.frame, dict(mesh.attributes, normals=words), mesh.source)
    before = continuity_audit(mesh.vertices, mesh.faces)
    expected_euler = sum(before['component_euler'])
    # Joining distinct components reduces their number by one. Closing two
    # ends of an already connected vascular network preserves its components.
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    edges = inverse[mesh.faces][:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
    graph = coo_matrix((np.ones(2*len(edges)), (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])), shape=(len(physical), len(physical))).tocsr()
    _, labels = connected_components(graph, directed=False)
    merge_components = int(labels[inverse[first_ids[0]]] != labels[inverse[second_ids[0]]])
    after = continuity_audit(points, mesh.faces, expected_mouths=before['mouths']-2,
                              expected_components=before['components']-merge_components, expected_euler=expected_euler)
    if not after['topological_continuity_pass']:
        raise ValueError('Joined source carrier fails its complete topology contract')
    return candidate, dict(complete=True, firstRing=first_ids.tolist(), secondRing=second_ids.tolist(),
             actualSourceNormalCosine=cosine, maximumSourceGap_m=gap, commonSectionArea_m2=float(new_area),
             relativeSectionAreaChanges=changes, maximumMotion_m=float(np.linalg.norm(points-mesh.vertices, axis=1).max()),
             movedRawIds=sorted(moved), topologyBefore=before, topologyAfter=after,
             allOriginalFacesRetained=True, lumenDisksCreated=0, sourceNormalTransport=normal_proof,
             commonEffectiveSeamNormals=True, anatomicalAttachmentEstablished=False, anatomicalAcceptance=False)
