"""Proper source-section alignment from literal adjacent exporter quads."""
from collections import Counter
import numpy as np


def proper_profile_rotation(source, target):
    source, target = np.asarray(source, float), np.asarray(target, float)
    if source.shape != target.shape or source.ndim != 2 or source.shape[1] != 3 or len(source) < 3:
        raise ValueError('Corresponding source and target section controls required')
    p, q = source-source.mean(0), target-target.mean(0)
    if not np.isfinite(p).all() or not np.isfinite(q).all():
        raise ValueError('Finite section controls required')
    u, singular, vt = np.linalg.svd(q.T@p)
    if singular[1] <= 1e-25:
        raise ValueError('Source section cannot determine a proper plane orientation')
    sign = np.eye(3)
    sign[-1, -1] = np.linalg.det(u@vt)
    rotation = u@sign@vt
    assert np.linalg.det(rotation) > 1-1e-12
    return rotation


def align_source_ring(mesh, topology, ring_id, *, maximum_rotation_degrees=60):
    inverse = np.asarray(topology['raw_to_physical'])
    physical = np.asarray(topology['physical_points'])
    np.testing.assert_array_equal(physical[inverse], mesh.vertices)
    ring = np.asarray(topology['rings'][ring_id], dtype=int)
    selected = set(ring)
    strips = [s for s in topology['strips'] if ring_id in s['rings']]
    if not 1 <= len(strips) <= 2:
        raise ValueError('A regular section needs one or two verified adjacent strips')
    source = physical[ring]
    energy = np.linalg.norm(source-source.mean(0))
    targets, neighbors = [], []
    faces = inverse[mesh.faces]
    for strip in strips:
        other = next(i for i in strip['rings'] if i != ring_id)
        other_ring = set(topology['rings'][other])
        mapping = {}
        source_faces = set(strip['source_face_ids'])
        for i in sorted(source_faces):
            if i % 2:
                continue
            assert i+1 in source_faces
            counts = Counter(tuple(sorted((int(a), int(b)))) for tri in faces[i:i+2]
                             for a, b in zip(tri, np.roll(tri, -1)))
            for (a, b), count in counts.items():
                if count != 1 or not ((a in selected and b in other_ring) or (b in selected and a in other_ring)):
                    continue
                a, b = (a, b) if a in selected else (b, a)
                assert a not in mapping or mapping[a] == b
                mapping[a] = b
        if set(mapping) != selected or len(set(mapping.values())) != len(ring):
            raise ValueError('Original quads do not give a bijective section correspondence')
        neighbor = physical[[mapping[int(i)] for i in ring]]
        centered = neighbor-neighbor.mean(0)
        if np.linalg.norm(centered) <= 1e-20:
            raise ValueError('Adjacent source section is collapsed')
        targets.append(centered*energy/np.linalg.norm(centered))
        neighbors.append(dict(ringId=other, correspondingPhysicalIds=[mapping[int(i)] for i in ring],
                              sourceFaceIds=strip['source_face_ids']))
    target = np.mean(targets, axis=0)
    rotation = proper_profile_rotation(source, target)
    angle = float(np.degrees(np.arccos(np.clip((np.trace(rotation)-1)/2, -1, 1))))
    if angle > maximum_rotation_degrees:
        raise ValueError('Section orientation correction exceeds the explicit rotation bound')
    ids = np.flatnonzero(np.isin(inverse, ring))
    center = source.mean(0)
    points = mesh.vertices.copy()
    points[ids] = center+(points[ids]-center)@rotation.T
    return points, dict(ringId=ring_id, rawIds=ids.tolist(), neighbors=neighbors,
                        rotation=rotation.tolist(), rotationDegrees=angle, center_m=center.tolist(),
                        profileShapeAreaAndCenterRetained=True, exactSourceAliasesIncluded=True,
                        sourceSectionsAssumedCircular=False, anatomicalAcceptance=False)
