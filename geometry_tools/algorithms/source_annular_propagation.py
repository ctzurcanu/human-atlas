"""Extend certified source section cycles through literal oriented annuli.

No exporter quad ordering is presumed. Each proposed next ring must retain
its seed's distinct control count N, a simple source-edge cycle and exactly 2N
source triangles forming an oppositely wound annulus. Forks/ambiguous strips
remain explicit unrecovered source faces. No source geometry moves.
"""
from collections import defaultdict, deque

import numpy as np

from .source_tube_topology import _annulus, _cycles


def propagate_source_annuli(vertices, faces, seed_raw_rings):
    v, raw = np.asarray(vertices, float), np.asarray(faces, int)
    points, representative, inv = np.unique(v, axis=0, return_index=True, return_inverse=True)
    f = inv[raw]
    directed = defaultdict(list)
    edge_faces = defaultdict(list)
    for i, tri in enumerate(f):
        for a, b, c in [(tri[0], tri[1], tri[2]), (tri[1], tri[2], tri[0]), (tri[2], tri[0], tri[1])]:
            directed[int(a), int(b)].append((i, int(c)))
            edge_faces[tuple(sorted((int(a), int(b))))].append(i)
    rings = []
    index = {}
    owner = {}
    for seed in seed_raw_rings:
        ring = list(map(int, inv[np.asarray(seed, int)]))
        key = frozenset(ring)
        if len(key) != len(ring) or len(ring) < 3:
            raise ValueError('Propagation requires at least three distinct literal seed controls')
        if key in index:
            continue
        if any(p in owner for p in ring):
            raise ValueError('Source seed sections overlap')
        index[key] = len(rings)
        for p in ring:
            owner[p] = len(rings)
        rings.append(ring)
    seed_count = len(rings)
    todo = deque(range(seed_count))
    strips = {}
    face_owner = {}
    rejected = 0
    while todo:
        first = todo.popleft()
        ring = rings[first]
        count = len(ring)
        for reverse in [False, True]:
            sequence = ring[::-1] if reverse else ring
            starting = []
            for a, b in zip(sequence, np.roll(sequence, -1)):
                alternatives = [(i, c) for i, c in directed[a, int(b)] if c not in ring]
                if len(alternatives) != 1:
                    starting = []
                    break
                starting.append(alternatives[0][0])
            if len(starting) != count:
                continue
            # Complete the one-sided vertex fan. Alternating source diagonals
            # may give duplicate opposite controls in the first N faces.
            # Flood only across non-ring edges and only through faces still
            # incident to this original ring; never traverse a second layer.
            patch = set()
            pending = list(starting)
            ring_set = set(ring)
            while pending and len(patch) <= 2*count:
                i = pending.pop()
                if i in patch:
                    continue
                patch.add(i)
                tri = f[i]
                for a, b in zip(tri, np.roll(tri, -1)):
                    if a in ring_set and b in ring_set:
                        continue
                    for other in edge_faces[tuple(sorted((int(a), int(b))))]:
                        if other not in patch and ring_set.intersection(map(int, f[other])):
                            pending.append(other)
            faces_here = sorted(patch)
            if len(faces_here) != 2*count:
                continue
            outside = set(map(int, f[faces_here].ravel()))-ring_set
            if len(outside) != count:
                continue
            edges = set()
            for i in faces_here:
                tri = f[i]
                for a, b in zip(tri, np.roll(tri, -1)):
                    if a in outside and b in outside:
                        edges.add(tuple(sorted((int(a), int(b)))))
            cycles = _cycles(edges)
            if len(cycles) != 1 or set(cycles[0]) != set(outside):
                continue
            second_ring = cycles[0]
            if not _annulus(f[faces_here], ring, second_ring):
                continue
            key = frozenset(second_ring)
            second = index.get(key)
            if second is None and any(p in owner for p in second_ring):
                rejected += 1
                continue
            pair = tuple(sorted((first, second))) if second is not None else None
            if any(i in face_owner and face_owner[i] != pair for i in faces_here):
                rejected += 1
                continue
            if second is None:
                second = len(rings)
                rings.append(second_ring)
                index[key] = second
                for p in second_ring:
                    owner[p] = second
                todo.append(second)
            pair = tuple(sorted((first, second)))
            strips[pair] = faces_here
            for i in faces_here:
                face_owner[i] = pair
    return dict(profiles=[representative[r] for r in rings],
                strips=[dict(rings=list(k), source_face_ids=ids) for k, ids in sorted(strips.items())],
                seedProfiles=seed_count, propagatedProfiles=len(rings)-seed_count,
                ambiguousStripsRejected=rejected, sourceFacesOwned=len(face_owner),
                unassignedSourceFaceIds=sorted(set(range(len(f)))-set(face_owner)),
                sourceCoordinatesUnchanged=True, exporterPairingAssumed=False,
                sourceSectionControlCounts=sorted({len(r)for r in rings}),
                anatomicalAcceptance=False)
