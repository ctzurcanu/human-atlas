"""Resolve duplicate source strips using an existing coherent normal field.

Only exactly coincident, same-winding facets can lose a redundant draw copy.
All non-normal corner fields must agree. Connected duplicate strips choose
one complete existing copy by its normal discontinuity against adjacent
nonduplicate source corners; no synthetic normal, position or weld is made.
"""
from collections import defaultdict
import numpy as np
from .tube_continuity import continuity_audit


def resolve_source_normal_duplicates(mesh):
    vertices, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
    physical = inverse[mesh.faces]
    buckets = defaultdict(list)
    phases = {}
    for i, face in enumerate(physical):
        phase = min(range(3), key=lambda k: tuple(np.roll(face, -k)))
        phases[i] = np.roll(np.arange(3), -phase)
        buckets[tuple(face[phases[i]])].append(i)
    groups = [ids for ids in buckets.values() if len(ids) > 1]
    if not groups:
        raise ValueError('No same-winding duplicate source facets')
    if any(len(ids) != 2 for ids in groups):
        raise ValueError('Only explicitly paired source copies supported')
    all_ids = set(sum(groups, []))
    normals = mesh.attributes['normals'].astype(float) * mesh.attributes['normalSigns'][:, None]
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    for ids in groups:
        a, b = ids
        first, second = mesh.faces[a, phases[a]], mesh.faces[b, phases[b]]
        np.testing.assert_array_equal(mesh.vertices[first], mesh.vertices[second])
        for key in mesh.attributes:
            if key not in {'normals', 'normalSigns'}:
                np.testing.assert_array_equal(mesh.attributes[key][first], mesh.attributes[key][second])
    # A strip is a connected set of duplicate geometric triangles.
    adjacent = [set() for _ in groups]
    by_edge = defaultdict(list)
    for gi, ids in enumerate(groups):
        face = physical[ids[0]]
        for a, b in zip(face, np.roll(face, -1)):
            by_edge[tuple(sorted((int(a), int(b))))].append(gi)
    for owners in by_edge.values():
        for a in owners:
            adjacent[a].update(set(owners) - {a})
    neighbors = defaultdict(list)
    for fi, face in enumerate(physical):
        if fi not in all_ids:
            for k, vertex in enumerate(face):
                neighbors[int(vertex)].append(normals[mesh.faces[fi, k]])
    remaining = set(range(len(groups)))
    removed = []
    strips = []
    while remaining:
        todo = [min(remaining)]
        component = set()
        while todo:
            i = todo.pop()
            if i not in component:
                component.add(i)
                todo += list(adjacent[i] - component)
        remaining -= component
        pairs = [sorted(groups[i]) for i in sorted(component)]
        # The ordinal copies must each form the same contiguous source strip.
        copies = [np.asarray([pair[k] for pair in pairs]) for k in range(2)]
        for copy in copies:
            audit = continuity_audit(vertices, physical[copy], expected_mouths=2,
                                     expected_components=1, expected_euler=0)
            if not audit['topological_continuity_pass']:
                raise ValueError('Duplicate copies must each be one complete source annulus')
        scores = []
        for copy in copies:
            values = []
            for fi in copy:
                for corner, vertex in enumerate(physical[fi]):
                    n = normals[mesh.faces[fi, corner]]
                    values += [float(1 - np.clip(n @ other, -1, 1)) for other in neighbors[int(vertex)]]
            if not values:
                raise ValueError('Duplicate strip has no neighboring source normal field')
            scores.append(float(np.mean(values)))
        winner = int(np.argmin(scores))
        retained, discarded = copies[winner], copies[1 - winner]
        removed.extend(discarded.tolist())
        correspondence = []
        for a, b in zip(retained, discarded):
            order = [int(np.flatnonzero(physical[b] == vertex)[0]) for vertex in physical[a]]
            correspondence.append(dict(retainedSourceFace=int(a), removedSourceFace=int(b), removedCornerOrder=order))
        strips.append(dict(retainedSourceFaceIds=retained.tolist(), removedSourceFaceIds=discarded.tolist(),
                           sourceNormalDiscontinuityScores=scores, selectedCopy=winner,
                           correspondence=correspondence))
    keep = np.setdiff1d(np.arange(len(mesh.faces)), removed)
    return keep, dict(complete=True, strips=strips, removedSourceFaceIds=sorted(removed),
                      sourceOrientedSurfaceSupportUnchanged=True, originalControlsAndFieldsUnchanged=True,
                      normalChoice='One existing source strip, minimum discontinuity to nonduplicate source corners',
                      anatomicalAcceptance=False)
