"""Bound a local mesh step using literal changed-face contact counts.

Every pair that can change has an incident moving face. Unchanged pairs cancel
from the comparison. Existing crossings remain unaccepted; a passing step is
only a non-regression check, not a tissue, lumen or anatomical certificate.
"""
from collections import Counter

import numpy as np

from geometry_tools.contact import pair_contacts
from geometry_tools.model import Mesh


def _pairs(first, second):
    arrays, proof = pair_contacts(first, second, 0., max_tests=250000, max_pairs=100000)
    if not proof['complete_within_pair_budget']:
        raise RuntimeError('Incomplete local triangle transition; no passing result')
    return arrays['pairs'][np.isin(arrays['relation'], ['crossing', 'coplanar-overlap'])]


def check_transition(before, after, face_members, fixed=()):
    """Compare all changed pairs, grouped by physical mesh owners.

    ``fixed`` contains (unique explicit name, Mesh) pairs in the same frame.
    Inter-member contacts are counted once; source face order is preserved.
    A pair may migrate across facets within an already crossed owner pair,
    but its count cannot rise. A previously clear owner pair must stay clear.
    """
    if before.frame != after.frame or not np.array_equal(before.faces, after.faces):
        raise ValueError('Transition requires one frame and literal corresponding faces')
    members = np.asarray(face_members)
    if members.shape != (len(before.faces),):
        raise ValueError('Every source face requires its original mesh owner')
    moved = np.any(before.vertices != after.vertices, axis=1)
    changed = np.any(moved[before.faces], axis=1)
    ids = np.flatnonzero(changed)
    if not len(ids):
        return dict(complete=True, passed=True, changedFaces=0, rows=[])
    snapshots = [before, after]
    primary = [Mesh(m.vertices, m.faces[ids], m.frame) for m in snapshots]
    triangles = np.concatenate([m.triangles for m in primary])
    low, high = triangles.min((0, 1))-1e-9, triangles.max((0, 1))+1e-9
    candidate_triangles = [m.triangles for m in snapshots]
    keep = np.zeros(len(before.faces), bool)
    for t in candidate_triangles:
        keep |= np.all(t.max(1) >= low, axis=1) & np.all(t.min(1) <= high, axis=1)
    peers = np.flatnonzero(keep)
    counts = []
    for p, m in zip(primary, snapshots):
        pairs = _pairs(p, Mesh(m.vertices, m.faces[peers], m.frame))
        first, second = ids[pairs[:, 0]], peers[pairs[:, 1]]
        # Ignore the same face; count pairs with two moving faces only once.
        take = (first != second) & (~changed[second] | (first < second))
        counts.append(Counter(tuple(sorted((str(members[a]), str(members[b]))))
                              for a, b in zip(first[take], second[take])))
    fixed_names = set()
    for name, mesh in fixed:
        if name in fixed_names or mesh.frame != before.frame:
            raise ValueError('Fixed meshes require unique names and the same frame')
        fixed_names.add(name)
        t = mesh.triangles
        keep = np.all(t.max(1) >= low, axis=1) & np.all(t.min(1) <= high, axis=1)
        if not keep.any():
            continue
        recipient = Mesh(mesh.vertices, mesh.faces[keep], mesh.frame)
        for counter, p in zip(counts, primary):
            pairs = _pairs(p, recipient)
            counter.update((str(members[ids[i]]), 'fixed:'+name) for i in pairs[:, 0])
    rows = [dict(owners=list(key), before=counts[0][key], after=counts[1][key])
            for key in sorted(counts[0].keys() | counts[1].keys())]
    return dict(complete=True, passed=all(r['after'] <= r['before'] for r in rows),
                changedFaces=len(ids), fixedMeshes=len(fixed_names), rows=rows,
                anatomicalAcceptance=False)
