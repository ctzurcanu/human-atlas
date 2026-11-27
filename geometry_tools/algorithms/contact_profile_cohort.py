"""Transport existing source contacts together without merging their meshes.

These constraints preserve an existing contact during a placement stage. They
do not establish a valid anatomical junction or excuse artery-vein overlap.
"""
import numpy as np


def contact_profile_cohort(group_counts, face_groups, contacts, aliases=(), fixed=()):
    offsets = np.r_[0, np.cumsum(group_counts)]
    parent = np.arange(int(offsets[-1]))

    def root(i):
        i = int(i)
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = int(parent[i])
        return i

    def merge(entries):
        entries = list(entries)
        if not entries:
            raise ValueError('An empty contact cannot couple source profiles')
        first = root(entries[0])
        for entry in entries[1:]:
            parent[root(entry)] = first

    for first, second, pairs in contacts:
        for i, j in np.asarray(pairs, dtype=int).reshape(-1, 2):
            merge(list(offsets[first] + face_groups[first][i]) +
                  list(offsets[second] + face_groups[second][j]))
    for entries in aliases:
        merge(offsets[actor] + group for actor, group in entries)
    labels = np.array([root(i) for i in range(len(parent))])
    _, labels = np.unique(labels, return_inverse=True)
    mappings = [labels[offsets[i]:offsets[i+1]] for i in range(len(group_counts))]
    protected = sorted({int(mappings[actor][group]) for actor, group in fixed})
    return dict(clusterOfGroup=[v.tolist() for v in mappings],
                clusters=int(labels.max())+1, fixedClusters=protected,
                sourceMeshesMerged=False, anatomicalInterfaceEstablished=False)
