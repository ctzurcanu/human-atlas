"""Least-squares projection of a shared displacement onto rigid source profiles.

Exact-position property aliases and intersecting profile ownership are coupled.
No ring is rescaled or discarded. This is a discrete proposal, not a certified
continuous map; emitted faces still need full interface/contact/quality checks.
"""
import numpy as np


def project_rigid_profiles(meshes, proposed, rings, fixed_masks, *, allow_rotation=False):
    if not meshes or not (len(meshes) == len(proposed) == len(rings) == len(fixed_masks)):
        raise ValueError('Corresponding source meshes, proposals, profiles and fixed masks required')
    sizes = np.cumsum([0]+[len(m.vertices) for m in meshes])
    original = np.concatenate([m.vertices for m in meshes])
    target = np.concatenate(proposed)
    fixed = np.concatenate(fixed_masks)
    if target.shape != original.shape or fixed.shape != (len(original),):
        raise ValueError('Every original control needs a proposal and explicit fixed policy')
    if not np.isfinite(target).all():
        raise ValueError('Finite source proposals required')
    physical, inverse = np.unique(original, axis=0, return_inverse=True)
    parents = np.arange(len(physical))
    def find(i):
        while parents[i] != i:
            parents[i] = parents[parents[i]]
            i = parents[i]
        return i
    profiles = 0
    for offset, mesh_rings in zip(sizes[:-1], rings):
        for ring in mesh_rings:
            ids = inverse[offset+np.asarray(ring, int)]
            if len(ids) < 3 or len(np.unique(ids)) != len(ids):
                raise ValueError('Each declared source profile must contain distinct controls')
            for j in ids[1:]:
                parents[find(int(j))] = find(int(ids[0]))
            profiles += 1
    _, owner = np.unique([find(i) for i in range(len(physical))], return_inverse=True)
    raw_owner = owner[inverse]
    delta = target-original
    physical_count = np.bincount(inverse)
    physical_delta = np.column_stack([np.bincount(inverse, weights=delta[:, k])/physical_count for k in range(3)])
    mismatch = float(np.linalg.norm(delta-physical_delta[inverse], axis=1).max())
    if mismatch > 1e-12:
        raise ValueError('Coincident source controls received inconsistent proposed displacements')
    count = np.bincount(owner)
    displacement = np.column_stack([np.bincount(owner, weights=physical_delta[:, k])/count for k in range(3)])
    fixed_groups = np.unique(raw_owner[fixed])
    displacement[fixed_groups] = 0
    result = original+displacement[raw_owner]
    if allow_rotation:
        # Orthogonal Procrustes: the closest proper rigid fit to the saved
        # guide, instead of holding the old profile orientation during a turn.
        projected = physical.copy()
        for group in range(len(count)):
            ids = np.flatnonzero(owner == group)
            if group in fixed_groups:
                continue
            p = physical[ids]
            q = p+physical_delta[ids]
            first, second = p.mean(0), q.mean(0)
            if len(ids) < 3:
                projected[ids] = p+displacement[group]
                continue
            u, _, vt = np.linalg.svd((p-first).T@(q-second))
            correction = np.diag([1., 1., np.linalg.det(u@vt)])
            rotation = u@correction@vt
            projected[ids] = (p-first)@rotation+second
        result = projected[inverse]
    assert np.array_equal(result[fixed], original[fixed])
    emitted = [result[a:b] for a,b in zip(sizes[:-1], sizes[1:])]
    return emitted, dict(sourceProfiles=profiles, coupledRigidGroups=len(count),
                         physicalControls=len(physical), rawControls=len(original),
                         fixedGroups=len(fixed_groups), coincidentProposalDiscrepancy_m=mismatch,
                         allSourceProfilesMovedRigidlyWithoutScaling=True,
                         properRotationsAllowed=allow_rotation,
                         exactPositionAliasesCoupled=True, continuousMapCertified=False)


def transport_corner_normals(mesh, vertices):
    """Keep original shading through the proper rotation of geometric normals."""
    _, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
    normals = []
    for points in (mesh.vertices, vertices):
        t = points[mesh.faces]
        area = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
        n = np.zeros((int(inverse.max())+1, 3))
        for k in range(3):
            np.add.at(n, inverse[mesh.faces[:, k]], area)
        length = np.linalg.norm(n, axis=1)
        if np.any(length <= 1e-20):
            raise ValueError('Ambiguous source geometric normal cannot be transported')
        normals.append((n/length[:, None])[inverse])
    first, second = normals
    axis = np.cross(first, second)
    cosine = np.einsum('ij,ij->i', first, second)
    if np.any(cosine <= -1+1e-8):
        raise ValueError('A reversed geometric normal requires separate source review')
    attributes = {k:v.copy() for k,v in mesh.attributes.items()}
    n = attributes['normals'].astype(float)*attributes['normalSigns'][:, None]
    n /= np.linalg.norm(n, axis=1)[:, None]
    n += np.cross(axis, n)+np.cross(axis, np.cross(axis, n))/(1+cosine[:, None])
    n /= np.linalg.norm(n, axis=1)[:, None]
    moved = np.any(vertices != mesh.vertices, axis=1)
    attributes['normals'][moved] = np.rint(n[moved]*32767).astype('<i2')
    attributes['normalSigns'][moved] = 1
    return attributes
