"""Unfold an explicitly selected free source tip without changing its shape.

No endpoint is selected by proximity. The caller owns the source/recipient
dependency inventory and must check whole-mesh contacts after this proposal.
"""
from collections import Counter
import numpy as np
from scipy.optimize import linprog

from .source_tube_topology import recover_regular_source_annuli
from ..oriented_rays import boundary_loops
from ..model import Mesh


def terminal_profile_pose(mesh, tip_ids):
    tip = np.asarray(tip_ids, int)
    if tip.ndim != 1 or len(tip) < 6 or np.any(tip < 0) or np.any(tip >= len(mesh.vertices)):
        raise ValueError('Explicit actual source tip controls required')
    topology = recover_regular_source_annuli(mesh.vertices, mesh.faces,
                                           np.arange(len(mesh.faces)//2*2).reshape(-1,2))
    inverse = np.asarray(topology['raw_to_physical']); physical = np.asarray(topology['physical_points'])
    selected = set(inverse[tip]); assert len(selected) == 6
    assert set(tip) == set(np.flatnonzero(np.isin(inverse, list(selected)))), 'Every exact source alias belongs to the tip'
    edges, _, _, _ = boundary_loops(Mesh(physical, inverse[mesh.faces], mesh.frame))
    assert selected <= set(edges.ravel()), 'An internal profile is not a free terminal'
    rid = next(i for i, ring in enumerate(topology['rings']) if set(ring) == selected)
    strips = [s for s in topology['strips'] if rid in s['rings']]; assert len(strips) == 1
    strip = strips[0]; other_id = next(i for i in strip['rings'] if i != rid)
    ring = np.asarray(topology['rings'][rid]); other = set(topology['rings'][other_id])
    mapping = {}; f = inverse[mesh.faces]; face_ids = set(strip['source_face_ids'])
    for i in range(0, len(f)-1, 2):
        if i not in face_ids: continue
        assert i+1 in face_ids
        counts = Counter(tuple(sorted((int(a),int(b)))) for tri in f[i:i+2]
                         for a,b in zip(tri,np.roll(tri,-1)))
        for (a,b), count in counts.items():
            if count == 1 and ((a in selected and b in other) or (b in selected and a in other)):
                a,b = (a,b) if a in selected else (b,a)
                assert a not in mapping or mapping[a] == b; mapping[a] = b
    assert set(mapping) == selected and len(set(mapping.values())) == 6
    p = physical[ring]; q = physical[[mapping[int(i)] for i in ring]]
    diameter = lambda v: float(np.linalg.norm(v[:,None]-v[None],axis=2).max())
    assert diameter(p) <= .00025 and diameter(q) > 2*diameter(p), 'Scope is the original collapsed tip only'
    center = p.mean(0); u,_,vt = np.linalg.svd((q-q.mean(0)).T@(p-center))
    orientation = np.eye(3); orientation[-1,-1] = np.linalg.det(u@vt)
    rotation = u@orientation@vt
    assert np.linalg.det(rotation) > 1-1e-12
    v = mesh.vertices.copy(); v[tip] = center+(v[tip]-center)@rotation.T
    return v, dict(tipIds=tip.tolist(), rotation=rotation.tolist(), center_m=center.tolist(),
                   orderedTipRawIds=[int(np.flatnonzero(inverse==i)[0]) for i in ring],
                   correspondingAdjacentRawIds=[int(np.flatnonzero(inverse==mapping[int(i)])[0]) for i in ring],
                   diameter_m=diameter(p), adjacentDiameter_m=diameter(q),
                   sourceFaceIds=strip['source_face_ids'], profileShapePreserved=True,
                   exactPositionAliasesIncluded=True, clinicalAttachmentEstablished=False)


def reconstructed_terminal_profile(mesh, tip_ids):
    """Replace a demonstrably malformed free tip by its source neighbour's shape.

    The original centre, maximum diameter, source correspondence and remote
    anatomy are fixed. The invalid original area is reported, not claimed
    conserved. No disk, vessel branch or subdivision is added.
    """
    _, pose = terminal_profile_pose(mesh, tip_ids)
    tip = np.asarray(tip_ids, int)
    ordered = np.asarray(pose['orderedTipRawIds'])
    adjacent = np.asarray(pose['correspondingAdjacentRawIds'])
    p, q = mesh.vertices[ordered], mesh.vertices[adjacent]
    center = p.mean(0)
    def projected_shape(points):
        z = points - points.mean(0)
        _, singular, axes = np.linalg.svd(z, full_matrices=False)
        xy = z @ axes[:2].T
        edges = np.roll(xy, -1, axis=0) - xy
        second = np.roll(edges, -1, axis=0)
        turns = edges[:,0]*second[:,1] - edges[:,1]*second[:,0]
        convex = bool(np.all(turns > 1e-20) or np.all(turns < -1e-20))
        area = float(abs(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1))) / 2)
        return dict(convex=convex, singularValues_m=singular.tolist(), projectedArea_m2=area)
    before, neighbor = projected_shape(p), projected_shape(q)
    if before['convex'] and before['singularValues_m'][1] / before['singularValues_m'][0] >= .05:
        raise ValueError('A valid source tip must not be reconstructed')
    if not neighbor['convex']:
        raise ValueError('The adjoining source profile is not a valid convex template')
    scale = pose['diameter_m'] / pose['adjacentDiameter_m']
    points = center + (q - q.mean(0)) * scale
    result = mesh.vertices.copy()
    for source, target in zip(p, points):
        aliases = np.flatnonzero(np.all(mesh.vertices == source, axis=1))
        assert set(aliases) <= set(tip)
        result[aliases] = target
    np.testing.assert_allclose(points.mean(0), center, atol=5e-16, rtol=0)
    diameter = float(np.linalg.norm(points[:,None]-points[None],axis=2).max())
    assert abs(diameter-pose['diameter_m']) < 5e-16
    return result, dict(method='homothetic adjoining source profile', tipIds=tip.tolist(),
                        orderedTipRawIds=ordered.tolist(), correspondingAdjacentRawIds=adjacent.tolist(),
                        sourceFaceIds=pose['sourceFaceIds'], center_m=center.tolist(),
                        diameter_m=diameter, adjacentDiameter_m=pose['adjacentDiameter_m'],
                        scale=scale, originalProfile=before, adjoiningProfile=neighbor,
                        rebuiltProfile=projected_shape(points), originalInvalidAreaPreserved=False,
                        centreAndMaximumDiameterPreserved=True, clinicalAttachmentEstablished=False)


def unfolded_tip_words(old, vertices, tip):
    """Reconstruct only changed stars; leave pre-existing remote faults explicit."""
    active = np.unique(old.faces[np.any(np.isin(old.faces,tip),axis=1)])
    t = vertices[old.faces]; area = np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
    lengths = np.linalg.norm(area,axis=1); assert np.all(lengths > 0)
    facet = area/lengths[:,None]; sums = np.zeros_like(vertices)
    for k in range(3): np.add.at(sums,old.faces[:,k],area)
    signs = old.attributes['normalSigns']; words = old.attributes['normals'].copy(); solves = 0
    for vertex in active:
        incident = facet[np.any(old.faces==vertex,axis=1)]
        direction = sums[vertex]/np.linalg.norm(sums[vertex])
        if np.min(incident@direction) <= 1e-4:
            fit = linprog([0,0,0,-1],A_ub=np.c_[-incident,np.ones(len(incident))],
                          b_ub=np.zeros(len(incident)),bounds=[(-1,1)]*3+[(0,1)],
                          method='highs',options={'time_limit':.5})
            assert fit.success and fit.x[3] > 1e-4, 'No strictly outward common normal cone'
            direction = fit.x[:3]/np.linalg.norm(fit.x[:3]); solves += 1
        words[vertex] = np.rint(direction*32767*signs[vertex]).astype('<i2')
        assert np.all(incident@(words[vertex].astype(float)*signs[vertex]) > 0)
    def bad(v,w):
        t = v[old.faces]; n = np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        effective = w.astype(float)*signs[:,None]
        return set(np.flatnonzero(np.einsum('ij,ij->i',n,effective[old.faces].mean(1))<=0))
    before = bad(old.vertices,old.attributes['normals']); after = bad(vertices,words)
    assert after <= before, 'A previously outward source facet cannot become inward'
    return words, dict(method='final incident-facet field with strictly outward local cones',
                       changedStars=active.tolist(), outwardConeSolves=solves,
                       preexistingInwardAverageFields=len(before), remainingInwardAverageFields=len(after),
                       unchangedStarsLiteral=True, newInwardAverageFields=0)
