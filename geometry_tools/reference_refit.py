"""Source-bound refits of existing open reference sheets, never new tissues."""
import numpy as np
from .audit import topology
from .common import digest_file
from .model import Mesh


def verify_reference_refit(original, candidate, policy):
    if policy.get('schema') == 'open-source-aperture/v1':
        return verify_reference_aperture(original, candidate, policy)
    if policy.get('schema') != 'open-source-z-refit/v1':
        raise ValueError('An explicit open-reference refit policy is required')
    if original.frame != candidate.frame or original.geometry_sha != policy['source_geometry_sha256']:
        raise ValueError('Reference refit source or frame changed')
    for key in ['refined_source', 'lineage', 'controls']:
        if digest_file(policy[key]) != policy[key+'_sha256']:
            raise ValueError('Reference refit evidence changed: '+key)
    baseline = Mesh.load(policy['refined_source'])
    with np.load(policy['lineage'], allow_pickle=False) as z:
        parents = z['parent_faces'].copy()
        original_ids = z['original_vertices'].copy()
    with np.load(policy['controls'], allow_pickle=False) as z:
        pinned = z['pinned'].copy()
    if baseline.frame != original.frame or not np.array_equal(original_ids, np.arange(len(original.vertices))):
        raise ValueError('Refined source must retain the literal original control prefix')
    if not np.array_equal(baseline.vertices[original_ids], original.vertices):
        raise ValueError('Source refinement moved original controls')
    if (not np.array_equal(candidate.faces, baseline.faces)
        or candidate.vertices.shape != baseline.vertices.shape
        or parents.shape != (len(baseline.faces),) or parents.dtype.kind not in 'iu'
        or np.any(parents < 0) or np.any(parents >= len(original.faces))
        or pinned.dtype.kind != 'b' or pinned.shape != (len(baseline.vertices),)):
        raise ValueError('Refit topology, parent lineage or pinned controls changed')
    # Positive coplanar child areas must cover every source parent once.
    t = original.triangles[parents]; q = baseline.triangles
    u, v = t[:,1]-t[:,0], t[:,2]-t[:,0]
    normal = np.cross(u, v); denominator = np.einsum('ij,ij->i', normal, normal)
    if np.any(denominator <= 0): raise ValueError('Degenerate reference parent')
    d = q-t[:,None,0]
    plane = np.abs(np.einsum('fci,fi->fc', d, normal))/np.sqrt(denominator[:,None])
    a = np.einsum('fci,fi->fc', np.cross(d, v[:,None]), normal)/denominator[:,None]
    b = np.einsum('fci,fi->fc', np.cross(u[:,None], d), normal)/denominator[:,None]
    oriented = np.einsum('ij,ij->i', np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]), normal)/denominator
    sums = np.bincount(parents, weights=oriented, minlength=len(original.faces))
    if np.any(plane > 1e-12) or np.any(a < -1e-10) or np.any(b < -1e-10) or np.any(a+b > 1+1e-10) or np.any(oriented <= 0) or not np.allclose(sums,1,atol=1e-9,rtol=0):
        raise ValueError('Refinement does not preserve every original triangle')
    if 'uvs' in original.attributes:
        if not np.array_equal(baseline.attributes.get('uvs',[])[:len(original_ids)], original.attributes['uvs']) or not np.array_equal(candidate.attributes.get('uvs'), baseline.attributes.get('uvs')):
            raise ValueError('Reference refit changed source UVs')
    tolerance = float(policy['encoding_tolerance_m']); maximum = float(policy['maximum_z_shift_m'])
    if not (0 <= tolerance <= 1e-7 and 0 < maximum <= .01):
        raise ValueError('Bound the declared encoding error and local displacement')
    delta = candidate.vertices-baseline.vertices
    if np.max(np.abs(delta[:,:2]),initial=0) > tolerance or np.max(np.abs(delta[:,2]),initial=0) > maximum+tolerance or np.max(np.abs(delta[pinned]),initial=0) > tolerance:
        raise ValueError('Refit exceeds its declared field or fixed controls')
    if not np.array_equal(candidate.vertices[original_ids[pinned[original_ids]]], original.vertices[pinned[original_ids]]):
        raise ValueError('A fixed original reference control changed')
    topologies = [topology(m) for m in [original,baseline,candidate]]
    signature = lambda t: sorted((c['euler'],bool(c['boundary_edges'])) for c in t['components'])
    if not topologies[0]['open_edges'] or any(signature(t)!=signature(topologies[0]) or t['nonmanifold_edges'] or t['inconsistent_winding_edges'] or t['invalid_vertex_links'] or t['zero_area_faces'] or t['duplicate_faces'] for t in topologies):
        raise ValueError('An existing open-sheet refit must preserve its source topology')
    return dict(passed=True, existing_open_reference_only=True, source_geometry_sha256=original.geometry_sha,
        original_parent_coverage_replayed=True, source_uvs_preserved=True, original_fixed_controls_literal=True,
        maximum_z_shift_m=float(np.max(np.abs(delta[:,2]),initial=0)), new_finite_tissue_claim=False)


def verify_reference_aperture(original, candidate, policy):
    """Remove one declared source-face disk, backed by a finite annular collar.

    This policy does not turn the remaining source sheet into a finite vessel.
    Geometry and every retained attribute corner outside the window are literal;
    all original mouths survive. The new boundary belongs to a separately
    closed, source-bound collar with a real hole, never a solid opening marker.
    """
    from .algorithms.tube_continuity import boundary_loops
    from .contact import pair_contacts
    if policy.get('schema') != 'open-source-aperture/v1':
        raise ValueError('An explicit source-aperture policy is required')
    if original.frame != candidate.frame or original.geometry_sha != policy['source_geometry_sha256']:
        raise ValueError('Aperture source or frame changed')
    for key in ['evidence', 'collar']:
        if digest_file(policy[key]) != policy[key+'_sha256']:
            raise ValueError('Aperture evidence changed: '+key)
    import json
    with open(policy['evidence']) as source:
        evidence=json.load(source)
    if evidence.get('schema') != 'open-source-aperture-evidence/v1':
        raise ValueError('Source-bound aperture evidence is required')
    if evidence.get('source_geometry_sha256') != original.geometry_sha:
        raise ValueError('Aperture evidence belongs to another source')
    removed=np.asarray(evidence.get('removed_source_faces',[]))
    if (removed.ndim != 1 or removed.dtype.kind not in 'iu' or not len(removed)
        or len(np.unique(removed)) != len(removed) or removed.min()<0 or removed.max()>=len(original.faces)):
        raise ValueError('One explicit nonempty source-face window is required')
    retained=np.ones(len(original.faces),bool);retained[removed]=False
    if not retained.any() or not np.array_equal(candidate.triangles,original.triangles[retained]):
        raise ValueError('Every face outside the aperture must remain literal and in source order')
    if original.attributes.keys() != candidate.attributes.keys():
        raise ValueError('Aperture changes source attributes')
    for key, values in original.attributes.items():
        if not np.array_equal(candidate.attributes[key][candidate.faces],values[original.faces[retained]]):
            raise ValueError('Aperture changes retained source attribute corners: '+key)
    patch=Mesh(original.vertices,original.faces[removed],original.frame)
    topologies=[topology(m) for m in [original,candidate,patch]]
    for t in topologies:
        if (t['nonmanifold_edges'] or t['inconsistent_winding_edges'] or t['invalid_vertex_links']
            or t['zero_area_faces'] or t['duplicate_faces']):
            raise ValueError('Invalid source/aperture topology')
    if not topologies[0]['open_edges'] or len(topologies[0]['components'])!=1 or len(topologies[1]['components'])!=1:
        raise ValueError('Apertures require an existing connected open source sheet')
    if len(topologies[2]['components'])!=1 or topologies[2]['components'][0]['euler']!=1:
        raise ValueError('Removed source window must be one topological disk')
    loops=[boundary_loops(m.faces) for m in [original,candidate,patch]]
    signature=lambda mesh,loop:tuple(sorted(map(tuple,mesh.vertices[loop])))
    before={signature(original,loop) for loop in loops[0]}
    after={signature(candidate,loop) for loop in loops[1]}
    patchloops={signature(patch,loop) for loop in loops[2]}
    if len(patchloops)!=1 or after != before | patchloops or len(after)!=len(before)+1:
        raise ValueError('Aperture must preserve every old mouth and add only the declared source boundary')
    collar=Mesh.load(policy['collar'])
    if collar.frame!=original.frame or collar.geometry_sha!=evidence.get('collar_geometry_sha256'):
        raise ValueError('Aperture collar source or frame changed')
    ct=topology(collar)
    if (not ct['closed_oriented_manifold'] or len(ct['components'])!=1
        or ct['components'][0]['closed_orientable_genus']!=1):
        raise ValueError('The local aperture collar must be a finite annular tissue cell')
    # A collection of matching points does not establish a communicating seam.
    # Require every literal source-window edge to be owned by the collar's
    # outer chart and its inward thickness rim, with opposite directed uses.
    from collections import defaultdict
    edges=defaultdict(list)
    for fi,triangle in enumerate(collar.triangles):
        for a,b in zip(triangle,np.roll(triangle,-1,axis=0)):
            edges[tuple(sorted((tuple(a),tuple(b))))].append((fi,tuple(a),tuple(b)))
    seam=[]
    for a,b in zip(loops[2][0],np.roll(loops[2][0],-1)):
        aa,bb=tuple(patch.vertices[a]),tuple(patch.vertices[b])
        uses=edges.get(tuple(sorted((aa,bb))),[])
        if len(uses)!=2 or {(u[1],u[2])for u in uses}!={(aa,bb),(bb,aa)}:
            raise ValueError('Aperture boundary is not a literal shared collar edge')
        source_faces=np.flatnonzero(np.any(patch.faces==a,axis=1)&np.any(patch.faces==b,axis=1))
        if len(source_faces)!=1:raise ValueError('Source-window boundary edge ownership is ambiguous')
        q=patch.triangles[source_faces[0]];normal=np.cross(q[1]-q[0],q[2]-q[0]);normal/=np.linalg.norm(normal)
        distance=[(collar.triangles[u[0]]-q[0])@normal for u in uses]
        outer=[i for i,d in enumerate(distance)if np.max(abs(d))<=1e-12 and uses[i][1:]==(aa,bb)]
        if len(outer)!=1 or np.min(distance[1-outer[0]])>=-1e-12 or np.max(distance[1-outer[0]])>1e-12:
            raise ValueError('Aperture edge lacks its source-tangent outer chart and inward thickness rim')
        seam.append(dict(source_edge=[int(a),int(b)],collar_outer_face=int(uses[outer[0]][0]),
                         collar_rim_face=int(uses[1-outer[0]][0])))
    from .contact import self_contact_summary
    collar_self,_=self_contact_summary(collar,max_tests=2000000)
    if (not collar_self['complete_within_pair_budget'] or collar_self['crossing_face_pairs']
        or collar_self['coplanar_overlap_face_pairs']):
        raise ValueError('Aperture finite collar self-intersects')
    _,contact=pair_contacts(candidate,collar,0.,max_pairs=250000,max_tests=2000000)
    if (not contact['complete_within_pair_budget'] or contact['crossing_face_pairs']
        or contact['coplanar_overlap_face_pairs']):
        raise ValueError('Aperture collar overlaps or crosses the retained source sheet')
    return dict(passed=True,existing_open_reference_only=True,source_geometry_sha256=original.geometry_sha,
        source_aperture_policy=True,removed_source_faces=removed.tolist(),retained_faces_literal=int(retained.sum()),
        retained_attribute_corners_literal=True,old_source_mouths_preserved=len(before),new_source_boundaries=1,
        finite_annular_collar_geometry_sha256=collar.geometry_sha,collar_reference_contact=contact,
        literal_source_window_edge_ownership=seam,collar_self_contact=collar_self,
        remaining_parent_is_open_reference=True,new_whole_parent_finite_tissue_claim=False)
