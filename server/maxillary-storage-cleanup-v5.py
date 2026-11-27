"""Bounded link-condition collapse of diagnosed new CSG storage slivers.

This helper has no CLI side effects. It cannot move an original source control
or collapse an edge longer than the explicit ten-nanometre physical budget.
"""
from collections import Counter
import numpy as np
from geometry_tools.audit import topology
from geometry_tools.model import Mesh


def cleanup(vertices, faces, normals, uvs, owners, ancestry, failed, frame):
    maximum = 1e-8
    v = vertices.copy(); f = faces.copy(); n = normals.copy(); uv = uvs.copy()
    prefix = len(ancestry['source_position_prefix'])
    original_v, original_n, original_uv = v[:prefix].copy(), n[:prefix].copy(), uv[:prefix].copy()
    origin = np.array([0., 1.671875, 0.])
    quantized = (v - origin).astype('<f4').astype(float) + origin
    triangles = quantized[f]
    zero = np.flatnonzero(np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0],
                                                 triangles[:, 2] - triangles[:, 0]), axis=1) == 0)
    bad = set(map(int, zero)) | {int(i) for pair in failed['self_intersections']['examples'] for i in pair}
    gv, inverse = np.unique(v, axis=0, return_inverse=True)
    gf = inverse[f]
    old_topology = topology(Mesh(gv, gf, frame))
    protected = set(map(int, inverse[:prefix]))
    support_triangles = ancestry['parent_triangles'][ancestry['face_kinds'] == 'junction_faces']
    assert len(support_triangles), 'Explicit constructed junction supports are required'
    low, high = support_triangles.min(1) - maximum, support_triangles.max(1) + maximum
    def supported(point):
        return bool(np.any(np.all(point >= low, axis=1) & np.all(point <= high, axis=1)))
    witness_controls = np.unique(gf[sorted(bad)])
    edges = np.unique(np.sort(np.vstack([gf[:, [0, 1]], gf[:, [1, 2]], gf[:, [2, 0]]]), axis=1), axis=0)
    adjacent = np.isin(edges, witness_controls).any(1)
    lengths = np.linalg.norm(gv[edges[:, 0]] - gv[edges[:, 1]], axis=1)
    candidates = [(float(length), int(a), int(b)) for (a, b), length, near in zip(edges, lengths, adjacent)
                  if near and length <= maximum and not ({int(a), int(b)} <= protected)
                  and supported(gv[a]) and supported(gv[b])]
    candidates.sort()
    retained_faces = np.arange(len(faces))
    groups = {i: {i} for i in range(len(gv))}
    mapping = np.arange(len(gv))
    collapsed, skipped = [], []
    for raw_length, aa, bb in candidates:
        a, b = int(mapping[aa]), int(mapping[bb])
        if a == b:
            continue
        if a in protected and b in protected:
            skipped.append(dict(edge=[a, b], reason='both source controls protected')); continue
        target, removed = (a, b) if a in protected else (b, a) if b in protected else (min(a, b), max(a, b))
        touching_a, touching_b = np.any(gf == a, axis=1), np.any(gf == b, axis=1)
        edge_faces = np.flatnonzero(touching_a & touching_b)
        if len(edge_faces) != 2:
            skipped.append(dict(edge=[a, b], reason='not a two-face interior edge')); continue
        neighbors_a = set(map(int, gf[touching_a].ravel())) - {a}
        neighbors_b = set(map(int, gf[touching_b].ravel())) - {b}
        opposite = set(map(int, gf[edge_faces].ravel())) - {a, b}
        if neighbors_a & neighbors_b != opposite or len(opposite) != 2:
            skipped.append(dict(edge=[a, b], reason='link condition fails')); continue
        members = groups[target] | groups[removed]
        aliases = np.flatnonzero(np.isin(inverse, list(members)))
        destination = quantized[np.flatnonzero(inverse == target)[0]]
        error = np.linalg.norm(vertices[aliases] - destination, axis=1)
        if error.max(initial=0) > maximum:
            skipped.append(dict(edge=[a, b], reason='physical error exceeds ten nanometres')); continue
        changed_aliases = aliases[np.linalg.norm(v[aliases] - destination, axis=1) > 0]
        if np.any(changed_aliases < prefix):
            skipped.append(dict(edge=[a, b], reason='would move an original source position')); continue
        replacement = gf.copy(); replacement[replacement == removed] = target
        degenerate = np.any(np.sort(replacement, axis=1)[:, 1:] == np.sort(replacement, axis=1)[:, :-1], axis=1)
        if set(np.flatnonzero(degenerate)) != set(edge_faces):
            skipped.append(dict(edge=[a, b], reason='unexpected face incidence collapse')); continue
        kept = ~degenerate
        canonical_faces = np.sort(replacement[kept], axis=1)
        if len(np.unique(canonical_faces, axis=0)) != len(canonical_faces):
            skipped.append(dict(edge=[a, b], reason='would create a duplicate face')); continue
        deleted = retained_faces[edge_faces].copy()
        affected = retained_faces[np.any(np.isin(gf, list(members)), axis=1)].copy()
        v[aliases] = destination
        # Source words/UVs remain exact. Only new intersection aliases at the
        # collapsed location receive a single normal; UV charts stay distinct.
        protected_aliases = aliases[aliases < prefix]
        words = np.unique(n[protected_aliases] if len(protected_aliases) else n[aliases], axis=0)
        assert not len(protected_aliases) or len(words) == 1
        if len(words) == 1:
            word = words[0]
        else:
            word = words.astype(float).sum(0); assert np.linalg.norm(word) > 0
            word = np.rint(word / np.linalg.norm(word) * 32767).astype(np.int16)
        n[aliases] = word
        gf, retained_faces = replacement[kept], retained_faces[kept]
        groups[target] = members; groups.pop(removed)
        mapping[np.isin(mapping, list(members))] = target
        collapsed.append(dict(edge_before=[a, b], target_control=target, removed_control=removed,
            raw_edge_length_m=raw_length, maximum_physical_error_m=float(error.max()),
            changed_new_chart_vertex_ids=changed_aliases.tolist(),
            removed_original_candidate_face_ids=deleted.tolist(),
            affected_original_candidate_face_ids=affected.tolist(),
            link_condition_passed=True, strict_junction_support_passed=True))
    assert np.array_equal(v[:prefix], original_v)
    assert np.array_equal(n[:prefix], original_n) and np.array_equal(uv[:prefix], original_uv)
    f, owners = faces[retained_faces], owners[retained_faces]
    raw_topology = topology(Mesh(v, f, frame))
    assert not any(raw_topology[k] for k in ['nonmanifold_edges', 'inconsistent_winding_edges',
                                           'invalid_vertex_links', 'duplicate_faces', 'zero_area_faces'])
    assert [(c['euler'], c['boundary_edges']) for c in raw_topology['components']] == [(c['euler'], c['boundary_edges']) for c in old_topology['components']]
    report = dict(schema='storage-aware-maxillary-junction-cleanup/v5', complete=True,
        executed_single_deterministic_edge_pass=True, maximum_authorized_error_m=maximum,
        initial_zero_faces=zero.tolist(), initial_crossing_witness_face_ids=sorted(bad - set(map(int, zero))),
        candidate_short_edges=len(candidates), collapsed_edges=collapsed, skipped_edges=skipped,
        all_original_position_normal_UV_prefix_literal=True, strict_junction_supports=len(support_triangles),
        maximum_physical_error_m=max([c['maximum_physical_error_m'] for c in collapsed], default=0.),
        removed_candidate_face_ids=sorted(set(range(len(faces))) - set(map(int, retained_faces))),
        retained_candidate_face_ids=retained_faces.tolist(), raw_topology_preserved=True,
        actual_Float32_gate_passed=False, installed=False)
    return v, f, n, uv, owners, report
