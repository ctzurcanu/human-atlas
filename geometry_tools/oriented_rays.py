"""Batched oriented all-hit ray diagnostics on literal, possibly open meshes.

Native Float32 rays only propose face IDs. Every event is replayed against the
original Float64 triangle. Open boundary shadows, grazing/edge events and ray
disagreement remain uncertain; no closure or source mutation is performed.
"""
import time
import numpy as np
from .open3d_backend import SurfaceIndex, library, _points


def boundary_loops(mesh):
    f = mesh.faces
    edges = np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]])
    keys, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0,
                                      return_inverse=True, return_counts=True)
    balance = np.bincount(inverse, weights=np.where(edges[:, 0] < edges[:, 1], 1, -1))
    boundary = keys[counts == 1]
    graph = {}
    for a, b in boundary:
        graph.setdefault(int(a), []).append(int(b)); graph.setdefault(int(b), []).append(int(a))
    bad = {k: len(v) for k, v in graph.items() if len(v) != 2}
    unseen = set(graph); loops = []; networks = []
    while unseen:
        first = min(unseen); pending = [first]; vertices = set()
        while pending:
            node = pending.pop()
            if node in vertices: continue
            vertices.add(node); unseen.discard(node); pending.extend(graph[node])
        ids = np.asarray(sorted(vertices), dtype=np.int64)
        simple = all(len(graph[node]) == 2 for node in vertices)
        if simple:
            a = first; b = graph[a][0]; loop = [a]
            while b != first:
                loop.append(b); a, b = b, next(v for v in graph[b] if v != a)
            loops.append(np.asarray(loop, dtype=np.int64))
        networks.append(dict(vertex_ids=ids, simple_cycle=simple))
    network_summary = [dict(network=i, boundary_vertices=len(net['vertex_ids']),
        boundary_edges=sum(len(graph[int(v)]) for v in net['vertex_ids'])//2,
        simple_cycle=net['simple_cycle'], bounds_m=np.array([
            mesh.vertices[net['vertex_ids']].min(0), mesh.vertices[net['vertex_ids']].max(0)]).tolist())
        for i, net in enumerate(networks)]
    return boundary, loops, networks, dict(boundary_edges=len(boundary), closed_boundary_loops=len(loops),
        connected_boundary_networks=len(networks), boundary_networks=network_summary,
        non_degree_two_boundary_vertices=len(bad), nonmanifold_edges=int((counts > 2).sum()),
        inconsistent_shared_edge_winding=int(((counts == 2) & (balance != 0)).sum()))


def literal_hits(origins, directions, triangles, *, bary_epsilon=1e-10, forward_m=1e-9):
    """Vectorized corresponding ray/triangle predicate, not an all-pairs loop."""
    a = triangles[:, 0]; e1 = triangles[:, 1] - a; e2 = triangles[:, 2] - a
    h = np.cross(directions, e2); det = np.einsum('ij,ij->i', e1, h)
    normal = np.cross(e1, e2); length = np.linalg.norm(normal, axis=1)
    valid = np.abs(det) > np.maximum(length * 1e-13, 1e-30)
    inv = np.zeros(len(det)); inv[valid] = 1 / det[valid]
    s = origins - a; u = inv * np.einsum('ij,ij->i', s, h)
    q = np.cross(s, e1); v = inv * np.einsum('ij,ij->i', directions, q)
    travel = inv * np.einsum('ij,ij->i', e2, q)
    bary = np.c_[1-u-v, u, v]
    valid &= (bary.min(1) >= -bary_epsilon) & (travel > forward_m)
    dot = np.einsum('ij,ij->i', normal, directions)
    cosine = np.divide(dot, length, out=np.zeros_like(dot), where=length > 0)
    return valid, travel, bary, np.sign(dot).astype(np.int8), cosine


class OrientedRayIndex:
    def __init__(self, mesh, *, nthreads=2, uncertainty_segments=None,
                 uncertainty_component_ids=None, boundary_evidence=None, boundary_tolerance_m=1e-6):
        self.mesh = mesh; self.index = SurfaceIndex(mesh, nthreads=nthreads)
        self.edges, self.loops, self.networks, self.topology = boundary_loops(mesh)
        self.boundary_evidence = boundary_evidence
        self.boundary_tolerance_m = float(boundary_tolerance_m)
        if not np.isfinite(self.boundary_tolerance_m) or self.boundary_tolerance_m < 0:
            raise ValueError('Finite nonnegative boundary tolerance required')
        if uncertainty_segments is not None:
            segments = np.asarray(uncertainty_segments, float)
            ids = np.asarray(uncertainty_component_ids)
            if segments.ndim != 3 or segments.shape[1:] != (2, 3) or len(ids) != len(segments) or not np.isfinite(segments).all():
                raise ValueError('Source-bound uncertainty segments [N,2,3] and component IDs required')
            self.networks = [dict(points=segments[ids == i].reshape(-1, 3), component=int(i)) for i in np.unique(ids)]
        else:
            self.networks = [dict(points=mesh.vertices[n['vertex_ids']], component=i) for i, n in enumerate(self.networks)]
        self.uncertainty_network_summary = [dict(component=n['component'], controls=len(n['points']),
            bounds_m=np.array([n['points'].min(0), n['points'].max(0)]).tolist()) for n in self.networks]

    def boundary_shadow(self, points, direction, epsilon_m):
        """Conservative projected AABBs of each actual boundary network.

        This constructs no cap. A half-ray entirely beyond a network cannot
        traverse its shadow. Branched networks affect only their own shadows.
        """
        axis = np.eye(3)[np.argmin(np.abs(direction))]
        u = np.cross(direction, axis); u /= np.linalg.norm(u); v = np.cross(direction, u)
        p = np.c_[points @ u, points @ v]; shadow = np.zeros(len(p), bool)
        near = np.zeros(len(p), bool)
        for network in self.networks:
            verts = network['points']; xy = np.c_[verts @ u, verts @ v]
            lo = xy.min(0)-epsilon_m; hi = xy.max(0)+epsilon_m
            forward = points @ direction <= (verts @ direction).max() + epsilon_m
            shadow |= np.all((p >= lo) & (p <= hi), axis=1) & forward
        return shadow, near

    def classify(self, points, directions, *, near_surface_m=1e-5, min_reliable_rays=2):
        began = time.perf_counter(); p = _points(points); d = _points(directions, 'directions')
        lengths = np.linalg.norm(d, axis=1)
        if np.any(lengths == 0): raise ValueError('Nonzero ray directions required')
        d = d / lengths[:, None]; n = len(p); k = len(d)
        local = ((p-self.index.origin)/self.index.scale).astype(np.float32)
        query_rounding = np.linalg.norm(local.astype(float)*self.index.scale+self.index.origin-p, axis=1)
        eps_m = max(1e-7, self.boundary_tolerance_m, 8*(self.index.vertex_rounding_m + float(query_rounding.max(initial=0))))
        winding = np.zeros((n, k), np.int16); counts = np.zeros((n, k), np.int16)
        numerical = np.zeros((n, k), bool); shadow = np.zeros((n, k), bool)
        edge = np.zeros((n, k), bool); native_hits = 0; invalid_hits = 0; events = []
        for j, direction in enumerate(d):
            rays = np.c_[local, np.broadcast_to(direction, (n, 3))].astype(np.float32)
            raw = {key: val.numpy() for key, val in self.index.scene.list_intersections(
                library().core.Tensor(rays), nthreads=self.index.nthreads).items()}
            ids = raw['ray_ids'].astype(np.int64); faces = raw['primitive_ids'].astype(np.int64)
            native_hits += len(ids)
            valid, travel, bary, sign, cosine = literal_hits(p[ids], np.broadcast_to(direction, (len(ids), 3)), self.mesh.triangles[faces])
            bad = ~valid; invalid_hits += int(bad.sum()); numerical[ids[bad], j] = True
            sensitive = (bary.min(1) <= 1e-5) | (np.abs(cosine) <= 1e-6) | (travel <= near_surface_m)
            edge[ids[sensitive], j] = True
            counts[:, j] = np.bincount(ids[valid], minlength=n).astype(np.int16)
            winding[:, j] = np.bincount(ids[valid], weights=sign[valid], minlength=n).astype(np.int16)
            shadow[:, j], projected_near = self.boundary_shadow(p, direction, eps_m)
            numerical[:, j] |= projected_near
            events.append(dict(direction=direction.tolist(), proposed_events=len(ids),
                rejected_literal_events=int(bad.sum()), edge_or_grazing_event_rays=int(edge[:, j].sum()),
                boundary_shadow_rays=int(shadow[:, j].sum())))
        nearest, _ = self.index.closest_points(p, certify=False)
        near_surface = nearest['distance_m'] <= near_surface_m + eps_m
        reliable = ~(numerical | edge | shadow)
        votes = reliable.sum(1); low = np.min(np.where(reliable, winding, 32767), axis=1)
        high = np.max(np.where(reliable, winding, -32768), axis=1)
        parity = counts % 2
        consistent = (low == high) & np.all(~reliable | (parity == (winding % 2)), axis=1)
        classification = np.full(n, -1, np.int8)
        nonparallel = np.zeros(n, bool)
        for a in range(k):
            for b in range(a+1, k):
                if abs(float(d[a] @ d[b])) < 1-1e-6:
                    nonparallel |= reliable[:, a] & reliable[:, b]
        usable = (votes >= min_reliable_rays) & nonparallel & consistent & ~near_surface
        if self.topology['nonmanifold_edges'] or self.topology['inconsistent_shared_edge_winding']:
            usable[:] = False
        classification[usable & (low == 0)] = 0
        classification[usable & (low == 1)] = 1
        arrays = dict(classification=classification, ray_winding=winding, ray_crossing_count=counts,
            ray_parity=parity.astype(np.int8), ray_reliable=reliable, ray_boundary_shadow=shadow,
            ray_edge_or_grazing=edge, ray_numerical_uncertain=numerical,
            reliable_rays=votes.astype(np.int8), ray_disagreement=(low != high) & (votes > 1),
            near_surface=near_surface, nonparallel_clean_rays=nonparallel, surface_distance_upper_bound_m=nearest['distance_m'],
            nearest_source_face_id=nearest['face_ids'], query_rounding_m=query_rounding)
        proof = dict(schema='oriented-all-intersection-diagnostic/v1', **self.index.evidence(),
            topology=self.topology, queries=n, directions=d.tolist(), native_proposed_events=native_hits,
            effective_uncertainty_networks=self.uncertainty_network_summary,
            boundary_evidence=self.boundary_evidence,
            rejected_literal_events=invalid_hits, event_statistics=events,
            inside_consensus=int((classification == 1).sum()), outside_consensus=int((classification == 0).sum()),
            uncertain=int((classification == -1).sum()), disagreement=int(arrays['ray_disagreement'].sum()),
            near_surface=int(near_surface.sum()), any_boundary_shadow=int(shadow.any(1).sum()),
            maximum_query_rounding_m=float(query_rounding.max(initial=0)), boundary_projection_epsilon_m=eps_m,
            min_reliable_rays=min_reliable_rays, near_surface_m=near_surface_m,
            all_native_events_replayed_on_literal_float64_triangles=True,
            all_hit_completeness_certified=False, source_open_surface_not_declared_closed=bool(len(self.edges)),
            boundary_policy='Each connected original boundary network has its own conservative projected AABB shadow plus numerical tolerance; branched networks are localized. No cap is constructed.',
            two_nonparallel_clean_rays_required=True,
            classification_scope='Multi-ray oriented winding/parity consensus; uncertainty is retained. Not anatomical acceptance.',
            source_modified=False, elapsed_seconds=time.perf_counter()-began)
        return arrays, proof

    def replay_all_triangles(self, points, directions, arrays, *, indices):
        """Bounded independent Float64 all-triangle intersection completeness check."""
        p = _points(points); d = _points(directions); d /= np.linalg.norm(d, axis=1)[:, None]
        mismatches = []; compared = 0; edge_cases = 0
        tri = self.mesh.triangles; count = len(tri)
        for row in np.asarray(indices, int):
            for j, direction in enumerate(d):
                hit, _, bary, signs, cosine = literal_hits(np.broadcast_to(p[row], (count, 3)),
                    np.broadcast_to(direction, (count, 3)), tri)
                if np.any(hit & ((bary.min(1) <= 1e-5) | (np.abs(cosine) <= 1e-6))):
                    edge_cases += 1; continue
                compared += 1; c = int(hit.sum()); w = int(signs[hit].sum())
                if c != arrays['ray_crossing_count'][row, j] or w != arrays['ray_winding'][row, j]:
                    mismatches.append(dict(query_index=int(row), ray=j, literal_count=c, literal_winding=w,
                        native_replayed_count=int(arrays['ray_crossing_count'][row, j]),
                        native_replayed_winding=int(arrays['ray_winding'][row, j])))
        return dict(sample_queries=len(indices), compared_clean_rays=compared,
            edge_or_grazing_rays_kept_uncertain=edge_cases, mismatches=mismatches,
            status='pass-bounded-completeness-replay' if not mismatches else 'fail-completeness-replay',
            entire_mesh_all_hit_completeness_proven=False)
