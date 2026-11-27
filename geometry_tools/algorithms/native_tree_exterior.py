"""Source-bound outer casts for explicitly connected native tube components.

``construct`` preserves native exterior mouths and source corner fields.  It
unions closed construction cells, using caller-specified conical continuations
at internal mouths, then removes every temporary exterior classifier disk.
New continuation faces retain separate constructed ancestry; they are never
reported as original source triangles.

This constructs a macroscopic exterior, not a finite vessel wall, a lumen, a
perfusion territory or an accepted anatomical interface.  The caller must bind
the mouth/anchor anatomy, inspect actual surfaces and sections, audit fixed
context, and separately verify Float32 storage before any local installation.
Use a separately supervised worker as well as the native 50-second deadline.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from contextlib import contextmanager
import signal
import threading
import time

import manifold3d as md
import numpy as np

from geometry_tools.algorithms.finite_tissue_boolean import (
    barycentric,
    make_solid,
    provenance,
)
from geometry_tools.audit import topology, winding_contains
from geometry_tools.contact import self_contact_summary
from geometry_tools.model import Mesh


class NativeTreeExteriorError(ValueError):
    """The supplied construction or its exact source bindings are invalid."""


def _require(condition, message):
    if not condition:
        raise NativeTreeExteriorError(message)


@contextmanager
def _wall_clock_limit(seconds):
    _require(0 < seconds <= 50, "A native deadline of at most 50 seconds is required")
    _require(threading.current_thread() is threading.main_thread(),
             "Run this bounded constructor in the main thread of a supervised worker")
    started = time.monotonic()
    previous_handler = signal.getsignal(signal.SIGALRM)
    previous_delay, previous_interval = signal.getitimer(signal.ITIMER_REAL)

    def expired(*_):
        raise TimeoutError("Native tree exterior construction exceeded its wall-clock deadline")

    signal.signal(signal.SIGALRM, expired)
    signal.setitimer(signal.ITIMER_REAL,
                     min(seconds, previous_delay) if previous_delay > 0 else seconds)
    try:
        yield
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)
        if previous_delay > 0:
            remaining = previous_delay - (time.monotonic() - started)
            signal.setitimer(signal.ITIMER_REAL, max(remaining, 1e-9), previous_interval)


def _boundary(faces):
    edges = np.vstack([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    _, first, counts = np.unique(np.sort(edges, axis=1), axis=0,
                                 return_index=True, return_counts=True)
    return edges[first[counts == 1]]


def _loops(faces):
    edges = _boundary(faces)
    following = {int(a): int(b) for a, b in edges}
    _require(len(following) == len(edges)
             and len({int(b) for _, b in edges}) == len(edges),
             "A boundary must be a disjoint union of directed cycles")
    result = []
    while following:
        start = min(following)
        cycle = [start]
        current = following.pop(start)
        while current != start:
            _require(current in following, "Open or repeated boundary vertex")
            cycle.append(current)
            current = following.pop(current)
        _require(len(cycle) >= 3, "A mouth requires at least three controls")
        result.append(np.asarray(cycle, dtype=np.int64))
    return result


def _signature(points):
    return frozenset(map(tuple, points))


def _oriented(triangle):
    return min(tuple(map(tuple, np.roll(triangle, k, axis=0))) for k in range(3))


def _quantize(normals):
    lengths = np.linalg.norm(normals, axis=-1)
    _require(np.isfinite(normals).all() and np.all(lengths > 0),
             "Zero or nonfinite constructed/interpolated normal")
    return np.rint(np.clip(normals / lengths[..., None], -1, 1) * 32767).astype(np.int16)


def _cone_fields(points, tip, faces, *, generated_normal_mode='mean'):
    _require(generated_normal_mode in ('mean', 'face'),
             "Generated continuation normal mode must be mean or face")
    triangles = np.vstack([points, tip])[faces]
    face_normals = np.cross(triangles[:, 1] - triangles[:, 0],
                            triangles[:, 2] - triangles[:, 0])
    count = len(points)
    rim = np.zeros((count, 3))
    for index, face in enumerate(faces):
        for vertex in face:
            if vertex < count:
                rim[vertex] += face_normals[index]
    lengths = np.linalg.norm(rim, axis=1)
    apex = face_normals.sum(0)
    _require(np.all(lengths > 0) and np.linalg.norm(apex) > 0,
             "Degenerate continuation normal field")
    rim /= lengths[:, None]
    apex /= np.linalg.norm(apex)
    normals = np.vstack([rim, apex])
    distances = np.linalg.norm(np.roll(points, -1, axis=0) - points, axis=1)
    _require(np.all(distances > 0), "Repeated continuation rim control")
    parameters = np.r_[0, np.cumsum(distances[:-1])] / distances.sum()
    uvs = []
    for a, b, _ in faces:
        ua, ub = float(parameters[a]), float(parameters[b])
        if abs(ua - ub) > .5:
            if ua < ub:
                ua += 1
            else:
                ub += 1
        uvs.append([[ua, 0.], [ub, 0.], [(ua + ub) / 2, 1.]])
    corner_normals = normals[faces] if generated_normal_mode == 'mean' else np.broadcast_to(
        face_normals[:, None, :], (len(faces), 3, 3))
    return _quantize(corner_normals), np.asarray(uvs, dtype=np.float32)


def construct(mesh: Mesh, component_face_ids: Sequence | Mapping,
              mouths: Sequence[Mapping], anchors: Mapping, *,
              wall_seconds: float = 50., max_pairs: int = 1_000_000,
              max_tests: int = 3_000_000, generated_normal_mode: str = 'mean',
              source_corner_normal_words=None, source_corner_uvs=None,
              expected_genus: int = 0):
    """Return ``(arrays, proof)`` without writing files or modifying ``mesh``.

    ``component_face_ids`` partitions *every* input face exactly once.  A mapping
    preserves its insertion order and records its keys; a sequence uses numeric
    component keys.  ``mouths`` gives every actual component boundary once, with
    ``positions_m``, integer ``source_control_ids`` into ``mesh.vertices``, and
    ``kind``.  Unanchored exterior mouths must have a kind beginning ``original``.
    ``anchors`` maps each internal mouth's sequence index to an explicit 3D tip.
    No missing mouth, anatomical branch identity or anchor is inferred.
    ``generated_normal_mode='face'`` assigns outward constant normals only to
    new conical primitive charts. Native fields, UVs and geometry are unchanged.
    Explicit source corner fields preserve distinct property charts at one
    physical control; they override the corresponding vertex-field lookup.
    A nonzero expected_genus requires the caller's explicit anatomical/source
    contract; the constructor never infers permission for an extra tunnel.

    Hard source/field/provenance violations raise ``NativeTreeExteriorError``.
    Completed geometry audits return an explicit ``proof['passed']``; callers
    must reject false or incomplete results.  The deadline raises ``TimeoutError``.
    Full context, wall/lumen, anatomical and storage acceptance stay false.
    """
    started = time.monotonic()
    with _wall_clock_limit(float(wall_seconds)):
        return _construct(mesh, component_face_ids, mouths, anchors,
                          max_pairs=max_pairs, max_tests=max_tests, started=started,
                          generated_normal_mode=generated_normal_mode,
                          source_corner_normal_words=source_corner_normal_words,
                          source_corner_uvs=source_corner_uvs, expected_genus=expected_genus)


def _construct(mesh, component_face_ids, mouths, anchors, *, max_pairs, max_tests, started,
               generated_normal_mode, source_corner_normal_words=None, source_corner_uvs=None,
               expected_genus=0):
    _require(isinstance(mesh, Mesh) and len(mesh.faces) > 0, "A nonempty source Mesh is required")
    _require(isinstance(anchors, Mapping), "Explicit mouth-index to point anchors are required")
    _require(generated_normal_mode in ('mean', 'face'),
             "Generated continuation normal mode must be mean or face")
    _require(isinstance(expected_genus, int) and not isinstance(expected_genus, bool)
             and 0 <= expected_genus <= 32, "Explicit expected genus must be an integer from zero to 32")
    normal_source = mesh.attributes.get('normals')
    uv_source = mesh.attributes.get('uvs')
    _require(normal_source is not None and normal_source.shape == (len(mesh.vertices), 3)
             and normal_source.dtype == np.dtype('int16'), "Literal int16 source normal words are required")
    _require(uv_source is not None and uv_source.shape == (len(mesh.vertices), 2)
             and uv_source.dtype == np.dtype('float32') and np.isfinite(uv_source).all(),
             "Literal finite float32 source UV fields are required")
    normal_charts = (normal_source[mesh.faces] if source_corner_normal_words is None
                     else np.asarray(source_corner_normal_words))
    uv_charts = (uv_source[mesh.faces] if source_corner_uvs is None
                 else np.asarray(source_corner_uvs))
    _require(normal_charts.shape == (len(mesh.faces), 3, 3)
             and normal_charts.dtype == np.dtype('int16')
             and np.all(np.linalg.norm(normal_charts, axis=2) > 0),
             "One nonzero int16 source normal word is required per face corner")
    _require(uv_charts.shape == (len(mesh.faces), 3, 2)
             and uv_charts.dtype == np.dtype('float32') and np.isfinite(uv_charts).all(),
             "One finite float32 source UV is required per face corner")
    if isinstance(component_face_ids, Mapping):
        component_keys = list(component_face_ids)
        components = list(component_face_ids.values())
    else:
        components = list(component_face_ids)
        component_keys = list(range(len(components)))
    _require(bool(components), "Explicit source components are required")
    components = [np.asarray(ids) for ids in components]
    for ids in components:
        _require(ids.ndim == 1 and len(ids) > 0 and np.issubdtype(ids.dtype, np.integer),
                 "Each component must contain actual integer source face IDs")
    all_faces = np.concatenate(components)
    _require(len(all_faces) == len(mesh.faces)
             and np.array_equal(np.sort(all_faces), np.arange(len(mesh.faces))),
             "Components must partition every source face exactly once")
    mouths = list(mouths)
    mouth_positions, mouth_controls, mouth_lookup = [], [], {}
    for index, record in enumerate(mouths):
        points = np.asarray(record['positions_m'], dtype=float)
        controls = np.asarray(record['source_control_ids'])
        _require(points.ndim == 2 and points.shape[1:] == (3,) and len(points) >= 3
                 and np.isfinite(points).all(), "Finite actual mouth positions are required")
        _require(controls.shape == (len(points),) and np.issubdtype(controls.dtype, np.integer)
                 and controls.min() >= 0 and controls.max() < len(mesh.vertices)
                 and len(set(controls.tolist())) == len(controls), "Actual distinct source mouth control IDs are required")
        _require(np.array_equal(mesh.vertices[controls], points), "Mouth positions do not match literal source controls")
        key = _signature(points)
        _require(len(key) == len(points) and key not in mouth_lookup, "A mouth is repeated or has repeated positions")
        mouth_lookup[key] = index
        mouth_positions.append(points)
        mouth_controls.append(controls.astype(np.int64))
    actual_anchors = {}
    for index, point in anchors.items():
        _require(isinstance(index, (int, np.integer)) and 0 <= index < len(mouths),
                 "Anchor key must be an actual mouth index")
        point = np.asarray(point, dtype=float)
        _require(point.shape == (3,) and np.isfinite(point).all(), "Finite explicit 3D anchor required")
        _require(not str(mouths[index].get('kind', '')).startswith('original'),
                 "An original exterior mouth cannot be replaced by a continuation")
        actual_anchors[int(index)] = point
    for index, record in enumerate(mouths):
        _require(index in actual_anchors or str(record.get('kind', '')).startswith('original'),
                 "Every internal mouth needs an explicit anchor")

    solids, primitives = [], {}
    derived_triangles, derived_normals, derived_uvs = [], [], []
    derived_components, derived_mouths = [], []
    seen_mouths = []
    for component, ids in enumerate(components):
        used = np.unique(mesh.faces[ids])
        mapping = np.full(len(mesh.vertices), -1, dtype=np.int64)
        mapping[used] = np.arange(len(used))
        vertices, faces = mesh.vertices[used].copy(), mapping[mesh.faces[ids]].copy()
        source_parent = ids.copy()
        derived_parent = np.full(len(faces), -1, dtype=np.int64)
        classifier = np.zeros(len(faces), dtype=np.uint8)
        for cycle in _loops(faces):
            key = _signature(vertices[cycle])
            _require(key in mouth_lookup, "A component boundary has no explicitly supplied mouth")
            mouth_index = mouth_lookup[key]
            seen_mouths.append(mouth_index)
            internal = mouth_index in actual_anchors
            tip = actual_anchors[mouth_index] if internal else vertices[cycle].mean(0)
            tip_id, ring = len(vertices), vertices[cycle].copy()
            vertices = np.vstack([vertices, tip])
            new_faces = np.asarray([[cycle[(j + 1) % len(cycle)], cycle[j], tip_id]
                                    for j in range(len(cycle))], dtype=np.int64)
            new_parent = np.full(len(new_faces), -1, dtype=np.int64)
            if internal:
                local_faces = np.asarray([[(j + 1) % len(cycle), j, len(cycle)]
                                          for j in range(len(cycle))], dtype=np.int64)
                normals, uvs = _cone_fields(ring, tip, local_faces,
                                             generated_normal_mode=generated_normal_mode)
                first = len(derived_triangles)
                derived_triangles.extend(vertices[new_faces])
                derived_normals.extend(normals)
                derived_uvs.extend(uvs)
                derived_components.extend([component] * len(new_faces))
                derived_mouths.extend([mouth_index] * len(new_faces))
                new_parent = np.arange(first, first + len(new_faces), dtype=np.int64)
            faces = np.vstack([faces, new_faces])
            source_parent = np.r_[source_parent, np.full(len(new_faces), -1, dtype=np.int64)]
            derived_parent = np.r_[derived_parent, new_parent]
            classifier = np.r_[classifier, np.full(len(new_faces), 0 if internal else 1, dtype=np.uint8)]
        centered = vertices - vertices.mean(0)
        volume = float(np.einsum('ij,ij->i', centered[faces[:, 0]],
                                  np.cross(centered[faces[:, 1]], centered[faces[:, 2]])).sum() / 6)
        _require(volume > 0, "Nonpositive closed source/continuation cell volume")
        solid, original_id = make_solid(vertices, faces)
        solids.append(solid)
        primitives[original_id] = dict(source_parent=source_parent, derived_parent=derived_parent,
                                       classifier=classifier, component=component)
    _require(sorted(seen_mouths) == list(range(len(mouths))), "Mouths must cover the exact component boundaries once")

    dv = np.asarray(derived_triangles, dtype=float).reshape(-1, 3)
    df = np.arange(len(dv), dtype=np.int64).reshape(-1, 3)
    dn = np.asarray(derived_normals, dtype=np.int16).reshape(-1, 3, 3)
    duv = np.asarray(derived_uvs, dtype=np.float32).reshape(-1, 3, 2)
    dcomponents = np.asarray(derived_components, dtype=np.int64)
    dmouths = np.asarray(derived_mouths, dtype=np.int64)
    union = md.Manifold.batch_boolean(solids, md.OpType.Add)
    _require(union.status() == md.Error.NoError, f"Boolean status: {union.status()}")
    packed = union.to_mesh64()
    vertices, faces = np.asarray(packed.vert_properties)[:, :3], np.asarray(packed.tri_verts, dtype=np.int64)
    original, parent = provenance(packed)
    source_parent = np.asarray([primitives[int(a)]['source_parent'][b] for a, b in zip(original, parent)])
    derived_parent = np.asarray([primitives[int(a)]['derived_parent'][b] for a, b in zip(original, parent)])
    classifiers = np.asarray([primitives[int(a)]['classifier'][b] for a, b in zip(original, parent)])
    component_ids = np.asarray([primitives[int(a)]['component'] for a in original])
    closed_face_count = len(faces)
    if actual_anchors:
        internal_centers = np.asarray([mouth_positions[index].mean(0) for index in actual_anchors])
        tips = np.asarray(list(actual_anchors.values()))
        probes = np.vstack([internal_centers, (internal_centers + tips) / 2, tips])
        communication = winding_contains(Mesh(vertices, faces, mesh.frame), probes)
        positive = bool(all(communication['inside']) and not any(communication['boundary_ambiguous']))
    else:
        probes = np.empty((0, 3))
        communication = {'inside': [], 'boundary_ambiguous': [], 'qualification': 'No internal continuation probes requested'}
        positive = True
    keep = classifiers == 0
    faces, source_parent = faces[keep], source_parent[keep]
    derived_parent, component_ids = derived_parent[keep], component_ids[keep]
    _require(np.all((source_parent >= 0) ^ (derived_parent >= 0)),
             "Unclassified or double-parent output surface")
    used, inverse = np.unique(faces, return_inverse=True)
    vertices, faces = vertices[used].copy(), inverse.reshape(-1, 3)
    positions = {tuple(point): index for index, point in enumerate(vertices)}
    exterior_indices = [i for i in range(len(mouths)) if i not in actual_anchors]
    ports = []
    for index in exterior_indices:
        points = mouth_positions[index]
        _require(all(tuple(point) in positions for point in points), "Original exterior mouth moved or disappeared")
        ports.append(np.asarray([positions[tuple(point)] for point in points], dtype=np.int64))
    expected_edges = {tuple(sorted((int(cycle[j]), int(cycle[(j + 1) % len(cycle)]))))
                      for cycle in ports for j in range(len(cycle))}
    actual_edges = {tuple(sorted(map(int, edge))) for edge in _boundary(faces)}
    _require(actual_edges == expected_edges, "Original exterior mouth edges were not retained exactly")

    source_barycentric = np.zeros((len(faces), 3, 3))
    derived_barycentric = source_barycentric.copy()
    normal_words = np.empty((len(faces), 3, 3), dtype=np.int16)
    uvs = np.empty((len(faces), 3, 2), dtype=np.float32)
    native = np.zeros((len(faces), 3), dtype=bool)
    replay, inverted, bad_normals = 0., 0, 0
    parent_faces = {}
    kinds = (derived_parent >= 0).astype(np.uint8)
    for index, face in enumerate(faces):
        derived = bool(kinds[index])
        chart = dv[df[derived_parent[index]]] if derived else mesh.triangles[source_parent[index]]
        points = vertices[face]
        weights = barycentric(points, chart)
        error = float(np.linalg.norm(weights @ chart - points, axis=1).max())
        replay = max(replay, error)
        _require(error <= 2e-12 and weights.min() >= -1e-8 and weights.max() <= 1 + 1e-8,
                 "Output left its actual source or explicitly constructed parent chart")
        face_normal = np.cross(points[1] - points[0], points[2] - points[0])
        parent_normal = np.cross(chart[1] - chart[0], chart[2] - chart[0])
        inverted += int(face_normal @ parent_normal <= 0)
        if derived:
            derived_barycentric[index] = weights
            normal_words[index] = _quantize(weights @ dn[derived_parent[index]].astype(float))
            uvs[index] = (weights @ duv[derived_parent[index]]).astype(np.float32)
        else:
            source_barycentric[index] = weights
            controls = mesh.faces[source_parent[index]]
            normal_words[index] = _quantize(weights @ normal_charts[source_parent[index]].astype(float))
            uvs[index] = (weights @ uv_charts[source_parent[index]]).astype(np.float32)
            for corner, point in enumerate(points):
                match = np.flatnonzero(np.all(chart == point, axis=1))
                if len(match):
                    native[index, corner] = True
                    normal_words[index, corner] = normal_charts[source_parent[index], match[0]]
                    uvs[index, corner] = uv_charts[source_parent[index], match[0]]
            parent_faces.setdefault(int(source_parent[index]), []).append(index)
        bad_normals += int(np.sum(normal_words[index] @ face_normal <= 0))
    unchanged, changed = [], []
    for parent in range(len(mesh.faces)):
        ids = parent_faces.get(parent, [])
        target = unchanged if len(ids) == 1 and _oriented(vertices[faces[ids[0]]]) == _oriented(mesh.triangles[parent]) else changed
        target.append(parent)
    result = Mesh(vertices, faces, mesh.frame)
    top = topology(result)
    self_contacts, _ = self_contact_summary(result, max_pairs=max_pairs, max_tests=max_tests, example_limit=10000)
    failures = [key for key in ['nonmanifold_edges', 'inconsistent_winding_edges', 'invalid_vertex_links',
                                'duplicate_faces', 'zero_area_faces'] if top[key]]
    expected_euler = 2 - 2*expected_genus - len(ports)
    if len(top['components']) != 1 or top['components'][0]['euler'] != expected_euler:
        failures.append('wrong_declared_genus_topology')
    if len(_loops(faces)) != len(ports):
        failures.append('wrong_exterior_mouth_count')
    if self_contacts['defect_pairs'] or not self_contacts['complete_within_pair_budget']:
        failures.append('self_crossing_or_incomplete_budget')
    if inverted:
        failures.append('parent_fragment_orientation')
    if bad_normals:
        failures.append('nonforward_interpolated_corner_normals')
    if not positive:
        failures.append('positive_connected_cast_interior_probe')
    arrays = dict(vertices=vertices, faces=faces, frame=np.asarray(mesh.frame),
                  source_vertices=mesh.vertices, source_faces=mesh.faces,
                  source_normal_words=normal_source, source_uvs=uv_source,
                  source_corner_normal_words=normal_charts, source_corner_uvs=uv_charts,
                  face_parent_kind=kinds, source_parent_face_ids=source_parent,
                  source_parent_face_id=source_parent, source_component_ids=component_ids,
                  source_component_keys=np.asarray(list(map(str, component_keys))),
                  corner_source_barycentric=source_barycentric, source_parent_barycentric=source_barycentric,
                  derived_parent_face_ids=derived_parent, corner_derived_barycentric=derived_barycentric,
                  derived_vertices=dv, derived_faces=df, derived_component_ids=dcomponents,
                  derived_mouth_ids=dmouths, derived_corner_normal_words=dn, derived_corner_uvs=duv,
                  corner_normal_words=normal_words, face_corner_normal_words=normal_words,
                  corner_uvs=uvs, face_corner_uvs=uvs, native_corner=native,
                  unchanged_source_parent_face_ids=np.asarray(unchanged, dtype=np.int64),
                  changed_source_parent_face_ids=np.asarray(changed, dtype=np.int64),
                  communication_probe_points_m=probes)
    for number, (mouth_index, port) in enumerate(zip(exterior_indices, ports)):
        arrays[f'external_boundary_{number}_vertex_ids'] = port
        arrays[f'external_boundary_{number}_positions_m'] = mouth_positions[mouth_index]
        arrays[f'external_boundary_{number}_source_control_ids'] = mouth_controls[mouth_index]
    proof = dict(schema='source-bound-native-tree-exterior/v1', complete=True, passed=not failures,
                 source_geometry_sha256=mesh.geometry_sha, source_attributes_sha256=mesh.attributes_sha,
                 input_components=len(components), input_mouths=len(mouths), internal_mouths=len(actual_anchors),
                 expected_exterior_mouths=len(ports), expected_boundary_edges=len(expected_edges),
                 expected_euler=expected_euler, vertices=len(vertices), faces=len(faces),
                 expected_genus=expected_genus,
                 closed_construction_faces=closed_face_count, removed_classifier_fragments=int((~keep).sum()),
                 temporary_classifier_faces_remaining=0, all_exterior_profiles_exact=True,
                 source_course_motion=False, source_calibre_change=False,
                 native_fragments=int((kinds == 0).sum()), derived_fragments=int(kinds.sum()),
                 derived_primitive_faces=len(df), source_original_faces_unchanged=len(unchanged),
                 source_original_faces_changed=len(changed), maximum_parent_replay_error_m=replay,
                 source_corner_words_literal_at_retained_controls=True,
                 generated_normal_mode=generated_normal_mode,
                 inverted_fragments=inverted, nonforward_corner_normals=bad_normals,
                 topology=top, self_intersections=self_contacts,
                 positive_cast_interior_probes=communication, failures=failures,
                 full_context_passed=False, anatomical_acceptance=False,
                 finite_wall_or_physiological_patency_claim=False, storage_acceptance=False,
                 installed=False, elapsed_seconds=time.monotonic() - started)
    return arrays, proof
