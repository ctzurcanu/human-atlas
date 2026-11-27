"""Adopt complementary PLC cells without a second coincident-surface Boolean.

The caller supplies anatomical ownership. A replacement must reproduce the
complete literal source boundary, and every new internal triangle must have
exactly two opposed owners. Source normal/UV charts are replayed on fragments.
"""
from collections import defaultdict
import time
import numpy as np
import igl
from .finite_tissue_boolean import barycentric, make_solid
from .conforming_fields import replay_corner_fields
from .source_face_complement import conform_surface_boundary_edges


def adopt_partition(positions, surfaces, cells, normals, uvs, *, host,
                    partition_positions, partition_cells, original_control_ids,
                    wall_seconds=90.):
    start = time.monotonic()
    def deadline():
        if time.monotonic()-start > wall_seconds:
            raise TimeoutError('Direct shared-tissue adoption deadline')
    p = np.asarray(positions, float)
    v = np.asarray(partition_positions, float)
    ids = np.asarray(original_control_ids, int)
    if host not in cells or not 1 < len(partition_cells) <= 8 or host not in partition_cells:
        raise ValueError('Explicit existing host and two to eight replacement owners required')
    if any(name in cells and name != host for name in partition_cells):
        raise ValueError('A replacement may not overwrite an unrelated tissue owner')
    if not np.array_equal(v[:len(ids)], p[ids]):
        raise ValueError('PLC carrier controls differ from the literal source')
    q = np.r_[p, v[len(ids):]]
    mapping = np.r_[ids, np.arange(len(v)-len(ids))+len(p)]
    child = {name: mapping[np.asarray(f, int)] for name, f in partition_cells.items()}
    all_faces = np.concatenate(list(child.values()))
    keys = np.sort(all_faces, axis=1)
    _, inverse, counts = np.unique(keys, axis=0, return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError('Multiply duplicated partition faces')
    pieces = np.split(inverse, np.cumsum([len(f) for f in child.values()])[:-1])
    internal = {}
    boundary = {}
    for (name, f), index in zip(child.items(), pieces):
        boundary[name] = f[counts[index] == 1]
        for face, key in zip(f[counts[index] == 2], index[counts[index] == 2]):
            internal.setdefault(int(key), []).append((name, face))
    uses = cells[host]
    host_faces = np.concatenate([surfaces[u['surface']] if u['side'] == 1 else
                                 surfaces[u['surface']][:, [0, 2, 1]] for u in uses])
    parent_rows = [(u['surface'], i, u['side']) for u in uses
                   for i in range(len(surfaces[u['surface']]))]
    literal_faces = {tuple(sorted(map(int, f))): i for i, f in enumerate(host_faces)}
    source_support = defaultdict(set)
    for ordinal, face in enumerate(host_faces):
        for vertex in face:
            source_support[int(vertex)].add(ordinal)
    split = defaultdict(list)
    parent_ids = defaultdict(list)
    area_replay = np.zeros(len(host_faces))
    maximum_plane_error = 0.
    for name, f in boundary.items():
        deadline()
        tri = q[f]
        d, parent, _ = igl.point_mesh_squared_distance(tri.mean(1), p, host_faces)
        literal = np.array([literal_faces.get(tuple(sorted(map(int, face))), -1) for face in f])
        direct = literal >= 0
        parent[direct] = literal[direct]
        # Thin, unchanged source triangles can defeat a nearest-point
        # numerical predicate. Literal indexed identity is the stronger proof.
        d[direct] = 0
        pending = list(zip(f, parent))
        resolved = []
        for propagation in range(4):
            unresolved = []
            for face, proposed in pending:
                candidates = [int(proposed)]
                candidates.extend(sorted(set.union(*(source_support[int(v)] for v in face))-set(candidates)))
                found = None
                for ordinal in candidates:
                    source_face_ids = host_faces[int(ordinal)]
                    if tuple(sorted(map(int, face))) == tuple(sorted(map(int, source_face_ids))):
                        weights = np.array([[float(v == i) for i in source_face_ids] for v in face])
                    else:
                        weights = barycentric(q[face], p[source_face_ids])
                    error = np.max(np.linalg.norm(weights @ p[source_face_ids]-q[face], axis=1))
                    if error <= 2e-12 and weights.min() >= -2e-8 and weights.max() <= 1+2e-8:
                        found = ordinal, weights, error
                        break
                if found is None:
                    unresolved.append((face, proposed))
                    continue
                ordinal, weights, error = found
                resolved.append((face, ordinal, weights, error))
                for vertex in face:
                    source_support[int(vertex)].add(int(ordinal))
            if not unresolved:
                break
            if len(unresolved) == len(pending):
                raise ValueError('A PLC boundary fragment has no contained literal carrier chart')
            pending = unresolved
        else:
            raise ValueError('PLC boundary-parent propagation exceeded its finite budget')
        for face, ordinal, weights, error in resolved:
            sid, source_face, sign = parent_rows[int(ordinal)]
            maximum_plane_error = max(maximum_plane_error, float(error))
            split[sid, name].append(face if sign == 1 else face[[0, 2, 1]])
            parent_ids[sid, name].append(source_face)
            t = q[face]
            area_replay[ordinal] += np.linalg.norm(np.cross(t[1]-t[0], t[2]-t[0]))/2
    original_tri = p[host_faces]
    original_area = np.linalg.norm(np.cross(original_tri[:, 1]-original_tri[:, 0],
                                            original_tri[:, 2]-original_tri[:, 0]), axis=1)/2
    if np.max(abs(area_replay-original_area)) > 1e-14:
        raise ValueError('PLC boundary fragments do not exactly cover every original facet')
    owners = defaultdict(list)
    for name, refs in cells.items():
        for ref in refs:
            owners[ref['surface']].append((name, ref))
    result = {sid: f for sid, f in surfaces.items()
              if not any(pid == host for pid, _ in owners[sid])}
    output_cells = {name: [dict(ref) for ref in refs
                          if ref['surface'] in result] for name, refs in cells.items()}
    for name in partition_cells:
        output_cells.setdefault(name, [])
    words = {sid: normals[sid] for sid in result}
    tex = {sid: uvs[sid] for sid in result}
    for (sid, name), faces in split.items():
        deadline()
        new_id = sid+'/owned:'+name
        fs = np.asarray(faces, int)
        parent = np.asarray(parent_ids[sid, name], int)
        result[new_id] = fs
        pseudo = dict(positions=q, surfaces={new_id: fs}, charts={new_id: dict(
            source_surfaces=[sid]*len(fs), source_faces=parent,
            normal_signs=np.ones(len(fs), int), barycentrics=np.array([
                barycentric(q[face], p[surfaces[sid][int(i)]]) for face, i in zip(fs, parent)]))})
        ww, tt = replay_corner_fields(p, surfaces, normals, uvs, pseudo)
        words[new_id], tex[new_id] = ww[new_id], tt[new_id]
        old_use = next(ref for pid, ref in owners[sid] if pid == host)
        peers = [(pid, ref) for pid, ref in owners[sid] if pid != host]
        output_cells[name].append(dict(surface=new_id, side=old_use['side'],
                                       **({'neighbor': peers[0][0]} if peers else {})))
        for peer, ref in peers:
            output_cells[peer].append(dict(surface=new_id, side=ref['side'], neighbor=name))
    interface_groups = defaultdict(list)
    for pair in internal.values():
        if len(pair) != 2 or pair[0][0] == pair[1][0]:
            raise ValueError('An internal face has invalid tissue owners')
        a, b = sorted(pair, key=lambda r: r[0])
        if not any(np.array_equal(np.roll(a[1], k), b[1][[0, 2, 1]]) for k in range(3)):
            raise ValueError('Internal tissue faces are not literally opposed')
        interface_groups[a[0], b[0]].append(a[1])
    from geometry_tools.tissue import normal_words
    for (a, b), faces in interface_groups.items():
        sid = 'owned-interface:'+a+':'+b
        fs = np.asarray(faces, int)
        result[sid] = fs
        words[sid] = normal_words(q, {sid: fs})[sid]
        tex[sid] = np.zeros((len(fs), 3, 2), '<f4')
        output_cells[a].append(dict(surface=sid, side=1, neighbor=b))
        output_cells[b].append(dict(surface=sid, side=-1, neighbor=a))
    deadline()
    conformed, conformity = conform_surface_boundary_edges(q, result, output_cells,
                                      tolerance=2e-12, record_face_parents=True)
    ancestry = conformity.pop('sourceFaceParents')
    for sid, fs in conformed.items():
        if np.array_equal(fs, result[sid]):
            continue
        parent = ancestry[sid]
        pseudo = dict(positions=q, surfaces={sid: fs}, charts={sid: dict(
            source_surfaces=[sid]*len(fs), source_faces=parent,
            normal_signs=np.ones(len(fs), int), barycentrics=np.array([
                barycentric(q[face], q[result[sid][int(i)]]) for face, i in zip(fs, parent)]))})
        ww, tt = replay_corner_fields(q, {sid: result[sid]}, {sid: words[sid]}, {sid: tex[sid]}, pseudo)
        words[sid], tex[sid] = ww[sid], tt[sid]
    deadline()
    volumes, engine_volumes = {}, {}
    for name, refs in output_cells.items():
        fs = np.concatenate([conformed[r['surface']] if r['side'] == 1 else
                              conformed[r['surface']][:, [0, 2, 1]] for r in refs])
        solid, _ = make_solid(q, fs)
        tri = q[fs]
        volumes[name] = float(np.einsum('ij,ij->i', tri[:, 0],
                                  np.cross(tri[:, 1], tri[:, 2])).sum()/6)
        engine_volumes[name] = float(solid.volume())
        if volumes[name] <= 0:
            raise ValueError('A replacement has nonpositive occupied volume')
    source_tri = p[host_faces]
    source_volume = float(np.einsum('ij,ij->i', source_tri[:, 0],
                             np.cross(source_tri[:, 1], source_tri[:, 2])).sum()/6)
    error = sum(volumes[name] for name in partition_cells)-source_volume
    if abs(error) > max(1e-18, source_volume*1e-9):
        raise ValueError('Direct tissue adoption changes the original whole wall volume')
    return dict(positions=q, surfaces=conformed, cells=output_cells, normals=words, uvs=tex,
                evidence=dict(method='direct complementary PLC boundary adoption',
                  source_controls_literal=np.array_equal(q[:len(p)], p),
                  maximum_source_plane_error_m=maximum_plane_error,
                  maximum_source_face_area_error_m2=float(max(abs(area_replay-original_area))),
                  internal_interfaces_opposed=True, positive_occupied_volumes_m3=volumes,
                  boolean_engine_volumes_m3=engine_volumes,
                  volume_replay_uses_literal_input_triangles=True,
                  boolean_engine_simplification_is_not_emitted=True,
                  source_host_volume_m3=source_volume, replacement_volume_error_m3=error,
                  boundary_conformity=conformity, elapsed_seconds=time.monotonic()-start,
                  anatomical_acceptance=False))
