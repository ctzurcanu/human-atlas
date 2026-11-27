"""Partition selected finite cells, carrying both sides of every source chart.

The paired Boolean supplies the removed source fragments directly. Protected
neighbours can therefore become exposed without a second planar-complement
triangulation. A region name retains/reassigns the selected intersection;
without one the selected intersection is a computational void. No attachment
or anatomical boundary is inferred by this operator.
"""
from collections import defaultdict
import time
import numpy as np
from .finite_tissue_boolean import CommonPositionBank, make_solid, provenance, barycentric, same_oriented_triangles
from .shared_triangle_edges import reconcile_boundaries, refine_surfaces
from .conforming_fields import replay_corner_fields
from .source_face_complement import conform_surface_boundary_edges


def partition_mask(positions, surfaces, cells, normals, uvs, selected,
                   mask_positions, mask_faces, *, region=None, cap_faces=(),
                   tolerance=2e-12, wall_seconds=95., checkpoint=None, output_checkpoint=None,
                   mask_face_sources=None):
    started = time.monotonic()
    p, v = np.asarray(positions, float), np.asarray(mask_positions, float)
    cf = np.asarray(mask_faces, np.int64)
    if mask_face_sources is not None and len(mask_face_sources) != len(cf):
        raise ValueError('One explicit source-chart reference is required per mask face')
    selected = set(selected)
    if not selected or not selected <= set(cells) or region in cells:
        raise ValueError('Declare existing selected cells and a new optional region')
    def deadline():
        if time.monotonic() - started > wall_seconds:
            raise TimeoutError('Paired cell partition deadline')
    owners = defaultdict(list)
    for pid, uses in cells.items():
        for use in uses:
            owners[use['surface']].append((pid, use))
    bank = CommonPositionBank(np.r_[p, v], tolerance)
    cutter, cutter_id = make_solid(v, cf)
    fragments, cuts, volumes = {}, {}, {}
    source_boundary_replays = 0
    for pid in sorted(selected):
        deadline()
        uses = cells[pid]
        faces = np.concatenate([surfaces[u['surface']] if u['side'] == 1
                                else surfaces[u['surface']][:, [0, 2, 1]] for u in uses])
        parents = [(u['surface'], i, u['side']) for u in uses
                   for i in range(len(surfaces[u['surface']]))]
        host, host_id = make_solid(p, faces)
        inside, outside = host.split(cutter)
        volumes[pid] = dict(before=float(host.volume()), inside=float(inside.volume()),
                            outside=float(outside.volume()))
        if abs(volumes[pid]['before'] - volumes[pid]['inside'] - volumes[pid]['outside']) > max(1e-18, abs(host.volume()) * 1e-9):
            raise ValueError('Paired partition changed source volume: ' + pid)
        for side_name, solid in [('inside', inside), ('outside', outside)]:
            mm = solid.to_mesh64()
            origin, parent = provenance(mm)
            mapping = np.array([bank.insert(point) for point in mm.vert_properties[:, :3]], np.int64)
            ff = mapping[np.asarray(mm.tri_verts)]
            current_positions = bank.positions
            source_rows = defaultdict(list)
            source_parents = defaultdict(list)
            mask_rows, mask_parents = [], []
            for face, family, ordinal in zip(ff, origin, parent):
                if family == host_id:
                    sid, source_face, sign = parents[int(ordinal)]
                    source_rows[sid].append(face if sign == 1 else face[[0, 2, 1]])
                    source_parents[sid].append(source_face)
                elif family == cutter_id:
                    if int(ordinal) in cap_faces:
                        raise ValueError('A virtual end cap contributed to a real tissue cut')
                    ref = mask_face_sources[int(ordinal)] if mask_face_sources is not None else None
                    original_use = next((u for u in uses if ref and u['surface'] == ref['surface']), None)
                    if original_use is not None:
                        parent_triangle = p[surfaces[ref['surface']][int(ref['source_face'])]]
                        actual_triangle = current_positions[face]
                        expected_normal = np.cross(parent_triangle[1]-parent_triangle[0],
                            parent_triangle[2]-parent_triangle[0]) * original_use['side']
                        actual_normal = np.cross(actual_triangle[1]-actual_triangle[0],
                            actual_triangle[2]-actual_triangle[0])
                        if actual_normal @ expected_normal <= 0:
                            # This is the newly exposed side of the mask, not
                            # the original outward boundary of this host.
                            original_use = None
                    if original_use is not None:
                        # Coplanar Boolean boundaries can inherit the cutter's
                        # ancestry even when they are literally an existing
                        # shared tissue chart. Replay only an explicit matched
                        # parent; proximity alone never establishes identity.
                        sid, source_face = ref['surface'], int(ref['source_face'])
                        source_triangle = p[surfaces[sid][source_face]]
                        actual = current_positions[face]
                        weights = barycentric(actual, source_triangle)
                        if np.max(np.linalg.norm(weights @ source_triangle - actual, axis=1)) > tolerance*4:
                            raise ValueError('Declared coincident source chart is off its parent plane')
                        source_rows[sid].append(face if original_use['side'] == 1 else face[[0, 2, 1]])
                        source_parents[sid].append(source_face)
                        source_boundary_replays += 1
                    else:
                        mask_rows.append(face)
                        mask_parents.append(int(ordinal))
                else:
                    raise ValueError('Unknown paired Boolean ancestry')
            fragments[pid, side_name] = {sid: (np.asarray(fs, np.int64).reshape(-1, 3),
                np.asarray(source_parents[sid], np.int64)) for sid, fs in source_rows.items()}
            cuts[pid, side_name] = (np.asarray(mask_rows, np.int64).reshape(-1, 3),
                                    np.asarray(mask_parents, np.int64))
        if checkpoint:
            checkpoint(pid, volumes[pid])
    q = bank.positions
    out, charts, output_cells, plans = {}, {}, {pid: [] for pid in cells}, defaultdict(set)
    if region:
        output_cells[region] = []
    comparisons = []
    for sid, original in surfaces.items():
        deadline()
        changed = [pid for pid, _ in owners[sid] if pid in selected]
        if not changed:
            out[sid] = original
            charts[sid] = dict(source_surfaces=[sid] * len(original),
                source_faces=np.arange(len(original)), normal_signs=np.ones(len(original), int),
                barycentrics=np.broadcast_to(np.eye(3), (len(original), 3, 3)))
            for pid, use in owners[sid]:
                output_cells[pid].append(dict(use))
            continue
        for side_name in ['outside', 'inside']:
            fs, parent = fragments[changed[0], side_name].get(sid, (np.empty((0, 3), int), np.empty(0, int)))
            for other in changed[1:]:
                second, second_parent = fragments[other, side_name].get(sid, (np.empty((0, 3), int), np.empty(0, int)))
                if same_oriented_triangles(fs, second):
                    continue
                if set(parent) != set(second_parent):
                    raise ValueError('Reciprocal owners classified a source chart differently: ' + sid)
                def groups(values):
                    order = np.argsort(values, kind='stable')
                    keys, starts = np.unique(values[order], return_index=True)
                    return {int(key): ids for key, ids in zip(keys, np.split(order, starts[1:]))}
                first_groups, second_groups = groups(parent), groups(second_parent)
                for ordinal, ids in first_groups.items():
                    error = reconcile_boundaries(q, fs[ids],
                        second[second_groups[ordinal]], plans, tolerance=tolerance)
                    comparisons.append(error)
            use_map = {}
            for pid, use in owners[sid]:
                target = (pid if side_name == 'outside' else region) if pid in selected else pid
                if target:
                    use_map.setdefault(target, []).append(use['side'])
            # Internal boundaries of the reassigned intersection cancel.
            use_map = {pid: values[0] for pid, values in use_map.items() if sum(values)}
            if not len(fs) or not use_map:
                continue
            name = sid + '/mask-' + side_name
            out[name] = fs
            charts[name] = dict(source_surfaces=[sid] * len(fs), source_faces=parent,
                normal_signs=np.ones(len(fs), int),
                barycentrics=np.array([barycentric(q[face], p[original[int(i)]])
                                       for face, i in zip(fs, parent)]))
            for pid, sign in use_map.items():
                peers = [x for x in use_map if x != pid]
                use = dict(surface=name, side=sign)
                if peers:
                    assert len(peers) == 1
                    use['neighbor'] = peers[0]
                output_cells[pid].append(use)
    source_surfaces = dict(surfaces)
    source_normals, source_uvs = dict(normals), dict(uvs)
    key = 'COMPUTE:partition-mask'
    source_surfaces[key] = cf + len(p)
    tri = v[cf]
    nn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    nn /= np.linalg.norm(nn, axis=1)[:, None]
    source_normals[key] = np.repeat(np.rint(nn * 32767).astype('<i2')[:, None], 3, axis=1)
    source_uvs[key] = np.zeros((len(cf), 3, 2), '<f4')
    for pid in sorted(selected):
        fs, parent = cuts[pid, 'outside']
        if not len(fs):
            continue
        name = pid + '/mask-interface'
        out[name] = fs
        charts[name] = dict(source_surfaces=[key] * len(fs), source_faces=parent,
            normal_signs=-np.ones(len(fs), int),
            barycentrics=np.array([barycentric(q[face], v[cf[int(i)]]) for face, i in zip(fs, parent)]))
        output_cells[pid].append(dict(surface=name, side=1, **(dict(neighbor=region) if region else {})))
        if region:
            output_cells[region].append(dict(surface=name, side=-1, neighbor=pid))
    if plans:
        out, ancestry, _ = refine_surfaces(bank, out, plans)
        q = bank.positions
        for sid, fs in out.items():
            old = charts[sid]
            ids = ancestry[sid]
            families = np.asarray(old['source_surfaces'], object)[ids]
            parent = old['source_faces'][ids]
            charts[sid] = dict(source_surfaces=families.tolist(), source_faces=parent,
                normal_signs=old['normal_signs'][ids], barycentrics=np.array([
                    barycentric(q[face], np.r_[p, v][source_surfaces[family][int(i)]])
                    for face, family, i in zip(fs, families, parent)]))
    result = dict(positions=q, surfaces=out, cells=output_cells, charts=charts)
    words, tex = replay_corner_fields(np.r_[p, v], source_surfaces, source_normals, source_uvs, result)
    conformed, conformity = conform_surface_boundary_edges(q, out, output_cells,
        tolerance=tolerance, record_face_parents=True)
    ancestry = conformity.pop('sourceFaceParents')
    for sid, fs in conformed.items():
        if np.array_equal(fs, out[sid]):
            continue
        parent = ancestry[sid]
        pseudo = dict(positions=q, surfaces={sid: fs}, charts={sid: dict(
            source_surfaces=[sid] * len(fs), source_faces=parent,
            normal_signs=np.ones(len(fs), int), barycentrics=np.array([
                barycentric(q[face], q[out[sid][int(i)]]) for face, i in zip(fs, parent)]))})
        ww, tt = replay_corner_fields(q, {sid: out[sid]}, {sid: words[sid]}, {sid: tex[sid]}, pseudo)
        words[sid], tex[sid] = ww[sid], tt[sid]
    deadline()
    if output_checkpoint:
        output_checkpoint(dict(positions=q, surfaces=conformed, cells=output_cells,
            normals=words, uvs=tex, selected_volumes_m3=volumes))
    for pid, uses in output_cells.items():
        if not uses:
            raise ValueError('A whole tissue cell vanished: ' + pid)
        fs = np.concatenate([conformed[u['surface']] if u['side'] == 1
                             else conformed[u['surface']][:, [0, 2, 1]] for u in uses])
        try:
            solid, _ = make_solid(q, fs)
        except ValueError as exc:
            raise ValueError('Canonical paired-partition topology failed for ' + pid) from exc
        if pid == region:
            expected = sum(r['inside'] for r in volumes.values())
        elif pid in selected:
            expected = volumes[pid]['outside']
        else:
            old_faces = np.concatenate([surfaces[u['surface']] if u['side'] == 1
                else surfaces[u['surface']][:, [0, 2, 1]] for u in cells[pid]])
            expected = float(make_solid(p, old_faces)[0].volume())
        if abs(float(solid.volume()) - expected) > max(1e-18, abs(expected) * 1e-9):
            raise ValueError('Canonical paired-partition volume differs: ' + pid)
    return dict(positions=q, surfaces=conformed, cells=output_cells, normals=words, uvs=tex,
        evidence=dict(paired_boolean=True, selected_volumes_m3=volumes,
            protected_source_regions_preserved=True, source_prefix_literal=np.array_equal(q[:len(p)], p),
            classifier_cap_contributions=0, source_region_reassigned_to=region,
            explicitly_replayed_coincident_source_fragments=source_boundary_replays,
            maximum_reciprocal_boundary_error_m=max(comparisons, default=0.),
            boundary_conformity=conformity, elapsed_seconds=time.monotonic()-started,
            anatomical_acceptance=False))
