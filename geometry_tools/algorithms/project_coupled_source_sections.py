"""Project explicitly listed source sections onto bounded proper similarities.

Actual overlapping source cycles and verified paired mouths share one pose.
Nothing recruits an unlisted branch core or an unrelated coincident owner.
The output supplies constraints for a subsequent source-displacement proposal;
that later field still needs a fresh section, triangle and context audit.
"""
import time

import numpy as np

from geometry_tools.source_port_dependencies import check_source_mouth
from .source_offset_pairs import check_offset_pair


def propose_section_targets(meshes, interfaces, profiles, proposed_positions,
                            fixed, *, maximum_uniform_scale,
                            minimum_uniform_scale=1., offset_pairs=(), wall_seconds=3.):
    """Return ``(prescribed, expanded_fixed, proof)`` without changing sources.

    ``profiles[owner]`` contains explicit raw-ID cycles. ``fixed[owner]`` lists
    protected raw controls, including caller-owned attachment stars. Targets
    use the coupled-displacement API: ``rawIds`` and ``targetPoints_m``.
    Optional ``offset_pairs`` name complete opposed same-arity boundaries and
    explicit ``max_source_gap_m`` / ``max_source_relative_area_error`` bounds.
    These nonjoining pairs form one proper UNIT-scale cell: their original
    separated points are not welded and their gap vectors rotate together.
    They cannot share a physical control or projected cell with verified joins.
    """
    started = time.monotonic()
    minimum_scale, maximum_scale = float(minimum_uniform_scale), float(maximum_uniform_scale)
    if (not np.isfinite(wall_seconds) or not 0 < wall_seconds <= 5
            or not np.isfinite(minimum_scale) or not np.isfinite(maximum_scale)
            or not .5 <= minimum_scale <= maximum_scale <= 2):
        raise ValueError('Explicit wall limit at most five seconds and finite scale range 0.5 <= minimum <= maximum <= 2 required')
    if (not isinstance(meshes, dict) or not 1 <= len(meshes) <= 128
            or not all(isinstance(pid, str) for pid in meshes)
            or not isinstance(profiles, dict) or not isinstance(fixed, dict)
            or not isinstance(proposed_positions, dict)
            or set(proposed_positions) != set(meshes)
            or set(profiles) - set(meshes) or set(fixed) - set(meshes)
            or not isinstance(interfaces, (list, tuple)) or len(interfaces) > 2048
            or not isinstance(offset_pairs,(list,tuple))or len(offset_pairs)>2048):
        raise ValueError('Explicit bounded owners, source cycles, proposed coordinates and fixed controls required')

    def deadline():
        if time.monotonic() - started >= wall_seconds:
            raise TimeoutError('Coupled section-projection deadline; no complete targets emitted')

    def raw_ids(pid, values, minimum=0):
        ids = np.asarray(values)
        if ids.ndim == 1 and not len(ids) and minimum == 0:
            return np.empty(0, dtype=int)
        if (ids.ndim != 1 or not np.issubdtype(ids.dtype, np.integer)
                or len(ids) < minimum or ids.min(initial=0) < 0
                or ids.max(initial=0) >= len(meshes[pid].vertices)):
            raise ValueError('Explicit valid integer source control IDs required')
        return ids.astype(int, copy=False)

    deadline()
    owner_ids = sorted(meshes)
    raw_count = sum(len(mesh.vertices) for mesh in meshes.values())
    face_count = sum(len(mesh.faces) for mesh in meshes.values())
    if not 1 <= raw_count <= 150000 or not 1 <= face_count <= 600000:
        raise ValueError('At most 150,000 source controls and 600,000 source faces required')
    if len({mesh.frame for mesh in meshes.values()}) != 1:
        raise ValueError('One explicit source coordinate frame required')
    points, proposed, maps, edges_by_owner = [], [], {}, {}
    base = 150001
    cursor = 0
    for pid in owner_ids:
        deadline()
        mesh = meshes[pid]
        q = np.asarray(proposed_positions[pid], dtype=float)
        if not len(mesh.faces) or q.shape != mesh.vertices.shape or not np.isfinite(q).all():
            raise ValueError('Actual source triangles and corresponding finite proposed points required')
        physical, representative, inverse = np.unique(mesh.vertices, axis=0, return_index=True, return_inverse=True)
        # Exact source aliases are one physical observation, rather than extra
        # fitting weight from UV/property chart duplication.
        if not np.array_equal(q, q[representative][inverse]):
            raise ValueError('Exact source aliases have inconsistent proposed positions')
        maps[pid] = inverse + cursor
        points.append(physical)
        proposed.append(q[representative])
        faces = maps[pid][mesh.faces]
        edges = np.unique(np.sort(faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1), axis=0)
        edges = edges[edges[:, 0] != edges[:, 1]]
        edges_by_owner[pid] = np.unique(edges[:, 0] * base + edges[:, 1])
        cursor += len(physical)
    p, q = np.vstack(points), np.vstack(proposed)
    parent = np.arange(len(p))
    active = np.zeros(len(p), dtype=bool)

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return int(i)

    def merge(ids):
        ids = np.unique(ids)
        anchor = root(int(ids[0]))
        for i in ids[1:]:
            parent[root(int(i))] = anchor
        active[ids] = True

    profile_rows = []
    for pid, rings in profiles.items():
        if not isinstance(rings, (list, tuple)) or len(rings) > 30000:
            raise ValueError('A bounded explicit source-cycle list required')
        for index, ring in enumerate(rings):
            deadline()
            raw = raw_ids(pid, ring, minimum=3)
            ids = maps[pid][raw]
            if len(np.unique(ids)) != len(ids):
                raise ValueError('A source section cycle repeats a physical control')
            edges = np.sort(np.column_stack([ids, np.roll(ids, -1)]), axis=1)
            codes = edges[:, 0] * base + edges[:, 1]
            bank = edges_by_owner[pid]
            locations = np.searchsorted(bank, codes)
            if np.any(locations >= len(bank)) or not np.array_equal(bank[locations], codes):
                raise ValueError('Every supplied section cycle needs its actual original source edges')
            merge(ids)
            profile_rows.append(dict(owner=pid, index=index, rawIds=raw.tolist()))
    if len(profile_rows) > 30000:
        raise ValueError('At most 30,000 explicit source cycles required')

    interface_rows = [];verified_controls=set();offset_rows=[]
    for row in interfaces:
        deadline()
        first, second = row['first'], row['second']
        if first not in meshes or second not in meshes or first == second:
            raise ValueError('Declared source mouths need two distinct existing owners')
        a, b = raw_ids(first, row['firstRing'], 3), raw_ids(second, row['secondRing'], 3)
        if len(a) != len(b):
            raise ValueError('Only explicit same-arity source-mouth pairs supported')
        proof = check_source_mouth(meshes[first], meshes[second], a, b)
        deadline()
        if not proof['passed']:
            raise ValueError('Declared complete source mouth fails its source gate')
        merge(np.r_[maps[first][a], maps[second][b]])
        interface_rows.append(dict(first=first, second=second, firstRing=a.tolist(),
                                   secondRing=b.tolist(), sourceProof=proof))
        verified_controls.update(map(int,np.r_[maps[first][a],maps[second][b]]))
    for row in offset_pairs:
        deadline();first,second=row['first'],row['second']
        if first not in meshes or second not in meshes or first==second:
            raise ValueError('Offset pair needs two distinct existing source owners')
        a,b=raw_ids(first,row['firstRing'],3),raw_ids(second,row['secondRing'],3)
        if len(a)!=len(b):raise ValueError('Offset pairs require same-arity complete source boundaries')
        if set(map(int,np.r_[maps[first][a],maps[second][b]]))&verified_controls:
            raise ValueError('Offset-pair controls overlap verified-interface controls')
        proof=check_offset_pair(meshes[first],meshes[second],a,b,row);deadline()
        merge(np.r_[maps[first][a],maps[second][b]])
        offset_rows.append(dict(first=first,second=second,firstRing=a.tolist(),secondRing=b.tolist(),**proof))
    if offset_rows and not minimum_scale<=1.<=maximum_scale:
        raise ValueError('Nonjoining offset-pair cells require unit scale within the declared scale interval')

    fixed_raw = {pid: raw_ids(pid, values) for pid, values in fixed.items()}
    fixed_physical = np.zeros(len(p), dtype=bool)
    for pid, raw in fixed_raw.items():
        fixed_physical[maps[pid][raw]] = True
    actual = np.flatnonzero(active)
    actual_roots = np.array([root(int(i)) for i in actual])
    _, labels = np.unique(actual_roots, return_inverse=True)
    physical_cell = np.full(len(p), -1, dtype=int)
    physical_cell[actual] = labels
    offset_cells={int(physical_cell[maps[row['first']][row['firstRing'][0]]])for row in offset_rows}
    verified_cells={int(physical_cell[maps[row['first']][row['firstRing'][0]]])for row in interface_rows}
    if offset_cells&verified_cells:raise ValueError('Offset-pair cell overlaps a verified-interface cell')
    target = p.copy()
    fixed_cell = np.zeros(int(labels.max(initial=-1)) + 1, dtype=bool)
    cells = []
    order = np.argsort(labels, kind='stable')
    sorted_ids = actual[order]
    cuts = np.r_[0, np.flatnonzero(np.diff(labels[order])) + 1, len(actual)] if len(actual) else np.array([0])
    for cell, (start, stop) in enumerate(zip(cuts[:-1], cuts[1:])):
        deadline()
        ids = sorted_ids[start:stop]
        x, y = p[ids], q[ids]
        center, destination = x.mean(0), y.mean(0)
        is_fixed = bool(fixed_physical[ids].any())
        if is_fixed:
            rotation, scale, unconstrained_scale = np.eye(3), 1., 1.
            destination = center
            posed = x.copy()
            fixed_cell[cell] = True
        else:
            xc, yc = x - center, y - destination
            source_rank = np.linalg.svd(xc, compute_uv=False)
            if len(source_rank) < 2 or source_rank[1] <= 0:
                raise ValueError('A source section cell has no noncollinear shape')
            u, singular, vt = np.linalg.svd(xc.T @ yc)
            correction = np.eye(3)
            correction[-1, -1] = np.linalg.det(vt.T @ u.T)
            rotation = vt.T @ correction @ u.T
            unconstrained_scale = float((singular * np.diag(correction)).sum() / np.square(xc).sum())
            scale = 1. if cell in offset_cells else float(np.clip(unconstrained_scale, minimum_scale, maximum_scale))
            posed = scale * xc @ rotation.T + destination
            if np.linalg.det(rotation) <= 0 or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12, rtol=0):
                raise ValueError('Section projection failed to produce a proper rotation')
        target[ids] = posed
        cells.append(dict(cell=cell, physicalControls=len(ids), fixedBySourceStar=is_fixed,
                          scale=scale, unconstrainedScale=unconstrained_scale,
                          rotation=rotation.tolist(), sourceCenter_m=center.tolist(),
                          targetCenter_m=destination.tolist(), determinant=float(np.linalg.det(rotation)),
                          maximumProposedPointResidual_m=float(np.linalg.norm(posed - y, axis=1).max()),
                          maximumSourceMotion_m=float(np.linalg.norm(posed - x, axis=1).max())))

    prescribed, expanded_fixed = {}, {}
    maximum_edge_error = 0.
    for pid in owner_ids:
        deadline()
        raw_cell = physical_cell[maps[pid]]
        raw = np.flatnonzero(raw_cell >= 0)
        if len(raw):
            prescribed[pid] = dict(rawIds=raw.tolist(), targetPoints_m=target[maps[pid][raw]].tolist())
        fixed_mask = np.zeros(len(meshes[pid].vertices), dtype=bool)
        fixed_mask[fixed_raw.get(pid, [])] = True
        in_cells = raw_cell >= 0
        fixed_mask[in_cells] |= fixed_cell[raw_cell[in_cells]]
        expanded_fixed[pid] = np.flatnonzero(fixed_mask).tolist()
        if not np.array_equal(target[maps[pid][fixed_mask]], meshes[pid].vertices[fixed_mask]):
            raise ValueError('Expanded source-star section cell moved')
        codes = edges_by_owner[pid]
        a, b = codes // base, codes % base
        edge_cells = physical_cell[a]
        same = (edge_cells >= 0) & (edge_cells == physical_cell[b])
        if same.any():
            source_length = np.linalg.norm(p[a[same]] - p[b[same]], axis=1)
            target_length = np.linalg.norm(target[a[same]] - target[b[same]], axis=1)
            scales = np.array([row['scale'] for row in cells])[edge_cells[same]]
            maximum_edge_error = max(maximum_edge_error, float(abs(target_length - scales * source_length).max()))
    for row in profile_rows:
        row['cell'] = int(physical_cell[maps[row['owner']][row['rawIds'][0]]])
    maximum_pair_delta_error = 0.
    for row in interface_rows:
        a, b = maps[row['first']][row['firstRing']], maps[row['second']][row['secondRing']]
        row['cell'] = int(physical_cell[a[0]])
        delta_error = float(np.linalg.norm((target[a] - p[a]) - (target[b] - p[b]), axis=1).max())
        maximum_pair_delta_error = max(maximum_pair_delta_error, delta_error)
        row['maximumPrescribedDisplacementDisagreement_m'] = delta_error
    for row in offset_rows:
        a,b=maps[row['first']][row['firstRing']],maps[row['second']][row['secondRing']]
        cell=int(physical_cell[a[0]]);rotation=np.asarray(cells[cell]['rotation']);original=p[b]-p[a];changed=target[b]-target[a]
        row.update(cell=cell,commonProperRotation=rotation.tolist(),sectionScale=1.,
                   originalOffsetLengthPreserved=True,originalOffsetVectorPreserved=False,
                   maximumRotatedOffsetResidual_m=float(np.linalg.norm(changed-original@rotation.T,axis=1).max()),
                   maximumOffsetLengthRoundoff_m=float(abs(np.linalg.norm(changed,axis=1)-np.linalg.norm(original,axis=1)).max()))
    deadline()
    return prescribed, expanded_fixed, dict(
        complete=True, method='proper bounded similarity on explicit overlapping source-section cells',
        minimumUniformScale=minimum_scale, maximumUniformScale=maximum_scale,
        fixedCellScaleSourcePriority=True, sourceRawControls=raw_count, sourceFaces=face_count,
        listedSourceProfiles=profile_rows, interfaces=interface_rows, offsetPairProofs=offset_rows,cells=cells,
        rigidCells=len(cells), fixedCells=int(fixed_cell.sum()),
        projectedPhysicalControls=len(actual), unlistedCorePhysicalControls=len(p) - len(actual),
        maximumWithinCellSourceEdgeLengthError_m=maximum_edge_error,
        maximumPrescribedMouthDisplacementDisagreement_m=maximum_pair_delta_error,
        withinOwnerAliasesExactOnly=True, crossOwnerCellsFromVerifiedDeclaredMouthsOnly=not bool(offset_rows),
        crossOwnerCellsFromVerifiedInterfacesAndExplicitOffsetPairsOnly=True,
        originalSourceOffsetsRetained=True, fixedSourceStarPriority=True,
        unlistedBranchCoresIncluded=False, sourceFacesAndAttributesUntouched=True,
        fittingWeights='one observation per exact within-owner source control; original counterpart chart offsets retained',
        pairwiseSectionShapePreservedByProperSimilarity=True,
        targetsOnly=True, subsequentFieldReplayRequired=True,
        subsequentSectionAndTriangleAuditRequired=True, foreignContactsVerified=False,
        anatomicalAcceptance=False, wallSeconds=wall_seconds,
        elapsedSeconds=time.monotonic() - started)
