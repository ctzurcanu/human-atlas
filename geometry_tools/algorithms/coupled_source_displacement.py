"""Bounded source-edge displacement across explicitly matched source mouths.

Literal aliases, verified complete source-mouth pairs, and explicitly bounded
nonjoining offset pairs may couple variables. Pair members retain offsets;
the solver imposes a common displacement, not an inferred geometric weld.
This is a proposal field.  It certifies neither calibre nor usable anatomy.
"""
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import spsolve

from geometry_tools.source_port_dependencies import check_source_mouth
from .source_offset_pairs import check_offset_pair


def _triangle_fem_blocks(vertices, faces):
    """Original Float64 triangle ``area * grad(phi) @ grad(phi).T`` blocks.

    A local Gram matrix is positive semidefinite, including obtuse triangles
    whose individual cotangent edge weights need not be positive. No target
    geometry, normal words or absolute-position smoothing enter this operator.
    """
    p, f = np.asarray(vertices, dtype=float), np.asarray(faces)
    if (p.ndim != 2 or p.shape[1] != 3 or not np.isfinite(p).all()
            or len(p) > 150000 or f.ndim != 2 or f.shape[1] != 3
            or not np.issubdtype(f.dtype, np.integer) or not 1 <= len(f) <= 600000
            or f.min() < 0 or f.max() >= len(p)):
        raise ValueError('Finite bounded original triangles required for source FEM energy')
    tri = p[f]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    twice_area = np.hypot.reduce(cross, axis=1)
    if np.any(twice_area <= 0) or not np.isfinite(twice_area).all():
        raise ValueError('Source FEM requires nondegenerate original facets')
    normal = cross / twice_area[:, None]
    opposite_edges = np.stack([tri[:, 2] - tri[:, 1], tri[:, 0] - tri[:, 2],
                               tri[:, 1] - tri[:, 0]], axis=1)
    gradients = np.cross(normal[:, None], opposite_edges) / twice_area[:, None, None]
    blocks = (.5 * twice_area[:, None, None]) * (gradients @ gradients.transpose(0, 2, 1))
    if not np.isfinite(blocks).all():
        raise ValueError('Original source FEM gradient energy is not representable')
    return blocks, dict(facets=len(f), minimumSourceArea_m2=float(twice_area.min() / 2),
        maximumSourceArea_m2=float(twice_area.max() / 2),
        positiveOffDiagonalEntries=int(np.count_nonzero(blocks[:, [0, 0, 1], [1, 2, 2]] > 0)),
        positiveSemidefiniteByAreaWeightedGradientGram=True,
        sourceNormalsFromOriginalWorldFloat64Facets=True)


def solve_coupled_source_displacement(meshes, interfaces, *, prescribed, fixed,
                                      declared_target_tolerance_m,
                                      pair_mode='preserve-source-offsets',
                                      energy='inverse-edge-length',
                                      offset_pairs=(), declared_offset_target_tolerance_m=0.,
                                      wall_seconds=3.):
    """Return ``({owner: posed raw coordinates}, evidence)`` without mutation.

    ``prescribed[owner]`` has ``rawIds`` and ``targetPoints_m``; ``fixed[owner]``
    is a list of raw control IDs.  An interface has ``first``, ``second``,
    ``firstRing`` and ``secondRing``.  Every pair is independently checked with
    the existing source-mouth gate before its corresponding variables join.

    The caller must derive the explicit target tolerance from its prescribed
    map and original source-mouth gaps, plus numeric roundoff.  It is used only
    to reconcile targets on *already declared* cross-owner correspondences.
    Exact aliases within an owner and fixed constraints remain hard.  Native
    sparse calls additionally require the workspace's external job supervisor.
    ``energy`` selects either the original edge-graph displacement energy or
    area-weighted ORIGINAL triangle-gradient FEM energy. This choice changes
    neither the prescribed pose nor the source correspondence/fixed controls.

    Optional ``offset_pairs`` are geometric co-move declarations, never joins.
    Each additionally requires ``max_source_gap_m`` and
    ``max_source_relative_area_error``. They must not overlap verified controls.
    Their separate ``declared_offset_target_tolerance_m`` reconciles only their
    target displacements; it never relaxes verified-interface/fixed constraints.
    The resulting FEM field retains each ORIGINAL offset vector, whereas a
    prior rigid section projection may rotate that vector. That distinction
    and any target replay residual remain explicit in the returned evidence.
    """
    started = time.monotonic()
    tolerance = float(declared_target_tolerance_m)
    offset_tolerance=float(declared_offset_target_tolerance_m)
    if (not np.isfinite(wall_seconds) or not 0 < wall_seconds <= 5
            or not np.isfinite(tolerance) or tolerance < 0
            or not np.isfinite(offset_tolerance) or offset_tolerance < 0
            or pair_mode != 'preserve-source-offsets'
            or energy not in ['inverse-edge-length', 'source-triangle-fem']):
        raise ValueError('Explicit finite target tolerance, source-offset mode and wall limit at most five seconds required')
    if (not isinstance(meshes, dict) or not 1 <= len(meshes) <= 128
            or not all(isinstance(pid, str) for pid in meshes)
            or not isinstance(interfaces, (list, tuple)) or len(interfaces) > 2048
            or not isinstance(offset_pairs,(list,tuple)) or len(offset_pairs)>2048
            or not isinstance(prescribed, dict) or not isinstance(fixed, dict)
            or set(prescribed) - set(meshes) or set(fixed) - set(meshes)):
        raise ValueError('Bounded explicit source owners, interfaces and owned constraints required')

    def deadline():
        if time.monotonic() - started >= wall_seconds:
            raise TimeoutError('Coupled source-displacement deadline; no complete field emitted')

    def raw_ids(values, pid, *, minimum=0):
        raw = np.asarray(values)
        if len(raw.shape) == 1 and len(raw) == 0 and minimum == 0:
            return np.empty(0, dtype=int)
        if (raw.ndim != 1 or not np.issubdtype(raw.dtype, np.integer)
                or len(raw) < minimum or raw.min(initial=0) < 0
                or raw.max(initial=0) >= len(meshes[pid].vertices)):
            raise ValueError('Explicit valid integer source control IDs required')
        return raw.astype(int, copy=False)

    deadline()
    owner_ids = sorted(meshes)
    frames = {meshes[pid].frame for pid in owner_ids}
    if len(frames) != 1:
        raise ValueError('All owners require the same explicit source coordinate frame')
    raw_count = sum(len(meshes[pid].vertices) for pid in owner_ids)
    face_count = sum(len(meshes[pid].faces) for pid in owner_ids)
    if not 1 <= raw_count <= 150000 or not 1 <= face_count <= 600000:
        raise ValueError('At most 150,000 source controls and 600,000 original faces required')

    physical_points, physical_owners, maps, source_edges, source_weights = [], [], {}, [], []
    cursor = 0
    for number, pid in enumerate(owner_ids):
        deadline()
        mesh = meshes[pid]
        if len(mesh.vertices) < 3 or not len(mesh.faces):
            raise ValueError('Every declared owner needs actual source triangles')
        p, inverse = np.unique(mesh.vertices, axis=0, return_inverse=True)
        maps[pid] = inverse + cursor
        physical_points.append(p)
        physical_owners.append(np.full(len(p), number, dtype=int))
        faces = inverse[mesh.faces]
        edges = np.unique(np.sort(faces[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2), axis=1), axis=0)
        edges = edges[edges[:, 0] != edges[:, 1]]
        lengths = np.hypot.reduce(p[edges[:, 0]] - p[edges[:, 1]], axis=1)
        if np.any(lengths <= 0) or not np.isfinite(lengths).all():
            raise ValueError('Original source edges need finite positive lengths')
        source_edges.append(edges + cursor)
        with np.errstate(over='ignore'):
            weights = 1 / lengths
        if not np.isfinite(weights).all():
            raise ValueError('Original source inverse-edge weights must be finite')
        source_weights.append(weights)
        cursor += len(p)
    points = np.vstack(physical_points)
    owners = np.concatenate(physical_owners)
    parent = np.arange(len(points))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return int(i)

    checked = [];verified_controls=set();offset_controls=set();offset_checked=[]
    for row in interfaces:
        deadline()
        first, second = row['first'], row['second']
        if first not in meshes or second not in meshes or first == second:
            raise ValueError('Each declared interface requires two distinct existing owners')
        a = raw_ids(row['firstRing'], first, minimum=3)
        b = raw_ids(row['secondRing'], second, minimum=3)
        if len(a) != len(b):
            raise ValueError('Only explicit same-arity source-mouth correspondences supported')
        proof = check_source_mouth(meshes[first], meshes[second], a, b)
        deadline()
        if not proof['passed']:
            raise ValueError('Declared complete source-mouth correspondence failed its source gate')
        checked.append(dict(first=first, second=second, firstRing=a.tolist(),
                            secondRing=b.tolist(), sourceProof=proof))
        verified_controls.update(map(int,np.r_[maps[first][a],maps[second][b]]))
        for i, j in zip(maps[first][a], maps[second][b]):
            parent[root(int(j))] = root(int(i))
    for row in offset_pairs:
        deadline();first,second=row['first'],row['second']
        if first not in meshes or second not in meshes or first==second:
            raise ValueError('Offset pair needs two distinct existing source owners')
        a=raw_ids(row['firstRing'],first,minimum=3);b=raw_ids(row['secondRing'],second,minimum=3)
        if len(a)!=len(b):raise ValueError('Offset pairs require same-arity complete source boundaries')
        ids=set(map(int,np.r_[maps[first][a],maps[second][b]]))
        if ids&verified_controls:raise ValueError('Offset-pair controls overlap verified-interface controls')
        proof=check_offset_pair(meshes[first],meshes[second],a,b,row);deadline()
        offset_checked.append(dict(first=first,second=second,firstRing=a.tolist(),secondRing=b.tolist(),**proof))
        offset_controls.update(ids)
        for i,j in zip(maps[first][a],maps[second][b]):parent[root(int(j))]=root(int(i))
    deadline()
    roots = np.array([root(i) for i in range(len(points))], dtype=int)
    _, groups = np.unique(roots, return_inverse=True)
    group_count = int(groups.max()) + 1
    offset_groups=set(map(int,groups[np.array(sorted(offset_controls),int)]))if offset_controls else set()
    # Inconsistent correspondence cycles cannot silently merge distinct source
    # controls of one owner, including any complete source-ring corners.
    _, counts = np.unique(np.column_stack([groups, owners]), axis=0, return_counts=True)
    if np.any(counts > 1):
        raise ValueError('Cross-owner correspondence cycles merge distinct controls within one owner')

    fixed_rows = {}
    fixed_groups = set()
    for pid, values in fixed.items():
        raw = raw_ids(values, pid)
        fixed_rows[pid] = raw
        for i in raw:
            fixed_groups.add(int(groups[maps[pid][i]]))

    constraints = {}
    exact_alias_targets = {}
    for pid, declaration in prescribed.items():
        deadline()
        raw = raw_ids(declaration['rawIds'], pid)
        targets = np.asarray(declaration['targetPoints_m'], dtype=float)
        if targets.shape != (len(raw), 3) or not np.isfinite(targets).all():
            raise ValueError('Every prescribed source control needs one finite target point')
        for i, target in zip(raw, targets):
            physical = int(maps[pid][i])
            previous = exact_alias_targets.setdefault(physical, target)
            if not np.array_equal(previous, target):
                raise ValueError('Exact aliases within an owner have incompatible prescribed targets')
            group = int(groups[physical])
            delta = target - meshes[pid].vertices[i]
            if group in fixed_groups and not np.array_equal(delta, np.zeros(3)):
                raise ValueError('A prescribed target conflicts with a hard fixed source control')
            constraints.setdefault(group, []).append((pid, int(i), target, delta))

    displacement = np.zeros((group_count, 3))
    canonical = {}
    maximum_disagreement = 0.
    for group, rows in constraints.items():
        chosen = rows[0]
        errors = np.array([np.linalg.norm(row[3] - chosen[3]) for row in rows])
        worst = float(errors.max(initial=0))
        maximum_disagreement = max(maximum_disagreement, worst)
        group_tolerance=offset_tolerance if group in offset_groups else tolerance
        if worst > group_tolerance:
            raise ValueError('Coupled declared target displacements exceed the explicit agreement tolerance')
        if group not in fixed_groups:
            displacement[group] = chosen[3]
            canonical[group] = chosen

    deadline()
    edges = np.vstack(source_edges)
    weights = np.concatenate(source_weights)
    mapped = groups[edges]
    if np.any(mapped[:, 0] == mapped[:, 1]):
        raise ValueError('Explicit cross-owner correspondence collapses an original source edge')
    a, b = mapped.T
    adjacency = coo_matrix((np.r_[weights, weights], (np.r_[a, b], np.r_[b, a])),
                           shape=(group_count, group_count)).tocsr()
    _, labels = connected_components(adjacency, directed=False)
    boundary = set(fixed_groups) | set(constraints)
    protected_labels = labels[np.array(sorted(boundary), dtype=int)]
    literal_groups = np.flatnonzero(~np.isin(labels, protected_labels))
    boundary.update(map(int, literal_groups))
    orphan_labels, orphan_counts = np.unique(labels[literal_groups], return_counts=True)
    unanchored = [dict(component=int(label), physicalGroups=int(count))
                  for label, count in zip(orphan_labels, orphan_counts)]
    boundary = np.array(sorted(boundary), dtype=int)
    free = np.setdiff1d(np.arange(group_count), boundary)
    fem_evidence = []
    if energy == 'inverse-edge-length':
        laplacian = diags(np.asarray(adjacency.sum(1)).ravel()) - adjacency
    else:
        row_indices, column_indices, entries = [], [], []
        for pid in owner_ids:
            deadline()
            mesh = meshes[pid]
            blocks, detail = _triangle_fem_blocks(mesh.vertices, mesh.faces)
            face_groups = groups[maps[pid][mesh.faces]]
            row_indices.append(np.broadcast_to(face_groups[:, :, None], blocks.shape).ravel())
            column_indices.append(np.broadcast_to(face_groups[:, None, :], blocks.shape).ravel())
            entries.append(blocks.ravel())
            fem_evidence.append(dict(owner=pid, **detail))
            deadline()
        laplacian = coo_matrix((np.concatenate(entries),
                                (np.concatenate(row_indices), np.concatenate(column_indices))),
                               shape=(group_count, group_count)).tocsr()
        if not np.isfinite(laplacian.data).all():
            raise ValueError('Assembled original source FEM operator is not finite')
        deadline()
    residual = 0.
    if len(free):
        deadline()
        matrix = laplacian[free][:, free].tocsc()
        rhs = -laplacian[free][:, boundary] @ displacement[boundary]
        displacement[free] = spsolve(matrix, rhs)
        deadline()
        if not np.isfinite(displacement).all():
            raise ValueError('Coupled original-source displacement solve is not finite')
        residual = float(np.max(np.abs(matrix @ displacement[free] - rhs), initial=0))

    canonical_base = np.zeros((group_count, 3))
    canonical_target = np.zeros((group_count, 3))
    target_word_group = np.zeros(group_count, dtype=bool)
    for group, (pid, raw, target, delta) in canonical.items():
        if not np.array_equal(delta, np.zeros(3)):
            canonical_base[group] = meshes[pid].vertices[raw]
            canonical_target[group] = target
            target_word_group[group] = True
    posed = {}
    for pid in owner_ids:
        deadline()
        mesh = meshes[pid]
        raw_groups = groups[maps[pid]]
        posed[pid] = mesh.vertices + displacement[raw_groups]
        # Set canonical hard target words literally, retaining each paired
        # source offset.  This changes addition order, not the common field.
        mask = target_word_group[raw_groups]
        selected = raw_groups[mask]
        posed[pid][mask] = canonical_target[selected] + (mesh.vertices[mask] - canonical_base[selected])

    snapped, maximum_replay = [], 0.
    for group, rows in constraints.items():
        records = []
        for pid, raw, target, delta in rows:
            error = float(np.linalg.norm(posed[pid][raw] - target))
            maximum_replay = max(maximum_replay, error)
            records.append(dict(owner=pid, rawId=raw, targetPoint_m=target.tolist(),
                                targetDisplacement_m=delta.tolist(), replayResidual_m=error))
        group_tolerance=offset_tolerance if group in offset_groups else tolerance
        if max(row['replayResidual_m'] for row in records) > group_tolerance:
            raise ValueError('Literal target replay exceeds the explicit declared tolerance')
        if any(row['replayResidual_m'] > 0 for row in records):
            snapped.append(dict(group=group, fixedSourcePriority=group in fixed_groups,
                                canonicalOwner=rows[0][0], canonicalRawId=rows[0][1],
                                constraints=records))
    for pid, raw in fixed_rows.items():
        if not np.array_equal(posed[pid][raw], meshes[pid].vertices[raw]):
            raise ValueError('Hard fixed source coordinates moved')
    maximum_pair_delta_error = 0.
    for row in checked:
        first, second = row['first'], row['second']
        a, b = row['firstRing'], row['secondRing']
        delta_a = posed[first][a] - meshes[first].vertices[a]
        delta_b = posed[second][b] - meshes[second].vertices[b]
        maximum_pair_delta_error = max(maximum_pair_delta_error,
                                       float(np.linalg.norm(delta_a - delta_b, axis=1).max()))
    for row in offset_checked:
        first,second=row['first'],row['second'];a,b=row['firstRing'],row['secondRing']
        original=meshes[second].vertices[b]-meshes[first].vertices[a]
        actual=posed[second][b]-posed[first][a]
        row.update(originalOffsetVectorPreserved=True,offsetRotationApplied=False,
                   maximumOffsetVectorRoundoff_m=float(np.linalg.norm(actual-original,axis=1).max()))
    deadline()
    return posed, dict(complete=True, method='original inverse-edge-length coupled Dirichlet displacement' if energy == 'inverse-edge-length' else 'original source-triangle FEM coupled Dirichlet displacement',
        energy=energy, sourceTriangleFEM=fem_evidence,
        sourceEdgeAdjacencyRetainedForComponents=True,
        pairMode=pair_mode, owners=owner_ids, sourceRawControls=raw_count,
        sourceFaces=face_count, sourcePhysicalControls=len(points), coupledPhysicalGroups=group_count,
        freePhysicalGroups=len(free), originalSourceEdges=len(edges), interfaces=checked,
        offsetPairProofs=offset_checked,declaredOffsetTargetTolerance_m=offset_tolerance,
        fixedRawControlIds={pid: raw.tolist() for pid, raw in fixed_rows.items()},
        prescribedRawControlCount=sum(len(row['rawIds']) for row in prescribed.values()),
        declaredTargetTolerance_m=tolerance, maximumDeclaredDisplacementDisagreement_m=maximum_disagreement,
        maximumPrescribedReplayResidual_m=maximum_replay, snappedDeclaredGroups=snapped,
        prescribedCanonicalPolicy='First explicitly declared target supplies the common displacement; literal fixed source controls override and cannot conflict',
        maximumPairedDisplacementRoundoff_m=maximum_pair_delta_error,
        equationMaximumResidual=residual, unanchoredSourceComponentsRetainedLiteral=unanchored,
        withinOwnerAliasesExactOnly=True, crossOwnerCouplingFromDeclaredVerifiedMouthsOnly=not bool(offset_checked),
        crossOwnerCouplingFromVerifiedInterfacesAndExplicitOffsetPairsOnly=True,
        sourceMouthOffsetsRetained=True, sourceFacesAndAttributesUntouched=True,
        fixedSourceCoordinatesLiteral=True, sourceGraphWeightsOriginal=True,
        actualSourceMouthReplayRequired=True, triangleQualityAndCalibreVerified=False,
        foreignContactsVerified=False, anatomicalAcceptance=False,
        wallSeconds=wall_seconds, elapsedSeconds=time.monotonic() - started)
