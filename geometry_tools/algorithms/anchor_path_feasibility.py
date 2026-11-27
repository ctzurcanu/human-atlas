"""Necessary distance bounds for explicit anatomical attachment positions.

If every original mesh edge can lengthen by at most ``maximum_edge_extension``,
the distance between two new anchor positions cannot exceed that factor times
any connecting original edge path. The shortest source path gives the strongest
such bound. Passing this test never proves a usable surface, a feasible set of
section poses, collision freedom, preserved calibre, or anatomical correctness.

Only equal source coordinates are aliases. No nearby controls, anatomical
attachments, or mesh edges are inferred. All input arrays remain unchanged.
"""
from __future__ import annotations

import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, dijkstra


def anchor_path_feasibility(source_positions, source_edges, anchor_raw_ids,
                            target_points, *, maximum_edge_extension,
                            wall_seconds=5., maximum_anchors=128):
    """Return only a necessary source-edge stretch test, with literal witnesses.

    Call under the workspace's separately supervised job deadline. Native
    shortest-path calls process one source anchor at a time; the local wall
    deadline is checked before and after each call and each preparation stage.
    Graph and anchor limits also bound the work of any one native call. A
    timeout raises instead of returning an incomplete or passing result.
    Disconnected anchor pairs impose no bound and are reported as such.
    """
    started = time.monotonic()
    if (not np.isfinite(wall_seconds) or not 0 < wall_seconds <= 5 or
            not isinstance(maximum_anchors, int) or isinstance(maximum_anchors, bool)
            or not 2 <= maximum_anchors <= 128):
        raise ValueError('Explicit wall deadline at most five seconds and 2–128 anchors required')
    factor = float(maximum_edge_extension)
    if not np.isfinite(factor) or factor < 1:
        raise ValueError('Finite maximum source-edge extension factor at least one required')

    def deadline():
        if time.monotonic() - started >= wall_seconds:
            raise TimeoutError('Anchor source-path deadline reached; no feasibility result emitted')

    deadline()
    p = np.asarray(source_positions, dtype=float)
    edges = np.asarray(source_edges)
    raw = np.asarray(anchor_raw_ids)
    target = np.asarray(target_points, dtype=float)
    if (p.ndim != 2 or p.shape[1] != 3 or not len(p) or
            len(p) > 250000 or not np.isfinite(p).all()):
        raise ValueError('Finite source coordinates, at most 250,000 controls, required')
    if (edges.ndim != 2 or edges.shape[1] != 2 or
            not np.issubdtype(edges.dtype, np.integer) or len(edges) > 1500000 or
            (len(edges) and (edges.min() < 0 or edges.max() >= len(p)))):
        raise ValueError('At most 1,500,000 explicit integer source edges within the source bank required')
    if (raw.ndim != 1 or not np.issubdtype(raw.dtype, np.integer) or
            not 2 <= len(raw) <= maximum_anchors or raw.min() < 0 or
            raw.max() >= len(p) or target.shape != (len(raw), 3) or
            not np.isfinite(target).all()):
        raise ValueError('2–128 explicit valid anchor IDs and corresponding finite target points required')
    deadline()

    # Alias reduction is exact, including duplicate raw IDs in the declaration.
    declared_count = len(raw)
    physical, inverse = np.unique(p, axis=0, return_inverse=True)
    anchor_physical = inverse[raw]
    declared = {}
    for i, pid in enumerate(anchor_physical):
        old = declared.setdefault(int(pid), i)
        if not np.array_equal(target[old], target[i]):
            raise ValueError('Exact source-coordinate aliases have inconsistent target positions')
    chosen = np.array(list(declared.values()), dtype=int)
    anchor_physical = anchor_physical[chosen]
    raw = raw[chosen]
    target = target[chosen]
    mapped = np.sort(inverse[edges], axis=1)
    collapsed = int(np.count_nonzero(mapped[:, 0] == mapped[:, 1]))
    mapped = np.unique(mapped[mapped[:, 0] != mapped[:, 1]], axis=0)
    length = np.hypot.reduce(physical[mapped[:, 0]] - physical[mapped[:, 1]], axis=1)
    if not np.isfinite(length).all() or np.any(length <= 0):
        raise ValueError('Every noncollapsed source edge needs a representable positive finite length')
    graph = coo_matrix((np.r_[length, length],
                        (np.r_[mapped[:, 0], mapped[:, 1]],
                         np.r_[mapped[:, 1], mapped[:, 0]])),
                       shape=(len(physical), len(physical))).tocsr()
    _, labels = connected_components(graph, directed=False)
    deadline()

    checked = disconnected = violations = 0
    largest = None
    witnesses = []
    coordinate_scale = max(float(np.abs(p).max()), float(np.abs(target).max()))
    for i, source in enumerate(anchor_physical):
        if i + 1 >= len(anchor_physical):
            break
        same = np.flatnonzero(labels[anchor_physical[i + 1:]] == labels[source]) + i + 1
        disconnected += len(anchor_physical) - i - 1 - len(same)
        if not len(same):
            continue
        deadline()
        distances = dijkstra(graph, directed=False, indices=int(source))
        deadline()
        for j in same:
            source_path = float(distances[anchor_physical[j]])
            # A known graph path must remain finite in the numeric length type.
            if not np.isfinite(source_path) or source_path <= 0:
                raise ValueError('A connected distinct anchor pair has no representable positive source path')
            separation = float(np.hypot.reduce(target[i] - target[j]))
            limit = factor * source_path
            if not np.isfinite(separation) or not np.isfinite(limit):
                raise ValueError('Anchor distance or permitted source path exceeds finite numeric range')
            # Roundoff protection is coordinate-scale based, never an anatomical
            # tolerance or permission to merge close source controls.
            roundoff = 64 * np.finfo(float).eps * (coordinate_scale + separation + limit)
            ratio = separation / source_path
            witness = dict(firstRawId=int(raw[i]), secondRawId=int(raw[j]),
                           originalShortestSourceEdgePath_m=source_path,
                           desiredAnchorSeparation_m=separation,
                           requiredEdgeExtensionLowerBound=float(ratio),
                           maximumAllowedEdgeExtension=factor,
                           numericalRoundoff_m=float(roundoff),
                           violatesNecessaryBound=bool(separation > limit + roundoff))
            checked += 1
            if largest is None or ratio > largest['requiredEdgeExtensionLowerBound']:
                largest = witness
            if witness['violatesNecessaryBound']:
                violations += 1
                if len(witnesses) < 32:
                    witnesses.append(witness)
        deadline()

    return dict(complete=True, passedNecessaryBound=violations == 0,
                rejectedAsImpossibleUnderEdgeBound=violations > 0,
                maximumAllowedEdgeExtension=factor,
                declaredRawAnchors=declared_count,
                physicalAnchors=len(chosen), sourceRawControls=len(p),
                sourcePhysicalControls=len(physical), sourceEdges=len(mapped),
                collapsedExactAliasEdges=collapsed,
                checkedConnectedPairs=checked, disconnectedPairsIgnored=disconnected,
                violatingPairs=violations, violatingWitnesses=witnesses,
                violatingWitnessesTruncated=violations > len(witnesses),
                largestRequiredExtensionWitness=largest,
                exactCoordinateAliasesOnly=True, sourceGeometryChanged=False,
                sufficientFeasibilityEstablished=False, anatomicalAcceptance=False,
                limitation='Necessary source-edge distance bound only; disconnected pairs impose no constraint. Passing does not establish a feasible surface, calibre, source sections, contacts or anatomy.',
                elapsedSeconds=time.monotonic() - started)


def mesh_anchor_path_feasibility(mesh, anchor_raw_ids, target_points, **kwargs):
    """Use every actual source triangle edge, retaining raw anchor identities."""
    f = np.asarray(mesh.faces)
    edges = f[:, [0, 1, 1, 2, 2, 0]].reshape(-1, 2)
    return anchor_path_feasibility(mesh.vertices, edges, anchor_raw_ids,
                                   target_points, **kwargs)
