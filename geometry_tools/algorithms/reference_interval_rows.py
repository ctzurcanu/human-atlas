"""Keep gaps from a reviewed pose while measuring all motion from the source.

The SAT plane is chosen on the current triangles. The sparse inequalities are
expressed against the immutable reference triangles, so a subsequent fit cannot
silently reset either its displacement or strain allowance.
"""
import time
import numpy as np
from scipy import sparse
from .triangle_cross_candidates import TriangleCrossIndex
from .section_pair_separation import separated_triangle_axes, triangle_axis_sparse_rows


def reference_preservation_rows(first, second, reference_first, reference_second,
                                first_station, second_station, station_count, *,
                                maximum_motion_m, self_pair=False,
                                released_first_faces=None, released_second_faces=None,
                                deadline=None):
    for current, reference in [(first, reference_first), (second, reference_second)]:
        assert current.frame == reference.frame
        np.testing.assert_array_equal(current.faces, reference.faces)
        assert np.linalg.norm(current.vertices-reference.vertices, axis=1).max() <= maximum_motion_m+2e-9
    a, b = np.asarray(first_station), np.asarray(second_station)
    moving = np.flatnonzero(np.any(a[first.faces] >= 0, axis=1))
    blocks, bounds = [], []
    tested = retained = overlapping = released = 0
    # Both current and proposed controls are inside their original motion ball.
    padding = 2*maximum_motion_m*(2 if np.any(b >= 0) else 1)
    for local in TriangleCrossIndex(second.triangles).batches(
            first.triangles[moving], padding_m=padding, batch_size=128,
            maximum_proposals=2000000):
        if deadline is not None:
            assert time.monotonic() < deadline, 'Finite reference-interval stage'
        pair = np.c_[moving[local[:, 0]], local[:, 1]]
        if self_pair:
            reciprocal = np.isin(pair[:, 1], moving)
            pair = pair[~reciprocal | (pair[:, 0] < pair[:, 1])]
            shared = np.any(np.all(first.triangles[pair[:, 0]][:, :, None] ==
                                  second.triangles[pair[:, 1]][:, None], axis=3), axis=(1, 2))
            pair = pair[~shared]
        # A crossed repair patch may leave its present interpenetration. Every
        # currently healthy face on either actor retains its positive intervals.
        if released_first_faces is not None:
            both = np.isin(pair[:, 0], released_first_faces)
            if released_second_faces is not None:
                both &= np.isin(pair[:, 1], released_second_faces)
            released += int(both.sum())
            pair = pair[~both]
        if not len(pair):
            continue
        tested += len(pair)
        axes, gaps = separated_triangle_axes(first.triangles[pair[:, 0]], second.triangles[pair[:, 1]])
        keep = np.isfinite(gaps) & (gaps > 1e-12)
        overlapping += int((~keep).sum())
        retained += int(keep.sum())
        if not keep.any():
            continue
        A, rhs, _ = triangle_axis_sparse_rows(
            reference_first, reference_second, a, b, pair[keep], axes[keep],
            station_count=station_count, gap_m=np.minimum(2e-8, gaps[keep]*.5))
        live = np.asarray(A.getnnz(1)).ravel() > 0
        assert np.all(rhs[~live] <= 1e-10)
        blocks.append(A[live]); bounds.append(rhs[live])
    A = sparse.vstack(blocks, format='csr') if blocks else sparse.csr_matrix((0, 3*station_count))
    rhs = np.concatenate(bounds) if bounds else np.empty(0)
    return A, rhs, dict(complete=True, sourceMovingFaceIds=moving.tolist(),
                       testedPairs=tested, preservedIntervals=retained,
                       preexistingOverlappingOrTouchingPairs=overlapping,
                       releasedRepairPatchPairs=released, rows=len(rhs),
                       immutableReferenceCoordinates=True, reviewedPoseAxes=True,
                       wholeAnatomicalDisjointnessEstablished=False)
