"""Certified distance bounds between literal piecewise-linear source contours.

Unlike vertex matching, this metric is independent of unequal sampling along an
edge. It does not establish anatomical identity, wall ownership or patency.
"""
from __future__ import annotations

import heapq
import numpy as np


def directed_polyline_distance(source, target, *, source_closed=True,
                               target_closed=True, tolerance=1e-7,
                               max_evaluations=10000):
    """Bound max(source -> target) over complete segments, with a real witness.

    Distance to one convex segment is convex along a source segment. Therefore
    min_j max(d(a, segment_j), d(b, segment_j)) bounds every interior point's
    distance to the union. Subdivision closes this bound without assuming that
    the maximum occurs at an input vertex. Units follow the input coordinates.
    """
    s, t = np.asarray(source, dtype=float), np.asarray(target, dtype=float)
    if (s.ndim != 2 or t.ndim != 2 or s.shape[1] != t.shape[1] or s.shape[1] < 1
            or min(len(s), len(t)) < 2 or not np.isfinite(s).all()
            or not np.isfinite(t).all() or not np.isfinite(tolerance) or tolerance <= 0
            or max_evaluations < len(s)):
        raise ValueError("Finite matching-dimensional polylines and positive limits required")
    ta = t if target_closed else t[:-1]
    tb = np.roll(t, -1, axis=0) if target_closed else t[1:]
    edge = tb - ta
    length2 = np.einsum("ij,ij->i", edge, edge)
    evaluations = 0
    lower = -1.0
    witness = None

    def evaluate(point):
        nonlocal evaluations, lower, witness
        if evaluations >= max_evaluations:
            raise RuntimeError("Polyline distance evaluation cap reached; no certified result")
        evaluations += 1
        numerator = np.einsum("ij,ij->i", point-ta, edge)
        fraction = np.divide(numerator, length2, out=np.zeros_like(numerator), where=length2 > 0)
        closest = ta + np.clip(fraction, 0, 1)[:, None] * edge
        distance = np.linalg.norm(point - closest, axis=1)
        j = int(distance.argmin())
        if distance[j] > lower:
            lower = float(distance[j])
            witness = dict(sourcePoint=point.tolist(), targetPoint=closest[j].tolist(),
                           targetSegment=j, distance=lower)
        return distance

    distances = [evaluate(point) for point in s]
    heap = []
    serial = 0

    def push(a, b, da, db):
        nonlocal serial
        # A small outward rounding allowance protects the numeric certificate.
        bound = float(np.min(np.maximum(da, db)))
        coordinate_scale = max(float(np.max(np.abs(q))) for q in (a, b, ta, tb))
        allowance = np.finfo(float).eps * 32 * max(1., coordinate_scale)
        serial += 1
        heapq.heappush(heap, (-(bound+allowance), serial, a, b, da, db))

    for i in range(len(s) if source_closed else len(s)-1):
        j = (i+1) % len(s)
        push(s[i], s[j], distances[i], distances[j])
    while -heap[0][0] > lower + tolerance:
        _, _, a, b, da, db = heapq.heappop(heap)
        mid = (a+b)/2
        dm = evaluate(mid)
        push(a, mid, da, dm)
        push(mid, b, dm, db)
    return dict(lower=lower, upper=max(lower, -heap[0][0]), tolerance=tolerance,
                witness=witness, evaluations=evaluations, completeSegments=True)


def polyline_hausdorff(first, second, **kwargs):
    """Symmetric complete-segment distance, retaining both directed witnesses."""
    forward = directed_polyline_distance(first, second, **kwargs)
    reverse_args = dict(kwargs)
    if "source_closed" in kwargs or "target_closed" in kwargs:
        reverse_args["source_closed"] = kwargs.get("target_closed", True)
        reverse_args["target_closed"] = kwargs.get("source_closed", True)
    backward = directed_polyline_distance(second, first, **reverse_args)
    return dict(lower=max(forward["lower"], backward["lower"]),
                upper=max(forward["upper"], backward["upper"]),
                forward=forward, backward=backward, completeSegments=True)
