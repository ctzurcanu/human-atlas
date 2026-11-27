"""Provisional corresponding inner envelope; the supplied exterior never moves.

This is a geometric construction, not a measured thyroid capsule thickness or
anatomical acceptance. Supply local float32 controls and an explicit complete
geometry gate. The two interfaces reuse one face table and any source charts;
the helper creates no region, texture, material, opening, or runtime record.
"""
from __future__ import annotations

from itertools import combinations
import numpy as np
from scipy.optimize import linprog


class CapsuleFieldError(ValueError):
    """Construction failed; diagnostics retain the offending control or trial."""

    def __init__(self, message, diagnostics=None):
        super().__init__(message)
        self.diagnostics = diagnostics or {}


def _require(condition, message):
    if not condition:
        raise CapsuleFieldError(message)


def incident_normal_field(outer, faces, target_distance, minimum_progress=.05, *, allow_feasible_lp_fallback=False):
    """Nearest nominal inward offset subject to all incident normal inequalities.

    All optimization coordinates are dimensionless, so small metre distances
    never become a solver feasibility tolerance. A bounded LP establishes a
    strictly feasible local cone only when the nominal direction needs repair.
    The unique Euclidean projection is a three-variable convex QP. Enumerating
    independent active sets of at most three constraints avoids a poorly scaled
    nonlinear solve. Each candidate is rechecked against every constraint.

    ``minimum_progress`` is an angular/numerical construction reserve, not a
    thickness measurement. Narrow feasible cones use a smaller explicit reserve.
    No feasible cone, ill-conditioned projection, or excessive departure fails
    explicitly rather than silently retaining an outward/collapsed offset.
    """
    outer = np.asarray(outer)
    faces = np.asarray(faces)
    _require(outer.ndim == 2 and outer.shape[1] == 3 and len(outer) >= 4
             and np.isfinite(outer).all(), 'Finite outer controls [n,3] required')
    _require(faces.ndim == 2 and faces.shape[1] == 3 and len(faces)
             and np.issubdtype(faces.dtype, np.integer)
             and faces.min() >= 0 and faces.max() < len(outer),
             'Integer triangles indexing the outer controls required')
    _require(np.isfinite(target_distance) and target_distance > 0,
             'Explicit positive target distance in metres required')
    _require(np.isfinite(minimum_progress) and 0 < minimum_progress < 1,
             'Positive dimensionless normal progress below one required')
    triangles = outer[faces].astype(float)
    area = np.cross(triangles[:, 1]-triangles[:, 0],
                    triangles[:, 2]-triangles[:, 0])
    length = np.linalg.norm(area, axis=1)
    _require(np.all(length > 0), 'Degenerate outer faces cannot define cones')
    normals = area/length[:, None]
    summed = np.zeros(outer.shape, dtype=float)
    incidents = [[] for _ in outer]
    for corner in range(3):
        np.add.at(summed, faces[:, corner], area)
        for face, control in enumerate(faces[:, corner]):
            incidents[int(control)].append(face)
    _require(all(incidents), 'Unused outer controls cannot define cones')
    summed_length = np.linalg.norm(summed, axis=1)
    nominal = np.divide(-summed, summed_length[:, None],
                        out=np.zeros_like(summed),
                        where=summed_length[:, None] > 0)
    directions = nominal.copy()
    reserves = np.full(len(outer), minimum_progress)
    corrections = []
    worst_condition = 1.
    for control, face_ids in enumerate(incidents):
        n = normals[face_ids]
        if summed_length[control] > 0 and np.all(n@nominal[control] <= -minimum_progress):
            continue
        # Maximize the inward margin in a compact cube. Strict feasibility of
        # the cone is equivalent to a positive optimum; actual row products,
        # rather than the LP status alone, establish a feasible direction.
        lp = linprog(np.array([0., 0., 0., -1.]),
                     A_ub=np.c_[n, np.ones(len(n))], b_ub=np.zeros(len(n)),
                     bounds=[(-1., 1.)]*3+[(0., 1.)], method='highs',
                     options={'primal_feasibility_tolerance': 1e-9,
                              'dual_feasibility_tolerance': 1e-9})
        if not lp.success or not np.all(n@lp.x[:3] < 0):
            raise CapsuleFieldError('No strictly inward incident-normal cone',
                                    {'control': control, 'incidentFaces': face_ids,
                                     'lpStatus': int(lp.status),
                                     'lpMessage': lp.message,
                                     'normals': n.tolist()})
        seed = lp.x[:3]/np.linalg.norm(lp.x[:3])
        cone_margin = float(np.min(-n@seed))
        reserve = min(minimum_progress, cone_margin*.25)
        reserves[control] = reserve
        # Slightly tighten the *construction* bound before solving; the strict
        # original bound is tested without a permissive residual allowance.
        bound = reserve*(1+256*np.finfo(float).eps)
        q = nominal[control]
        candidates = []
        if np.all(n@q <= -reserve):
            candidates.append((0., q, (), 1.))
        for count in range(1, min(3, len(n))+1):
            for chosen in combinations(range(len(n)), count):
                active = n[list(chosen)]
                gram = active@active.T
                condition = float(np.linalg.cond(gram))
                if not np.isfinite(condition) or condition > 1e12:
                    continue
                multiplier = np.linalg.solve(gram, active@q+bound)
                if np.any(multiplier < 0):
                    continue
                candidate = q-active.T@multiplier
                if np.all(n@candidate <= -reserve):
                    candidates.append((float(np.sum((candidate-q)**2)),
                                       candidate, chosen, condition))
        used_fallback=False
        if not candidates and allow_feasible_lp_fallback and np.all(n@seed<=-reserve):
            # A strictly feasible unit LP direction is independently checked
            # above and has four times the requested cone reserve. Ill-
            # conditioned active-set projection arithmetic need not prevent a
            # bounded feasible field. Callers must explicitly opt into this
            # less-near-nominal construction; complete field gates still apply.
            candidates.append((float(np.sum((seed-q)**2)),seed,(),1.))
            used_fallback=True
        if not candidates:
            raise CapsuleFieldError('Conditioned incident-normal projection failed',
                                    {'control': control, 'incidentFaces': face_ids,
                                     'coneMargin': cone_margin,
                                     'minimumProgress': reserve})
        departure, candidate, active, condition = min(candidates, key=lambda x:x[0])
        if np.linalg.norm(candidate) > 2:
            raise CapsuleFieldError('Incident-normal projection exceeds bounded departure',
                                    {'control': control, 'offsetTargetRatio':
                                     float(np.linalg.norm(candidate))})
        directions[control] = candidate
        worst_condition = max(worst_condition, condition)
        corrections.append({'control': control, 'incidentFaces': face_ids,
                            'feasibleLpFallback': used_fallback,
                            'activeFaces': [face_ids[k] for k in active],
                            'coneMargin': cone_margin,
                            'minimumProgress': reserve,
                            'projectionSquaredDeparture': departure,
                            'activeGramCondition': condition})
    dot = np.einsum('fvi,fi->fv', directions[faces], normals)
    _require(np.all(dot < 0), 'Projected field is not strictly inward')
    return directions*target_distance, {
        'targetDistanceMeters': float(target_distance),
        'clinicallyMeasuredThickness': False,
        'normalProgressReserve': float(minimum_progress),
        'correctedControls': corrections,
        'correctedControlCount': len(corrections),
        'maximumIncidentNormalDotDirection': float(dot.max()),
        'minimumControlProgressReserve': float(reserves.min()),
        'maximumActiveGramCondition': worst_condition,
        'maximumDepartureFromNominalMeters':
            float(np.linalg.norm((directions-nominal)*target_distance, axis=1).max()),
    }


def build_inner_envelope(outer_local, faces, target_distance, geometry_gate, *,
                         minimum_progress=.05, minimum_target_fraction=.0625):
    """Build one stored inner interface, with explicit bounded backtracking.

    ``geometry_gate(outer_local, inner_local, faces)`` must validate complete
    outer/inner/shell topology, self/cross intersections and continuous paired
    prism extrema on the actual float32 arrays, returning ``passesGeometryGates``.
    It is required: local cones alone do not prove global injectivity. First try
    the nominal distance; reduce the entire projected field by powers of two if
    necessary. This is a conservative line search, not a global thickness optimum.
    """
    outer = np.asarray(outer_local)
    _require(outer.dtype == np.dtype('<f4'),
             'Supply actual local float32 outer controls; implicit rounding is forbidden')
    _require(callable(geometry_gate), 'Complete saved-float32 geometry gate required')
    _require(np.isfinite(minimum_target_fraction) and 0 < minimum_target_fraction <= 1,
             'Explicit positive minimum target fraction below or equal to one required')
    field, report = incident_normal_field(outer, faces, target_distance, minimum_progress)
    before = outer.copy()
    attempts = []
    fraction = 1.
    while fraction >= minimum_target_fraction:
        inner = (outer.astype(float)+field*fraction).astype('<f4')
        evidence = geometry_gate(outer, inner, faces)
        attempts.append({'targetFraction': fraction, 'evidence': evidence})
        _require(np.array_equal(outer, before), 'Geometry gate changed fixed exterior')
        if evidence.get('passesGeometryGates') is True:
            report.update({'selectedTargetFraction': fraction,
                           'selectedTargetDistanceMeters': float(target_distance*fraction),
                           'attempts': attempts,
                           'maximumFloat32OffsetRoundingMeters':
                           float(np.linalg.norm(inner.astype(float)-outer.astype(float)
                                                -field*fraction, axis=1).max()),
                           'minimumStoredOffsetMeters':
                           float(np.linalg.norm(inner.astype(float)-outer.astype(float), axis=1).min()),
                           'maximumStoredOffsetMeters':
                           float(np.linalg.norm(inner.astype(float)-outer.astype(float), axis=1).max()),
                           'outerControlsExact': True,
                           'installed': False, 'anatomicallyAccepted': False})
            return inner, report
        fraction *= .5
    report.update({'attempts': attempts, 'minimumTargetFraction': minimum_target_fraction,
                   'installed': False, 'anatomicallyAccepted': False})
    raise CapsuleFieldError('No finite inner envelope within bounded target fractions', report)
