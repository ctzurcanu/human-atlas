#!/usr/bin/env python3
"""Measured combined popliteal corridor; one final common-field contract."""
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def section(triangles, y):
    points, faces, corners, alphas = [], [], [], []
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        a, b = triangles[:, i], triangles[:, j]
        take = np.flatnonzero(((a[:, 1] < y) & (b[:, 1] > y)) | ((a[:, 1] > y) & (b[:, 1] < y)))
        if len(take):
            t = (y - a[take, 1]) / (b[take, 1] - a[take, 1])
            points.extend(a[take] + t[:, None] * (b[take] - a[take]))
            faces.extend(take)
            corners.extend([[i, j]] * len(take))
            alphas.extend(t)
    return np.asarray(points).reshape((-1, 3)), np.asarray(faces, dtype=np.int64), np.asarray(corners, dtype=np.int8).reshape((-1, 2)), np.asarray(alphas)


def enclosure(curve, ys, triangles, metas):
    """Max X/Z deviation over whole triangles, not just their vertices.

    At each fixed Y, linear triangle X/Z extrema lie on edges. Along an edge,
    subtract the PCHIP cubic and check endpoints, knots and derivative roots.
    """
    maximum, witnesses = np.zeros(2), [None, None]
    for fi, triangle in enumerate(triangles):
        for i, j in [(0, 1), (1, 2), (2, 0)]:
            a, b = triangle[i], triangle[j]
            lo, hi = sorted((a[1], b[1]))
            for axis_index, axis in enumerate([0, 2]):
                candidates = [float(lo), float(hi)]
                if hi > lo:
                    slope = (b[axis] - a[axis]) / (b[1] - a[1])
                    cuts = np.unique(np.r_[lo, ys[(ys > lo) & (ys < hi)], hi])
                    candidates.extend(cuts.tolist())
                    for left, right in zip(cuts[:-1], cuts[1:]):
                        mid = (left + right) / 2
                        if mid <= ys[0] or mid >= ys[-1]:
                            continue
                        k = np.clip(np.searchsorted(ys, mid, side='right') - 1, 0, len(ys) - 2)
                        coeff = curve[axis_index].c[:, k]
                        roots = np.roots([-3 * coeff[0], -2 * coeff[1], slope - coeff[2]])
                        candidates.extend(float(root.real + ys[k]) for root in roots if abs(root.imag) < 1e-10 and left < root.real + ys[k] < right)
                for y in candidates:
                    for alpha in ([0., 1.] if hi == lo else [(y - a[1]) / (b[1] - a[1])]):
                        point = a + alpha * (b - a)
                        clipped = np.clip(y, ys[0], ys[-1])
                        cc = float(curve[axis_index](clipped) + curve[axis_index].derivative()(clipped) * (y - clipped))
                        offset = abs(float(point[axis]) - cc)
                        if offset > maximum[axis_index]:
                            maximum[axis_index] = offset
                            witnesses[axis_index] = {**metas[fi], 'edgeCorners': [i, j], 'edgeAlpha': float(alpha), 'point_m': point.tolist(), 'curveCoordinate_m': cc, 'coordinate': 'X' if axis == 0 else 'Z', 'absoluteOffset_m': offset}
    return maximum, witnesses


def main():
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1'
    baseline = base / 'body-envelope-fit-v2'
    candidate = baseline / 'candidate-v1'
    first, previous = base / 'lower-limb-joints-v1', base / 'lower-limb-joints-v2'
    output = base / 'lower-limb-joints-v3'
    output.mkdir(exist_ok=False)
    paths = [baseline / 'contract.json', first / 'constraints.json', first / 'contract.json', first / 'vessel-sample-ownership.npz', previous / 'contract.json', previous / 'trial-v2-audit/comparison.json', previous / 'localized-skin-regression-v1/report.json', candidate / 'instance-ownership.json']
    bindings = {str(path): sha(path) for path in paths}
    contract = json.loads(paths[0].read_text())
    guides = json.loads(paths[1].read_text())
    old = json.loads(paths[2].read_text())
    rejected_comparison = json.loads(paths[5].read_text())
    rows = json.loads(paths[7].read_text())['instances']
    popliteal = {}
    for mesh_name in ['arterial', 'venous']:
        mesh_path, ancestry_path = candidate / (mesh_name + '.npz'), candidate / (mesh_name + '-ancestry.npz')
        bindings[str(mesh_path)], bindings[str(ancestry_path)] = sha(mesh_path), sha(ancestry_path)
        z, az = np.load(mesh_path), np.load(ancestry_path)
        vertices, faces = z['vertices'], z['faces']
        source_faces = az['candidate_face_source_face']
        for row in rows:
            if row['candidateMesh'] != mesh_name + '.npz' or row['name'] not in ['Right popliteal artery', 'Left popliteal artery', 'Right popliteal vein', 'Left popliteal vein']:
                continue
            start, count = row['candidateFaceStart'], row['candidateFaceCount']
            t = vertices[faces[start:start + count]]
            metas = [{'sourceInstanceId': row['sourceInstanceId'], 'name': row['name'], 'candidateMesh': mesh_name, 'candidateFace': start + i, 'sourceFace': int(source_faces[start + i])} for i in range(count)]
            popliteal[row['name']] = (t, metas)
        del z, az, vertices, faces, source_faces
    fields, side_reports, all_points, point_rows = [], [], [], []
    per_field_bounds = []
    for side, old_field in zip(['Right', 'Left'], old['post_fields']):
        vessels = [popliteal[side + ' popliteal artery'], popliteal[side + ' popliteal vein']]
        triangles = np.concatenate([item[0] for item in vessels])
        metas = [row for item in vessels for row in item[1]]
        saved_ys = sorted(set(q['sectionY_m'] for q in guides['poplitealActualSectionGuides'] if q['vessel'].startswith(side)))
        minimum, maximum = float(triangles[:, :, 1].min()), float(triangles[:, :, 1].max())
        # Additional actual near-end sections cover the full source course;
        # these are geometric field controls, never anatomical mouth points.
        ys = np.unique(np.r_[saved_ys, minimum + .0001, maximum - .0001])
        stations, station_reports = [], []
        for station_index, y in enumerate(ys):
            points, fs, corners, alphas = section(triangles, float(y))
            assert len(points), 'No actual source triangle section at required station'
            bounds = np.array([points.min(axis=0), points.max(axis=0)])
            centre = bounds.mean(axis=0)
            stations.append({'y_m': float(y), 'x_m': float(centre[0]), 'z_m': float(centre[2])})
            station_reports.append({'stationIndex': station_index, 'y_m': float(y), 'pointCount': len(points), 'combinedArterialVenousBounds_m': bounds.tolist(), 'midpointIsGeometricGuideNotAnatomicalLandmark': True})
            all_points.extend(points)
            point_rows.extend({'side': side, 'stationIndex': station_index, **metas[int(fi)], 'edgeCorners': corners[k].tolist(), 'edgeAlpha': float(alphas[k])} for k, fi in enumerate(fs))
        curves = [PchipInterpolator(ys, [s[key] for s in stations], extrapolate=False) for key in ['x_m', 'z_m']]
        offsets, extrema = enclosure(curves, ys, triangles, metas)
        core = np.array([.015, .055, .010])
        support = np.array([.08, .16, .040])
        for measured, dimension, transition in [(offsets[0], 0, .05), (offsets[1], 2, .03)]:
            if measured + 1e-9 >= core[dimension]:
                core[dimension] = np.ceil((measured + .001) * 1000) / 1000
                support[dimension] = max(support[dimension], core[dimension] + transition)
        d = np.asarray(old_field['translation_m'])
        lower = 1 - 2.109375 * abs(d[2]) / (support[2] - core[2])
        assert lower > 0 and offsets[0] + 1e-9 < core[0] and offsets[1] + 1e-9 < core[2]
        field = {'kind': 'separable-corridor-translation', 'centerline_stations': stations, 'core_axes_m': core.tolist(), 'support_axes_m': support.tolist(), 'center_y_m': old_field['center'][1], 'translation_m': old_field['translation_m']}
        fields.append(field)
        per_field_bounds.append(lower)
        side_reports.append({'side': side, 'fullPoplitealSourceFaces': len(triangles), 'fullSourceYBounds_m': [minimum, maximum], 'stationCount': len(stations), 'stations': station_reports, 'maximumCompleteTriangleOffsetsXZ_m': offsets.tolist(), 'offsetExtremaWitnesses': extrema, 'completePoplitealTriangleEnvelopeInsideFlatCoreXZ': True, 'flatCoreMarginsXZ_m': [float(core[0] - offsets[0]), float(core[2] - offsets[1])], 'originalRequestedCoreAxes_m': [.015, .055, .010], 'chosenCoreAxes_m': core.tolist(), 'chosenSupportAxes_m': support.tolist(), 'coreExpandedOnceFromMeasuredEnvelope': bool(core[0] != .015 or core[2] != .010), 'analyticDeterminantLowerBound': lower})
    proof = {'kind': 'ordered-composition', 'field_count': 2, 'analytic_determinant_lower_bound': float(np.prod(per_field_bounds)), 'fields': [{'kind': f['kind'], 'analytic_determinant_lower_bound': lower, 'fixed_outside_support': True} for f, lower in zip(fields, per_field_bounds)], 'pureZDeterminant': '1 + dZ * partialWeight/partialZ; X/Y gradients do not enter determinant', 'discrete_mesh_checks_pending': True}
    contract['schema'] = 'biodigital-post-body-common-knee-trial/v3'
    contract['stagePurpose'] = 'One final source-course-following separable posterior corridor translation after measured v2 transverse surface compression; no centre search'
    contract['post_fields'] = fields
    contract['post_field_proof'] = proof
    contract['analytic_proof'] = {'analytic_determinant_lower_bound': float(contract['analytic_proof']['analytic_determinant_lower_bound'] * proof['analytic_determinant_lower_bound']), 'composition': 'Jpost @ Jbody after unchanged Gnew; continuous positive bounds multiply', 'discrete_and_anatomical_checks_pending': True}
    contract['post_field_constraints'] = {**old['post_field_constraints'], 'candidateApplied': False, 'oneFinalTrialOnly': True, 'measuredCorridorEvidence': str(output / 'corridor-evidence.json'), 'sourceCurvePolicy': 'PCHIP X(Y),Z(Y) from measured combined section bounds midpoints; C1 linear endpoint-tangent extension: value(clippedY)+slope(clippedY)*(Y-clippedY); world axes unchanged', 'completeTriangleFlatCoreEnclosure': side_reports, 'v2RejectionEvidence': {'contractSha256': bindings[str(paths[4])], 'boneComparisonSha256': bindings[str(paths[5])], 'skinGuardSha256': bindings[str(paths[6])], 'reason': 'No new skin exits but major popliteal surface triangle areas compressed below half; not installed', 'poplitealAreaStrain': [r for r in rejected_comparison['actualTriangleAreaStrain'] if r['name'] in popliteal]}, 'acceptanceCriteria': ['Preserve all1724 source instances, exact retained face/vertex/UV ancestry and transported normals', 'No new zero-area or inverted retained source triangles', 'Frozen-source-face main knee strict crossings improve versus body-v2', 'No major popliteal retained surface face area ratio below0.5 versus body-v2; report full area-ratio distributions for all major genicular/femoral/tibial routes', 'Zero new baseline-inside to outside vertices in the repeated actual-skin regression guard', 'Actual skin crossings, adaptive clearance and bone/lumen evidence remain unaccepted until separately checked', 'All skin/skeleton fixed; no isolated vessel displacement or opening inferred from centroid']}
    contract['required_followups'] = ['Run this one final source-course corridor candidate; no optimizer or field-parameter search', 'Reuse frozen source-face lower-limb audit and localized actual-skin before/after guard', 'Reject if deformation or skin checks fail; preserve numeric evidence'] + contract['required_followups']
    (output / 'contract.json').write_text(json.dumps(contract, indent=2) + '\n')
    points = np.asarray(all_points)
    np.savez_compressed(output / 'actual-combined-section-points.npz', points_m=points, side=np.asarray([r['side'] for r in point_rows]), station_index=np.asarray([r['stationIndex'] for r in point_rows]), source_instance_id=np.asarray([r['sourceInstanceId'] for r in point_rows]), candidate_mesh=np.asarray([r['candidateMesh'] for r in point_rows]), candidate_face=np.asarray([r['candidateFace'] for r in point_rows]), original_source_face=np.asarray([r['sourceFace'] for r in point_rows]), edge_corners=np.asarray([r['edgeCorners'] for r in point_rows]), edge_alpha=np.asarray([r['edgeAlpha'] for r in point_rows]))
    report = {'schema': 'measured-source-course-common-knee-corridor/v1', 'inputBindings': bindings, 'sideEvidence': side_reports, 'sourceSectionPointCount': len(points), 'interpolation': 'SciPy PCHIP with C1 linear endpoint-tangent extension; combined artery+vein section X/Z bounds midpoint', 'continuousTriangleEnclosureMethod': 'For every actual source triangle edge: cubic PCHIP subtraction, interval endpoints, station knots and real derivative roots; fixed-Y linear-triangle extrema lie on edges. Endpoint-extension pieces are linear and use both interval endpoints. Float64 with1nm margin, no samples substituted for edge extrema.', 'postFieldProof': proof, 'fullBodyPostAnalyticBound': contract['analytic_proof']['analytic_determinant_lower_bound'], 'geometryChanged': False, 'newSpatialQuery': False, 'optimizerUsed': False, 'elapsedSeconds': time.monotonic() - began}
    (output / 'corridor-evidence.json').write_text(json.dumps(report, indent=2) + '\n')
    assert all(sha(path) == value for path, value in bindings.items()), 'Bound evidence changed'
    print(json.dumps({'contract': str(output / 'contract.json'), 'sides': [{'side': r['side'], 'stations': r['stationCount'], 'faces': r['fullPoplitealSourceFaces'], 'offsetsXZ_mm': (np.asarray(r['maximumCompleteTriangleOffsetsXZ_m']) * 1000).tolist(), 'coreAxes_m': r['chosenCoreAxes_m'], 'supportAxes_m': r['chosenSupportAxes_m'], 'lowerBound': r['analyticDeterminantLowerBound']} for r in side_reports], 'postProductBound': proof['analytic_determinant_lower_bound'], 'fullBodyPostBound': report['fullBodyPostAnalyticBound'], 'elapsedSeconds': time.monotonic() - began}))


if __name__ == '__main__':
    main()
