#!/usr/bin/env python3
"""One source-bound ellipsoidal mechanism correction; no candidate geometry."""
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.compact_positive_fields import compose_compact_radial_maps


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1'
    previous = base / 'lower-limb-joints-v1'
    out = base / 'lower-limb-joints-v2'
    out.mkdir(exist_ok=False)
    baseline_path = base / 'body-envelope-fit-v2/contract.json'
    old_contract_path = previous / 'contract.json'
    guard_path = previous / 'localized-skin-regression-v1/report.json'
    comparison_path = previous / 'trial-v1-audit/comparison.json'
    contract = json.loads(baseline_path.read_text())
    old = json.loads(old_contract_path.read_text())
    guard = json.loads(guard_path.read_text())
    comparison = json.loads(comparison_path.read_text())
    paths = [baseline_path, old_contract_path, guard_path, comparison_path, previous / 'popliteal-section-guides.npz', previous / 'vessel-sample-ownership.npz']
    bindings = {str(path): sha(path) for path in paths}
    assert sha(baseline_path) == old['post_field_constraints']['baselineContractSha256']
    fields = [{'kind': 'ellipsoidal-translation', 'center': field['center'], 'support_axes_m': [.08, .16, .024], 'core_fraction': .2, 'translation_m': field['translation_m']} for field in old['post_fields']]
    original_guides = np.load(previous / 'popliteal-section-guides.npz')['actual_source_section_points_m']
    centres = np.asarray([field['center'] for field in fields])
    posterior_failure_controls = centres.copy()
    posterior_failure_controls[:, 2] = -.0722
    outside_controls = np.vstack((posterior_failure_controls, centres + [0, 0, -.025], centres + [0, 0, .025], np.array([[0., 1.4, 0.], [0., 1., 0.], [0., 0., 0.]])))
    points = np.vstack((original_guides, centres, outside_controls))
    mapped, jacobian, proof = compose_compact_radial_maps(points, fields)
    assert np.array_equal(mapped[-len(outside_controls):], outside_controls), 'Outside ellipsoid controls changed'
    assert np.array_equal(jacobian[-len(outside_controls):], np.broadcast_to(np.eye(3), (len(outside_controls), 3, 3))), 'Outside ellipsoid Jacobians changed'
    assert np.linalg.det(jacobian).min() >= proof['analytic_determinant_lower_bound'] - 1e-12
    for i, field in enumerate(fields):
        # The contralateral support does not reach the other centre.
        assert np.array_equal(mapped[len(original_guides) + i], centres[i] + field['translation_m'])
    np.savez_compressed(out / 'field-control-checks.npz', input_points_m=points, mapped_points_m=mapped, jacobian=jacobian, determinant=np.linalg.det(jacobian), actual_source_guide_point_count=np.array(len(original_guides)), outside_control_start=np.array(len(original_guides) + len(centres)))
    contract['schema'] = 'biodigital-post-body-common-knee-trial/v2'
    contract['stagePurpose'] = 'One ellipsoidal support mechanism correction after measured spherical-field superficial-saphenous skin failure; no change of centres or deep-corridor translations'
    contract['post_fields'] = fields
    contract['post_field_proof'] = proof
    contract['analytic_proof'] = {'analytic_determinant_lower_bound': float(contract['analytic_proof']['analytic_determinant_lower_bound'] * proof['analytic_determinant_lower_bound']), 'composition': 'Jpost @ Jbody, after unchanged Gnew; continuous lower bounds multiply', 'discrete_and_anatomical_checks_pending': True}
    contract['post_field_constraints'] = {**old['post_field_constraints'], 'baselineContract': str(baseline_path), 'baselineContractSha256': bindings[str(baseline_path)], 'candidateApplied': False, 'oneTrialOnly': True, 'supportShapeCorrection': {'previousSphericalContract': str(old_contract_path), 'previousSphericalContractSha256': bindings[str(old_contract_path)], 'skinGuardEvidence': str(guard_path), 'skinGuardEvidenceSha256': bindings[str(guard_path)], 'priorBoneComparison': str(comparison_path), 'priorBoneComparisonSha256': bindings[str(comparison_path)], 'observedFailure': '14 newly outside small-saphenous/extension vertices, max outside distance1.407mm; broad spherical core moved superficial posterior band', 'mechanism': 'Same centres and translations; 24mm AP support gives literal identity beyond the posterior ellipsoid limit. No source-family mask and no isolated vein displacement.', 'worldSupportAxes_m': [.08, .16, .024], 'normalizedCoreFraction': .2, 'posteriorIdentityPlanesZ_m': [float(field['center'][2] - .024) for field in fields], 'failedBandRepresentativeZ_m': -.0722, 'failedBandControlIdentityVerified': True}, 'acceptanceCriteria': ['No new zero-area/inverted retained source triangles; preserve exact source bank/face/vertex/UV ancestry and normal transport', 'Frozen-source-face major knee crossings improve versus body-v2; retain nutrient/foraminal uncertainty rather than force clearing branches', 'Actual skin guard shows zero new baseline-inside to outside vertices, including small-saphenous/extension owners', 'Review full named major-trunk/genicular triangle-area strain and strict skin triangle crossings before acceptance', 'Central torso and all points outside the ellipsoids preserve current body-v2 coordinates literally']}
    contract['post_field_constraints']['sideMeasurements'] = [{key: value for key, value in row.items() if key in ['side', 'measuredFullPosteriorShiftForSectionBoneZMargin_m', 'mostRestrictiveSkinZShiftLowerBound_m', 'guidePointCount', 'fullShiftIsNotAnAnatomicalTarget', 'translation_m']} for row in old['post_field_constraints']['sideMeasurements']]
    contract['required_followups'] = ['Run this single support-correction candidate only; no centre-fit/optimizer loop', 'Reuse server/audit-biodigital-lower-limb-trial-v1.py with --candidate-stage lower-limb-joints-v2 --output trial-v2-audit', 'Repeat the exact localized skin regression guard; actual posterior superficial controls must remain fixed', 'Keep rendered skin/skeleton and all nonblood anatomy unchanged'] + contract['required_followups']
    (out / 'contract.json').write_text(json.dumps(contract, indent=2) + '\n')
    report = {'schema': 'ellipsoidal-knee-mechanism-correction/v1', 'status': 'contract-ready; candidate-not-applied', 'inputBindings': bindings, 'fields': fields, 'analyticPostProof': proof, 'fullBodyPostAnalyticBound': contract['analytic_proof']['analytic_determinant_lower_bound'], 'sourceGuideControlCount': len(original_guides), 'minimumControlDeterminant': float(np.linalg.det(jacobian).min()), 'outsideControlsLiteralIdentity': True, 'sameCentresAndTranslations': True, 'sourceFaceBenchmarkUnchanged': bindings[str(previous / 'vessel-sample-ownership.npz')] == old['post_field_constraints']['sameSourceFaceBenchmarkSha256'], 'renderedGeometryModified': False, 'newSpatialQuery': False, 'optimizerUsed': False, 'elapsedSeconds': time.monotonic() - began}
    (out / 'preparation.json').write_text(json.dumps(report, indent=2) + '\n')
    assert all(sha(path) == value for path, value in bindings.items()), 'Bound input changed during contract preparation'
    print(json.dumps({'contract': str(out / 'contract.json'), 'postBounds': [x['analytic_determinant_lower_bound'] for x in proof['fields']], 'postProductBound': proof['analytic_determinant_lower_bound'], 'fullBodyPostBound': report['fullBodyPostAnalyticBound'], 'minimumControlDeterminant': report['minimumControlDeterminant'], 'outsideIdentity': True, 'elapsedSeconds': time.monotonic() - began}))


if __name__ == '__main__':
    main()
