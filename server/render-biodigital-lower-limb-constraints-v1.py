#!/usr/bin/env python3
"""Render saved literal triangles and write a small common-field recipe."""
import json
import hashlib
import time
from pathlib import Path
import argparse

import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--no-render', action='store_true')
    args = parser.parse_args()
    began = time.monotonic()
    base = Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
    out = base / 'lower-limb-joints-v1'
    evidence = out / 'constraints.json'
    r = json.loads(evidence.read_text())
    fields = r['firstCorrection']['fields']
    baseline_contract_path = base / 'body-envelope-fit-v2/contract.json'
    trial_contract = json.loads(baseline_contract_path.read_text())
    trial_contract['schema'] = 'biodigital-post-body-common-knee-trial/v1'
    trial_contract['stagePurpose'] = 'One damped bilateral posterior popliteal-family common-field trial after unchanged body-envelope-v2; fixed skeleton and skin'
    trial_contract['post_fields'] = fields
    trial_contract['post_field_proof'] = r['firstCorrection']['trialGuideChecks']['analyticProof']
    trial_contract['analytic_proof'] = {'analytic_determinant_lower_bound': float(trial_contract['analytic_proof']['analytic_determinant_lower_bound'] * trial_contract['post_field_proof']['analytic_determinant_lower_bound']), 'composition': 'Jpost @ Jbody; both after Gnew similarity; full Jacobian normal transport', 'discrete_and_anatomical_checks_pending': True}
    trial_contract['post_field_constraints'] = {'baselineContract': str(baseline_contract_path), 'baselineContractSha256': hashlib.sha256(baseline_contract_path.read_bytes()).hexdigest(), 'measuredEvidence': str(evidence), 'measuredEvidenceSha256': hashlib.sha256(evidence.read_bytes()).hexdigest(), 'sameSourceFaceBenchmark': str(out / 'vessel-sample-ownership.npz'), 'sameSourceFaceBenchmarkSha256': hashlib.sha256((out / 'vessel-sample-ownership.npz').read_bytes()).hexdigest(), 'fixedContextSha256': r['inputBindings'][str((base / 'regional-interface-v1/fixed-context/full.npz').resolve())], 'sideMeasurements': r['firstCorrection']['sideMeasurements'], 'sourceClosure': str(out / 'source-closure.json'), 'scope': 'All1724 original instances receive the same spatial maps; nutrient/perforating owners are preserved, not force-cleared or individually displaced', 'candidateApplied': False, 'oneTrialOnly': True, 'acceptanceCriteria': ['No new zero-area/inverted retained source face; all source ancestry/UV banks and normal transport preserved', 'Same-source-face major knee crossing count improves; no claim that all nutrient/canal intersections are wrong', 'Full named major-trunk triangle area ratios and skin crossings/inside/adaptive-clearance replay reviewed before acceptance', 'Central torso and all points outside both finite knee supports preserve current v2 coordinates literally']}
    trial_contract['required_followups'] = ['Run the one named knee trial; do not optimize repeatedly', 'Replay saved exact lower-limb source-face samples, full major-family surface area ratios, actual skin and fixed bone evidence'] + trial_contract['required_followups']
    (out / 'contract.json').write_text(json.dumps(trial_contract, indent=2) + '\n')
    recipe = {'schema': 'post-body-envelope-common-knee-trial-recipe/v1', 'baselineCandidate': r['baselineCandidate'], 'baselineContractSha256': r['inputBindings'][str((base / 'body-envelope-fit-v2/contract.json').resolve())], 'fixedContextSha256': r['inputBindings'][str((base / 'regional-interface-v1/fixed-context/full.npz').resolve())], 'evidenceSha256': hashlib.sha256(evidence.read_bytes()).hexdigest(), 'sourceClosure': str(out / 'source-closure.json'), 'applicationOrder': ['source-v2 original world instances', 'Gnew registration', 'unchanged body-envelope-v2 common map', 'the following two compact common fields, with complete chain-rule Jacobian and inverse-transpose normals'], 'movingScope': 'All1724 source-v2 instances receive the same ordered common map; effective influence is confined to finite knee supports', 'fixedScope': 'All skin/skeleton/current nonblood anatomy', 'fields': fields, 'analyticProof': r['firstCorrection']['trialGuideChecks'].get('analyticProof'), 'measuredTrial': r['firstCorrection']['sideMeasurements'], 'tolerances': {'fieldAnalyticJacobianLowerBound': .738501832560733, 'newZeroAreaFaces': 0, 'newInvertedFaces': 0, 'commonOutsideSupportCoordinates': 'literal identity', 'sourceFaceAndVertexAncestry': 'exact retained IDs', 'requiredPostApplicationChecks': ['same saved strict bone witness triangles with full candidate transformed faces', 'full major popliteal/genicular/femoral/tibial family surface strain', 'actual skin triangle crossings plus multi-ray inside and adaptive local clearance', 'communicating route/mouth and lumen geometry review']}, 'qualifications': ['Section points are actual surface/plane intersections; bone-ray hits are actual fixed source triangles. Neither is labelled an anatomical mouth/foramen.', 'One damped trial only; residual section bone-corridor shift remains about7mm and is deliberately not driven tozero.', '274section-point actual skin distances are finite evidence, not full-surface acceptance.', 'Anterior tibial interosseous routes and nutrient openings need their own actual geometry constraints; no exemption from cortical crossings follows from name or centroid.'], 'engineNote': 'Use explicit post_fields after body_envelope_map, with Jpost@Jbody; preserve original body-v2 module/contract hashes and do not silently reuse older regional fields.', 'candidateApplied': False, 'geometryModified': False}
    (out / 'constraint-recipe.json').write_text(json.dumps(recipe, indent=2) + '\n')
    if args.no_render:
        print(json.dumps({'contract': str(out / 'contract.json'), 'recipe': str(out / 'constraint-recipe.json'), 'rendered': False, 'elapsedSeconds': time.monotonic() - began}))
        return
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    full = np.load(base / 'regional-interface-v1/fixed-context/full.npz')
    owners = full['face_owner_part_index']
    chosen = np.asarray([p['partIndex'] for p in r['fixedBoneParts']])
    bone_ids = np.flatnonzero(np.isin(owners, chosen))
    bt = full['vertices'][full['faces'][bone_ids]]
    local = (bt[:, :, 1].max(axis=1) >= .38) & (bt[:, :, 1].min(axis=1) <= .59)
    bt = bt[local]
    sample = np.load(out / 'vessel-sample-ownership.npz')
    vt, source_ids = sample['triangles_m'], sample['source_instance_id']
    knee_vessel = np.asarray([any(k in x.lower() for k in ['popliteal', 'genicular', 'femoral_artery', 'femoral_vein', 'tibial_artery', 'tibial_veins']) for x in source_ids])
    guides = np.load(out / 'popliteal-section-guides.npz')
    old, trial = guides['actual_source_section_points_m'], guides['trial_mapped_section_points_m']
    fig, axes = plt.subplots(1, 2, figsize=(11, 8), sharey=True)
    for ax, side, sign in zip(axes, ['Right', 'Left'], [-1, 1]):
        bone = bt[(bt[:, :, 0].mean(axis=1) * sign) > 0]
        ax.add_collection(PolyCollection(bone[:, :, [2, 1]] * 1000, facecolors='#a5aab0', edgecolors='#6c7379', linewidths=.22, alpha=.22, rasterized=True))
        side_mask = (vt[:, :, 0].mean(axis=1) * sign > 0) & knee_vessel
        for mesh_name, colour in [('arterial', '#bd3430'), ('venous', '#2768b0')]:
            t = vt[side_mask & (sample['candidate_mesh'] == mesh_name)]
            t = t[(t[:, :, 1].max(axis=1) >= .38) & (t[:, :, 1].min(axis=1) <= .59)]
            ax.add_collection(PolyCollection(t[:, :, [2, 1]] * 1000, facecolors=colour, edgecolors=colour, linewidths=.25, alpha=.32, rasterized=True))
        use = old[:, 0] * sign > 0
        ax.scatter(old[use, 2] * 1000, old[use, 1] * 1000, s=5, color='#292929', label='Measured popliteal surface sections', zorder=4)
        ax.scatter(trial[use, 2] * 1000, trial[use, 1] * 1000, s=7, marker='+', color='#008765', label='Proposed common-field section points', zorder=5)
        ax.set_xlim(-115, 50)
        ax.set_ylim(380, 590)
        ax.set_aspect('equal')
        ax.grid(alpha=.16)
        ax.set_title(side + ' knee: sagittal projection')
        ax.set_xlabel('Z (mm): posterior ← / → anterior')
    axes[0].set_ylabel('Y superior (mm)')
    axes[1].legend(loc='upper right', fontsize=7)
    fig.suptitle('Actual fixed bone and saved vascular triangles; common posterior knee trial', fontsize=12)
    fig.text(.5, .015, 'Orthographic triangle projections, not sections. Gray: fixed bone/cartilage; red: arterial; blue: venous. No anatomy changed.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .035, 1, .96))
    fig.savefig(out / 'knee-sagittal-projections.png', dpi=160)
    plt.close(fig)
    print(json.dumps({'recipe': str(out / 'constraint-recipe.json'), 'projection': str(out / 'knee-sagittal-projections.png'), 'elapsedSeconds': time.monotonic() - began}))


if __name__ == '__main__':
    main()
