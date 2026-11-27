"""One native C2Course call per declared whole-stem interval; no wall emission."""
from pathlib import Path
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.tube_c2_course import C2Course

BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-intrapulmonary-circulation-current-v1'
OUTPUT = BASE / 'whole-unbranched-root-stem-c2-v1'
CATALOGUE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def main():
    began = time.monotonic()
    contract_path = OUTPUT / 'before-construction-contract.json'
    assert digest(contract_path) == '5fa919fa6e9fd3798736f7d0d8d06cc65e04e671a8c760dd89fb1c089b21dea3'
    contract = json.loads(contract_path.read_text())
    assert contract['both_intervals_algebraically_compatible']
    assert digest(CATALOGUE) == contract['source_catalogue_sha256']
    helper_path = ROOT / contract['helper']
    assert digest(helper_path) == contract['helper_sha256']
    fixed_parent = BASE / 'collector-left-common-basal-v1/fixed-current-lipv.npz'
    assert digest(fixed_parent) == contract['fixed_parent_sha256']
    result_path = OUTPUT / 'guide-result.json'
    assert not result_path.exists(), 'One solve stage only; no implicit restart'
    execution = {
        'contract_sha256': digest(contract_path),
        'script_sha256': digest(__file__),
        'helper_sha256': digest(helper_path),
        'native_calls': 'At most one per interval,131 then60; stop at first failure.',
        'time_limit_seconds_per_call': 9,
        'total_solve_budget_seconds': 20,
        'no_parameter_restart': True,
    }
    write(OUTPUT / 'before-solve-binding.json', execution)
    result = {
        'schema': 'whole-unbranched-root-stem-C2-guide-result/v1',
        'contract_sha256': digest(contract_path), 'calls': [],
        'passed': False, 'guide_only': True, 'wall_emitted': False,
    }
    for spec in contract['specs']:
        assert spec['algebraic_preflight_passed']
        curve = spec['curve']
        course = C2Course(
            np.asarray(spec['reference_controls_mm']),
            spec['physical_parameter_length_mm'], spec['controls'],
            spec['minimum_bend_radius_mm'], speed_lower=spec['speed_lower'],
            trust_radius_mm=spec['trust_radius_mm'],
            boundary_condition='not-a-knot',
        )
        assert np.array_equal(course.directions, np.asarray(spec['reference_directions']))
        length = course.length
        tangents = spec['endpoint_unit_tangents']
        accelerations = spec['endpoint_second_derivative_jets_per_mm']
        equalities = [
            (0., 1, tangents[0], f'curve{curve}:literal-start-profile-frame'),
            (0., 2, accelerations[0], f'curve{curve}:incoming-estimated-profile-frame-curvature'),
            (length, 1, tangents[1], f'curve{curve}:literal-held-distal-profile-frame'),
            (length, 2, accelerations[1], f'curve{curve}:outgoing-estimated-profile-frame-curvature'),
        ]
        constraints = [
            (float(parameter), contract['fixed_basal_plane_distal_normal'], 0.,
             {'kind': 'fixed-basal-plane-centre', 'curve': curve, 'parameter_mm': float(parameter)})
            for parameter in np.linspace(0, length, 64)
        ]
        remaining = 20. - sum(call['elapsed_seconds'] for call in result['calls'])
        if remaining <= 0:
            result.update(status='TotalSolveBudgetExhausted', stopped_curve=curve)
            break
        controls, diagnostic = course.solve(
            constraints, time_limit=min(9., remaining), equalities=equalities,
        )
        result['calls'].append({'curve': curve, **diagnostic})
        if controls is None:
            result.update(status=diagnostic['status'], stopped_curve=curve,
                          qualification='One declared conic guide problem; no restart or global anatomical impossibility claim.')
            break
        path = OUTPUT / f'curve-{curve}-guide.npz'
        np.savez_compressed(path, controls_mm=controls, knots_mm=course.knots,
                            reference_mm=course.reference)
        result['calls'][-1]['guide_path'] = str(path)
        result['calls'][-1]['guide_sha256'] = digest(path)
        if not diagnostic['certificate']['continuous_curvature_bound_pass']:
            result.update(status='ContinuousGuideCertificateFailed', stopped_curve=curve)
            break
    else:
        result.update(status='BothContinuousGuideCertificatesPassed', passed=True)

    result['actual_profile_sweep_and_join_gates'] = 'Pending separate bounded stage; no wall or CSG authorized.'
    result['elapsed_seconds'] = time.monotonic() - began
    assert digest(CATALOGUE) == contract['source_catalogue_sha256']
    assert digest(fixed_parent) == contract['fixed_parent_sha256']
    write(result_path, result)
    print(json.dumps({'path': str(result_path), 'sha256': digest(result_path),
                      'status': result['status'], 'solve_calls': len(result['calls']),
                      'elapsed_seconds': result['elapsed_seconds']}), flush=True)


if __name__ == '__main__':
    main()
