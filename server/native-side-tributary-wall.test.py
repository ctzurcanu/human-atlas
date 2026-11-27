"""End-to-end nonanatomical three-port finite-wall regression.

Run with the geometry environment through geometry_tools.runtime --seconds 60.
The fixture and its receipts are retained locally for failed-case inspection.
No atlas, authority database, anatomical source or installation is modified.
"""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile
import time
import unittest

import igl
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.tube_continuity import boundary_loops  # noqa: E402
from geometry_tools.common import atomic_json, atomic_npz, digest_file  # noqa: E402
from geometry_tools.model import Mesh  # noqa: E402

PRODUCER = ROOT/'server/construct-native-side-tributary-wall.py'
FRAME = 'nonanatomical-test:orthogonal-three-port:metres'


def fixture_tube(path, centers, first_axis, second_axis, radius):
    """A regular polygonal open pipe; no anatomical geometry is consumed."""
    angles = np.arange(12)*2*np.pi/12 + .137
    profile = radius*(np.cos(angles)[:, None]*np.asarray(first_axis)
                      + np.sin(angles)[:, None]*np.asarray(second_axis))
    profiles = np.asarray(centers)[:, None] + profile
    indices = np.arange(profiles.shape[0]*12).reshape(-1, 12)
    faces = []
    for lower, upper in zip(indices[:-1], indices[1:]):
        for j in range(12):
            following = (j+1) % 12
            faces.extend([[lower[j], lower[following], upper[following]],
                          [lower[j], upper[following], upper[j]]])
    atomic_npz(path, vertices=profiles.reshape(-1, 3), faces=np.asarray(faces, dtype=np.int64),
               profiles=profiles, profile_vertex_indices=indices, frame=np.array(FRAME))
    return profiles.reshape(-1, 3)


class NativeSideTributaryWallTest(unittest.TestCase):
    def test_actual_cli_three_port_wall_and_axis_guard(self):
        began = time.monotonic()
        evidence_root = ROOT/'.local-models/geometry-workbench/testing'
        evidence_root.mkdir(exist_ok=True, parents=True)
        directory = Path(tempfile.mkdtemp(prefix='native-side-tributary-', dir=evidence_root))
        parent_path, branch_path = directory/'parent.npz', directory/'tributary.npz'
        parent_centers = np.zeros((13, 3))
        parent_centers[:, 2] = np.linspace(-.003, .003, 13)
        branch_centers = np.zeros((7, 3))
        branch_centers[:, 0] = np.linspace(0, .003, 7)
        parent = fixture_tube(parent_path, parent_centers, [1, 0, 0], [0, 1, 0], .0013)
        branch = fixture_tube(branch_path, branch_centers, [0, 1, 0], [0, 0, 1], .001)
        inventory = directory/'fixture-inventory.json'
        guard = directory/'fixture-guard.json'
        atomic_json(inventory, dict(schema='nonanatomical-test-inventory/v1',
            anatomical_terms=[], actual_objects=['synthetic parent pipe', 'synthetic side pipe'],
            source_ownership='Original analytical test data', anatomical_acceptance=False))
        atomic_json(guard, dict(nonanatomical=True, no_atlas_sources=True, no_installation=True))
        binding = lambda path: dict(path=str(path), sha256=digest_file(path))
        plan = dict(schema='native-side-tributary-finite-wall/v1',
            plan_validated_for_construction=True, anatomical_acceptance=False, frame=FRAME,
            parent=binding(parent_path), tributary=binding(branch_path),
            inventory=binding(inventory), guards=[binding(guard)], mode='native-overlap',
            parent_anchor=dict(station=6, fraction=0.), branch_mouth=0,
            max_axis_residual_m=1e-7, inner_radial_fraction=.88,
            junction_support_radius_m=.0025, storage_origin_m=[0., 0., 0.])
        plan_path = directory/'plan.json'
        atomic_json(plan_path, plan)
        output = directory/'positive'
        runs = []

        def run(stage, path=plan_path, destination=output, accepted=True, label=None):
            start = time.monotonic()
            # The outer unittest also has a 60-second supervisor. Each child
            # has its own finite cap and separate elapsed-time receipt.
            command = [sys.executable, '-B', '-m', 'geometry_tools.runtime', '--seconds', '10', '--',
                       sys.executable, '-B', str(PRODUCER), '--plan', str(path),
                       '--plan-sha256', digest_file(path), '--output', str(destination), '--stage', stage]
            process = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=12)
            row = dict(stage=stage, label=label or stage, returncode=process.returncode,
                       elapsed_seconds=time.monotonic()-start, command=command,
                       stdout=process.stdout, stderr=process.stderr)
            runs.append(row)
            atomic_json(directory/f'run-{row["label"]}.json', row)
            if accepted:
                self.assertEqual(process.returncode, 0, process.stdout + process.stderr)
            else:
                self.assertNotEqual(process.returncode, 0, 'Unsafe negative case unexpectedly succeeded')

        for stage in ('prepare', 'outer', 'inner', 'assemble'):
            run(stage)
        receipt = json.loads((output/'assemble.json').read_text())
        wall_audit = json.loads((output/'wall-audit.json').read_text())
        blood_audit = json.loads((output/'blood-audit.json').read_text())
        self.assertTrue(receipt['passed'])
        self.assertEqual(receipt['external_blood_mouths'], 3)
        self.assertEqual(receipt['closed_finite_wall_genus'], 2)
        self.assertEqual(receipt['closed_blood_classifier_genus'], 0)
        self.assertTrue(wall_audit['passed'])
        self.assertTrue(blood_audit['passed'])
        for checked, genus in ((wall_audit, 2), (blood_audit, 0)):
            self.assertEqual(len(checked['topology']['components']), 1)
            self.assertEqual(checked['topology']['components'][0]['closed_orientable_genus'], genus)
            self.assertTrue(checked['self_intersections']['complete_within_pair_budget'])
            self.assertEqual(checked['self_intersections']['defect_pairs'], 0)
        self.assertEqual(receipt['internal_classifier_caps'], 0)
        self.assertEqual(receipt['lumen_disks'], 0)
        self.assertFalse(receipt['anatomical_acceptance'])
        self.assertFalse(receipt['installed'])
        with np.load(output/'prepared.npz', allow_pickle=False) as prepared:
            np.testing.assert_array_equal(prepared['parent_outer'], parent)
            np.testing.assert_array_equal(prepared['tributary_outer'], branch)
        outer = Mesh.load(output/'outer-sheet.npz')
        inner = Mesh.load(output/'inner-sheet.npz')
        self.assertEqual(len(boundary_loops(outer.faces)), 3)
        self.assertEqual(len(boundary_loops(inner.faces)), 3)
        used_outer = {tuple(p) for p in outer.vertices[np.unique(outer.faces)]}
        far_native = np.r_[parent, branch]
        far_native = far_native[np.linalg.norm(far_native, axis=1) > .0025]
        self.assertGreater(len(far_native), 0)
        self.assertTrue(all(tuple(p) in used_outer for p in far_native),
                        'A protected remote native control moved or disappeared')
        ports = json.loads((output/'ports.json').read_text())['ports']
        with np.load(output/'surfaces-and-ancestry.npz', allow_pickle=False) as surfaces:
            for i, port in enumerate(ports):
                rim = surfaces[f'annulus_{i}']
                outside_ids, inside_ids = set(port['outer_indices']), set(port['inner_indices'])
                self.assertEqual(set(np.unique(rim)), outside_ids | inside_ids)
                self.assertTrue(all(set(face) & outside_ids and set(face) & inside_ids for face in rim),
                                'An external closure includes a lumen-spanning disk')
        # Independently classify actual emitted geometry along all three arms.
        # Winding is one in the blood domain and zero in tissue on the axis.
        z = np.linspace(-.00295, .00295, 61)
        parent_axis = np.column_stack([np.zeros_like(z), np.zeros_like(z), z])
        x = np.linspace(0., .00295, 31)
        branch_axis = np.column_stack([x, np.zeros_like(x), np.zeros_like(x)])
        queries = np.r_[parent_axis, branch_axis]
        blood = Mesh.load(output/'inner-classifier.npz')
        wall = Mesh.load(output/'finite-wall.npz')
        blood_winding = np.asarray(igl.winding_number(blood.vertices, blood.faces, queries))
        tissue_winding = np.asarray(igl.winding_number(wall.vertices, wall.faces, queries))
        np.testing.assert_allclose(blood_winding, 1., atol=1e-8, rtol=0.)
        np.testing.assert_allclose(tissue_winding, 0., atol=1e-8, rtol=0.)
        # Deliberately bind the same source mouth to the remote parent end.
        # The guard must fail before even a prepared candidate is emitted.
        negative_plan = dict(plan, parent_anchor=dict(station=0, fraction=0.))
        negative_path = directory/'negative-axis-plan.json'
        atomic_json(negative_path, negative_plan)
        negative_output = directory/'negative-axis'
        run('prepare', negative_path, negative_output, accepted=False, label='negative-axis')
        failure = json.loads((negative_output/'prepare-failure.json').read_text())
        self.assertIn('not at the declared parent axis anchor', failure['error'])
        self.assertFalse((negative_output/'prepared.npz').exists())
        self.assertFalse((negative_output/'finite-wall.npz').exists())
        evidence = dict(passed=True, nonanatomical_fixture=True, no_atlas_changes=True,
            producer=str(PRODUCER), producer_sha256=digest_file(PRODUCER),
            test_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            positive_receipt=str(output/'assemble.json'),
            positive_receipt_sha256=digest_file(output/'assemble.json'),
            negative_receipt=str(negative_output/'prepare-failure.json'),
            negative_receipt_sha256=digest_file(negative_output/'prepare-failure.json'),
            three_arm_axis_samples=len(queries), blood_axis_occupied=True, tissue_axis_empty=True,
            minimum_blood_winding=float(blood_winding.min()),
            maximum_absolute_tissue_winding=float(abs(tissue_winding).max()),
            native_controls_fixed=True, annular_external_rims_only=True,
            child_runs=[dict(stage=r['stage'], label=r['label'], returncode=r['returncode'],
                             elapsed_seconds=r['elapsed_seconds']) for r in runs],
            elapsed_seconds=time.monotonic()-began)
        atomic_json(directory/'test-evidence.json', evidence)
        print(json.dumps(dict(test_evidence=str(directory/'test-evidence.json'),
                              sha256=digest_file(directory/'test-evidence.json'),
                              elapsed_seconds=evidence['elapsed_seconds'])))


if __name__ == '__main__':
    unittest.main()
