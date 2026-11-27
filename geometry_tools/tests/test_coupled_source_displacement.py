"""Explicit source-mouth identity and original displacement-field invariants."""
import unittest
from unittest.mock import patch

import numpy as np

from geometry_tools.algorithms.coupled_source_displacement import solve_coupled_source_displacement
from geometry_tools.model import Mesh


def tube(z=0., *, alias=False):
    sides, rings = 6, 4
    angles = np.arange(sides) * 2 * np.pi / sides
    points = np.array([[.001 * np.cos(a), .001 * np.sin(a), z + ring * .001]
                       for ring in range(rings) for a in angles])
    faces = []
    for ring in range(rings - 1):
        for side in range(sides):
            a, b = ring * sides + side, ring * sides + (side + 1) % sides
            faces.extend([[a, b, b + sides], [a, b + sides, a + sides]])
    faces = np.array(faces)
    if alias:
        points = np.vstack([points, points[18]])
        location = np.argwhere(faces == 18)[0]
        faces[tuple(location)] = 24
    return Mesh(points, faces, 'test:metres', {'uv': np.arange(len(points) * 2).reshape(-1, 2)})


def pair(gap=0., *, alias=False):
    meshes = {'a': tube(alias=alias), 'b': tube(.003 + gap)}
    interface = dict(first='a', second='b', firstRing=list(range(18, 24)), secondRing=list(range(6)))
    return meshes, [interface]


def declaration(mesh, raw, shift):
    return dict(rawIds=list(raw), targetPoints_m=mesh.vertices[raw] + shift)


def solve(meshes, interfaces, *, prescribed=None, fixed=None, tolerance=1e-15, **kwargs):
    return solve_coupled_source_displacement(
        meshes, interfaces, prescribed=prescribed or {}, fixed=fixed or {},
        declared_target_tolerance_m=tolerance, pair_mode='preserve-source-offsets', **kwargs)


class CoupledSourceDisplacementTests(unittest.TestCase):
    def test_zero_motion_retains_all_original_words_including_declared_mouth_offsets(self):
        meshes, interfaces = pair(gap=2e-7, alias=True)
        originals = {pid: (mesh.vertices.copy(), mesh.faces.copy(), mesh.attributes['uv'].copy())
                     for pid, mesh in meshes.items()}
        q, proof = solve(meshes, interfaces)
        for pid, mesh in meshes.items():
            np.testing.assert_array_equal(q[pid], originals[pid][0])
            np.testing.assert_array_equal(mesh.vertices, originals[pid][0])
            np.testing.assert_array_equal(mesh.faces, originals[pid][1])
            np.testing.assert_array_equal(mesh.attributes['uv'], originals[pid][2])
        self.assertEqual(proof['pairMode'], 'preserve-source-offsets')
        self.assertTrue(proof['interfaces'][0]['sourceProof']['passed'])
        self.assertFalse(proof['anatomicalAcceptance'])
        self.assertFalse(proof['triangleQualityAndCalibreVerified'])
        self.assertFalse(proof['foreignContactsVerified'])

    def test_one_prescribed_shared_mouth_translates_both_complete_sources(self):
        meshes, interfaces = pair(gap=2e-7, alias=True)
        shift = np.array([.0001, -.00004, .00003])
        prescribed = {'a': declaration(meshes['a'], np.arange(18, 24), shift)}
        q, proof = solve(meshes, interfaces, prescribed=prescribed)
        for pid, mesh in meshes.items():
            np.testing.assert_allclose(q[pid], mesh.vertices + shift, atol=2e-18, rtol=0)
        np.testing.assert_array_equal(q['a'][18], q['a'][24])
        old_offset = meshes['b'].vertices[:6] - meshes['a'].vertices[18:24]
        np.testing.assert_allclose(q['b'][:6] - q['a'][18:24], old_offset, atol=2e-18, rtol=0)
        self.assertLess(proof['equationMaximumResidual'], 1e-12)
        self.assertEqual(proof['prescribedRawControlCount'], 6)

    def test_fixed_far_boundaries_and_common_mouth_pose_are_preserved(self):
        meshes, interfaces = pair(alias=True)
        shift = np.array([.0001, 0., 0.])
        prescribed = {'a': declaration(meshes['a'], np.arange(18, 24), shift),
                      'b': declaration(meshes['b'], np.arange(6), shift)}
        q, proof = solve(meshes, interfaces, prescribed=prescribed,
                         fixed={'a': list(range(6)), 'b': list(range(18, 24))})
        np.testing.assert_array_equal(q['a'][:6], meshes['a'].vertices[:6])
        np.testing.assert_array_equal(q['b'][18:24], meshes['b'].vertices[18:24])
        np.testing.assert_array_equal(q['a'][18:24], prescribed['a']['targetPoints_m'])
        np.testing.assert_array_equal(q['b'][:6], prescribed['b']['targetPoints_m'])
        for pid, raw in [('a', np.arange(6, 18)), ('b', np.arange(6, 18))]:
            displacement = q[pid][raw] - meshes[pid].vertices[raw]
            self.assertTrue(np.all((displacement[:, 0] > 0) & (displacement[:, 0] < shift[0])))
            np.testing.assert_array_equal(displacement[:, 1:], np.zeros((12, 2)))
        self.assertTrue(proof['fixedSourceCoordinatesLiteral'])

    def test_coincident_unrelated_owners_are_never_coupled(self):
        meshes = {'a': tube(), 'b': tube()}
        prescribed = {'a': declaration(meshes['a'], np.arange(6), np.array([.0001, 0, 0]))}
        q, proof = solve(meshes, [], prescribed=prescribed)
        np.testing.assert_array_equal(q['b'], meshes['b'].vertices)
        self.assertGreater(np.linalg.norm(q['a'] - meshes['a'].vertices), 0)
        self.assertEqual(proof['coupledPhysicalGroups'], proof['sourcePhysicalControls'])
        self.assertEqual(len(proof['unanchoredSourceComponentsRetainedLiteral']), 1)

    def test_original_inverse_lengths_choose_the_displacement_energy_minimum(self):
        # The only free vertex has edges of lengths 1 and 2 to fixed/prescribed
        # vertices. Its shift is (0/1 + .3/2)/(1/1 + 1/2) = .1, not the .15
        # produced by unit weights or a solve of absolute vertex positions.
        points = np.array([[0., 0., 0.], [2., 0., 0.], [2., 1., 0.], [0., 1., 0.]])
        mesh = Mesh(points, np.array([[0, 1, 2], [0, 2, 3]]), 'test:metres')
        q, _ = solve({'patch': mesh}, [], fixed={'patch': [0, 1]},
                     prescribed={'patch': declaration(mesh, [2], np.array([.3, 0, 0]))})
        np.testing.assert_allclose(q['patch'][3], [.1, 1., 0.], atol=1e-15, rtol=0)

    def test_only_declared_paired_targets_may_snap_within_explicit_tolerance(self):
        meshes, interfaces = pair(gap=2e-7)
        shift = np.array([.0001, 0, 0])
        prescribed = {'a': declaration(meshes['a'], np.arange(18, 24), shift),
                      'b': declaration(meshes['b'], np.arange(6), shift + [2e-9, 0, 0])}
        q, proof = solve(meshes, interfaces, prescribed=prescribed, tolerance=3e-9)
        self.assertEqual(len(proof['snappedDeclaredGroups']), 6)
        self.assertAlmostEqual(proof['maximumPrescribedReplayResidual_m'], 2e-9, places=16)
        np.testing.assert_allclose(q['b'][:6], meshes['b'].vertices[:6] + shift, atol=2e-18, rtol=0)
        self.assertTrue(all(len(group['constraints']) == 2 for group in proof['snappedDeclaredGroups']))
        with self.assertRaisesRegex(ValueError, 'agreement tolerance'):
            solve(meshes, interfaces, prescribed=prescribed, tolerance=1e-9)

    def test_fixed_constraint_is_hard_even_with_a_larger_declared_target_tolerance(self):
        meshes, interfaces = pair(gap=2e-7)
        prescribed = {'b': declaration(meshes['b'], np.arange(6), np.array([1e-10, 0, 0]))}
        with self.assertRaisesRegex(ValueError, 'hard fixed'):
            solve(meshes, interfaces, prescribed=prescribed, fixed={'a': [18]}, tolerance=1e-6)

    def test_exact_alias_targets_cannot_use_the_cross_owner_snap_tolerance(self):
        meshes, interfaces = pair(alias=True)
        targets = meshes['a'].vertices[[18, 24]].copy()
        targets[1, 0] += 1e-10
        with self.assertRaisesRegex(ValueError, 'Exact aliases'):
            solve(meshes, interfaces, prescribed={'a': dict(rawIds=[18, 24], targetPoints_m=targets)}, tolerance=1e-6)

    def test_incomplete_and_noncorresponding_mouths_are_rejected_before_coupling(self):
        meshes, interfaces = pair()
        incomplete = dict(interfaces[0], firstRing=list(range(18, 23)), secondRing=list(range(5)))
        with self.assertRaisesRegex(ValueError, 'complete external source mouth'):
            solve(meshes, [incomplete])
        wrong = dict(interfaces[0], secondRing=list(range(1, 6)) + [0])
        with self.assertRaisesRegex(ValueError, 'source gate'):
            solve(meshes, [wrong])
        unequal = dict(interfaces[0], secondRing=list(range(5)))
        with self.assertRaisesRegex(ValueError, 'same-arity'):
            solve(meshes, [unequal])

    def test_wall_deadline_does_not_emit_a_partial_field(self):
        meshes, interfaces = pair()
        with patch('geometry_tools.algorithms.coupled_source_displacement.time.monotonic', side_effect=[0., 10.]):
            with self.assertRaises(TimeoutError):
                solve(meshes, interfaces)

    def test_nonfinite_targets_invalid_ids_and_work_limits_are_rejected(self):
        meshes, interfaces = pair()
        for seconds in [0., 5.01, float('nan')]:
            with self.subTest(seconds=seconds), self.assertRaises(ValueError):
                solve(meshes, interfaces, wall_seconds=seconds)
        for tolerance in [-1., float('inf'), float('nan')]:
            with self.subTest(tolerance=tolerance), self.assertRaises(ValueError):
                solve(meshes, interfaces, tolerance=tolerance)
        bad = dict(rawIds=[0], targetPoints_m=[[float('nan'), 0, 0]])
        with self.assertRaises(ValueError):
            solve(meshes, interfaces, prescribed={'a': bad})
        with self.assertRaises(ValueError):
            solve(meshes, interfaces, fixed={'a': [100]})
        large = Mesh(np.zeros((150001, 3)), np.array([[0, 1, 2]]), 'test:metres')
        with self.assertRaisesRegex(ValueError, '150,000'):
            solve({'large': large}, [])


if __name__ == '__main__':
    unittest.main()
