"""Actual section shape, fixed-star closure and explicit counterpart ownership."""
import unittest
from unittest.mock import patch

import numpy as np

from geometry_tools.algorithms.project_coupled_source_sections import propose_section_targets
from geometry_tools.model import Mesh
from geometry_tools.tests.test_coupled_source_displacement import tube, pair


def pose(points, scale=1., shift=(.0003, -.0002, .0001)):
    angle = .23
    rotation = np.array([[np.cos(angle), 0., np.sin(angle)], [0., 1., 0.],
                         [-np.sin(angle), 0., np.cos(angle)]])
    # Compute once per actual source alias. Different SIMD batches may otherwise
    # round duplicate raw rows differently, which is not a consistent field.
    physical, inverse = np.unique(points, axis=0, return_inverse=True)
    return (scale * physical @ rotation.T + shift)[inverse], rotation


def project(meshes, profiles, proposed, *, interfaces=(), fixed=None, maximum=1.030452, minimum=1., **kwargs):
    return propose_section_targets(meshes, list(interfaces), profiles, proposed,
                                    fixed or {}, maximum_uniform_scale=maximum,
                                    minimum_uniform_scale=minimum, **kwargs)


def targets(mesh, declaration):
    q = mesh.vertices.copy()
    q[declaration['rawIds']] = declaration['targetPoints_m']
    return q


def pair_distances(points):
    return np.linalg.norm(points[:, None] - points[None], axis=2)


class ProjectCoupledSourceSectionsTests(unittest.TestCase):
    def test_known_proper_scaled_pose_retains_full_section_shape_and_aliases(self):
        mesh = tube(alias=True)
        q, _ = pose(mesh.vertices, scale=1.02)
        raw = np.arange(18, 24)
        prescribed, fixed, proof = project({'a': mesh}, {'a': [raw]}, {'a': q})
        result = targets(mesh, prescribed['a'])
        np.testing.assert_allclose(result[raw], q[raw], atol=2e-18, rtol=0)
        np.testing.assert_allclose(pair_distances(result[raw]), 1.02 * pair_distances(mesh.vertices[raw]), atol=2e-18, rtol=0)
        np.testing.assert_array_equal(result[18], result[24])
        self.assertEqual(set(prescribed['a']['rawIds']), set(raw) | {24})
        self.assertEqual(fixed['a'], [])
        self.assertEqual(proof['rigidCells'], 1)
        self.assertFalse(proof['anatomicalAcceptance'])
        self.assertFalse(proof['unlistedBranchCoresIncluded'])

    def test_uniform_scale_is_bounded_below_and_above_without_shape_distortion(self):
        mesh = tube()
        raw = np.arange(6, 12)
        for source_scale, wanted in [(.8, 1.), (1.2, 1.030452)]:
            with self.subTest(scale=source_scale):
                q, _ = pose(mesh.vertices, scale=source_scale)
                prescribed, _, proof = project({'a': mesh}, {'a': [raw]}, {'a': q})
                result = targets(mesh, prescribed['a'])
                np.testing.assert_allclose(pair_distances(result[raw]), wanted * pair_distances(mesh.vertices[raw]), atol=2e-18, rtol=0)
                np.testing.assert_allclose(result[raw].mean(0), q[raw].mean(0), atol=2e-18, rtol=0)
                self.assertEqual(proof['cells'][0]['scale'], wanted)

    def test_explicit_shrinking_range_retains_a_known_point_94_similarity(self):
        mesh = tube(alias=True)
        raw = np.arange(18, 24)
        q, _ = pose(mesh.vertices, scale=.94)
        prescribed, _, proof = project({'a': mesh}, {'a': [raw]}, {'a': q}, minimum=.9, maximum=1.)
        result = targets(mesh, prescribed['a'])
        np.testing.assert_allclose(result[raw], q[raw], atol=2e-18, rtol=0)
        np.testing.assert_allclose(pair_distances(result[raw]), .94 * pair_distances(mesh.vertices[raw]), atol=2e-18, rtol=0)
        self.assertAlmostEqual(proof['cells'][0]['scale'], .94)
        self.assertEqual(proof['minimumUniformScale'], .9)
        np.testing.assert_array_equal(result[18], result[24])

    def test_fixed_cells_keep_unit_scale_outside_the_free_cell_shrinking_range(self):
        mesh = tube()
        raw = np.arange(18, 24)
        q, _ = pose(mesh.vertices, scale=.94)
        prescribed, fixed, proof = project({'a': mesh}, {'a': [raw]}, {'a': q},
                                           fixed={'a': [18]}, minimum=.9, maximum=.95)
        np.testing.assert_array_equal(targets(mesh, prescribed['a'])[raw], mesh.vertices[raw])
        self.assertEqual(fixed['a'], raw.tolist())
        self.assertEqual(proof['cells'][0]['scale'], 1.)
        self.assertTrue(proof['fixedCellScaleSourcePriority'])

    def test_invalid_uniform_scale_intervals_are_rejected(self):
        mesh = tube()
        for minimum, maximum in [(.49, 1.), (1.1, 1.03), (.9, 2.01),
                                  (float('nan'), 1.), (.9, float('inf'))]:
            with self.subTest(minimum=minimum, maximum=maximum), self.assertRaisesRegex(ValueError, 'scale range'):
                project({'a': mesh}, {'a': [np.arange(6)]}, {'a': mesh.vertices.copy()},
                        minimum=minimum, maximum=maximum)

    def test_fixed_alias_star_keeps_its_whole_cell_literal_but_does_not_recruit_core(self):
        mesh = tube(alias=True)
        q, _ = pose(mesh.vertices, scale=1.02)
        rings = [np.arange(18, 24), np.arange(6, 12)]
        prescribed, fixed, proof = project({'a': mesh}, {'a': rings}, {'a': q}, fixed={'a': [24, 1]})
        result = targets(mesh, prescribed['a'])
        np.testing.assert_array_equal(result[18:24], mesh.vertices[18:24])
        np.testing.assert_array_equal(result[24], mesh.vertices[24])
        self.assertEqual(set(fixed['a']), {1, 18, 19, 20, 21, 22, 23, 24})
        np.testing.assert_allclose(result[6:12], q[6:12], atol=2e-18, rtol=0)
        self.assertNotIn(12, prescribed['a']['rawIds'])
        self.assertEqual(proof['fixedCells'], 1)

    def test_overlapping_listed_source_cycles_form_one_fixed_cell(self):
        mesh = Mesh(np.array([[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.], [0., 0., .001]]),
                    np.array([[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]]), 'test:metres')
        q, _ = pose(mesh.vertices)
        prescribed, fixed, proof = project({'a': mesh}, {'a': [[0, 1, 2], [0, 2, 3]]}, {'a': q}, fixed={'a': [1]})
        np.testing.assert_array_equal(targets(mesh, prescribed['a']), mesh.vertices)
        self.assertEqual(fixed['a'], [0, 1, 2, 3])
        self.assertEqual(proof['rigidCells'], 1)

    def test_explicit_paired_mouths_have_one_pose_and_keep_original_chart_offsets(self):
        meshes, interfaces = pair(gap=2e-7, alias=True)
        proposed = {pid: pose(mesh.vertices, scale=1.02)[0] for pid, mesh in meshes.items()}
        prescribed, _, proof = project(meshes, {'a': [np.arange(18, 24)]}, proposed, interfaces=interfaces)
        first, second = targets(meshes['a'], prescribed['a']), targets(meshes['b'], prescribed['b'])
        np.testing.assert_allclose(first[18:24], proposed['a'][18:24], atol=3e-18, rtol=0)
        np.testing.assert_allclose(second[:6], proposed['b'][:6], atol=3e-18, rtol=0)
        _, rotation = pose(meshes['a'].vertices)
        old_offset = meshes['b'].vertices[:6] - meshes['a'].vertices[18:24]
        np.testing.assert_allclose(second[:6] - first[18:24], 1.02 * old_offset @ rotation.T, atol=3e-18, rtol=0)
        self.assertEqual(proof['rigidCells'], 1)
        self.assertEqual(set(prescribed['b']['rawIds']), set(range(6)))
        self.assertGreater(proof['unlistedCorePhysicalControls'], 0)

    def test_unrelated_coincident_owners_keep_independent_section_cells(self):
        meshes = {'a': tube(), 'b': tube()}
        proposed = {'a': meshes['a'].vertices + [.0001, 0, 0],
                    'b': meshes['b'].vertices + [0, -.0002, 0]}
        profiles = {pid: [np.arange(6)] for pid in meshes}
        prescribed, _, proof = project(meshes, profiles, proposed)
        for pid in meshes:
            np.testing.assert_allclose(np.array(prescribed[pid]['targetPoints_m']), proposed[pid][:6], atol=1e-18, rtol=0)
        self.assertEqual(proof['rigidCells'], 2)

    def test_reflected_nonplanar_target_cannot_make_an_improper_section_pose(self):
        mesh = Mesh(np.array([[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.], [0., 0., .001]]),
                    np.array([[0, 1, 2], [0, 2, 3], [0, 3, 1], [1, 3, 2]]), 'test:metres')
        reflected = mesh.vertices * [-1, 1, 1]
        prescribed, _, proof = project({'a': mesh}, {'a': [[0, 1, 2], [0, 2, 3]]}, {'a': reflected})
        result = targets(mesh, prescribed['a'])
        self.assertGreater(proof['cells'][0]['determinant'], .999999999999)
        np.testing.assert_allclose(pair_distances(result), pair_distances(mesh.vertices), atol=2e-18, rtol=0)
        self.assertGreater(np.linalg.norm(result - reflected), .0001)

    def test_invalid_cycle_incomplete_mouth_and_alias_conflicts_are_rejected(self):
        meshes, interfaces = pair(alias=True)
        proposed = {pid: mesh.vertices.copy() for pid, mesh in meshes.items()}
        with self.assertRaisesRegex(ValueError, 'actual original source edges'):
            project(meshes, {'a': [[0, 2, 4]]}, proposed)
        bad = dict(interfaces[0], firstRing=list(range(18, 23)), secondRing=list(range(5)))
        with self.assertRaisesRegex(ValueError, 'complete external source mouth'):
            project(meshes, {}, proposed, interfaces=[bad])
        proposed['a'][24, 0] += 1e-12
        with self.assertRaisesRegex(ValueError, 'aliases'):
            project(meshes, {'a': [np.arange(18, 24)]}, proposed)

    def test_empty_section_scope_and_native_deadline_are_explicit(self):
        mesh = tube()
        prescribed, fixed, proof = project({'a': mesh}, {}, {'a': mesh.vertices.copy()}, fixed={'a': [2]})
        self.assertEqual(prescribed, {})
        self.assertEqual(fixed, {'a': [2]})
        self.assertEqual(proof['rigidCells'], 0)
        with patch('geometry_tools.algorithms.project_coupled_source_sections.time.monotonic', side_effect=[0., 10.]):
            with self.assertRaises(TimeoutError):
                project({'a': mesh}, {'a': [np.arange(6)]}, {'a': mesh.vertices.copy()})


if __name__ == '__main__':
    unittest.main()
