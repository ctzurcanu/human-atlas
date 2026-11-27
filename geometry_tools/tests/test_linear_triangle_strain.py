"""Closed-form field stretches against independent deformed-triangle SVDs."""
import unittest

import numpy as np

from geometry_tools.algorithms.linear_triangle_strain import LinearTriangleStrain


def independent_stretches(points, faces, displacement, weight):
    """Refactor every actual source and deformed triangle for this weight."""
    values = []
    deformed = points + weight * displacement
    for face in faces:
        source = (points[face[1:]] - points[face[0]]).T
        target = (deformed[face[1:]] - deformed[face[0]]).T
        _, singular, vh = np.linalg.svd(source, full_matrices=False)
        # Map a source orthonormal tangent basis into the new 3-D triangle.
        tangent_map = target @ vh.T @ np.diag(1 / singular)
        values.append(np.sort(np.linalg.svd(tangent_map, compute_uv=False)))
    return np.array(values)


class LinearTriangleStrainTests(unittest.TestCase):
    def setUp(self):
        self.points = np.array([
            [-.3, .2, .1], [1.4, -.1, .4], [.2, 1.1, -.4],
            [.5, .4, 1.8], [-.7, .9, .6], [1.2, .8, 1.1],
        ])
        self.faces = np.array([[0, 1, 2], [0, 3, 4], [1, 3, 5], [2, 4, 5]])

    def test_nontrivial_linear_fields_agree_with_fresh_3_by_2_svds(self):
        shear = self.points @ np.array([[.1, .2, 0], [-.15, 0, .08], [.02, -.12, .07]])
        nonlinear = np.column_stack([
            .08 * np.sin(self.points[:, 1]),
            .06 * self.points[:, 0] * self.points[:, 2],
            -.09 * np.cos(self.points[:, 0]),
        ])
        random = np.random.default_rng(183).normal(0, .08, self.points.shape)
        for displacement in [shear, nonlinear, random]:
            metric = LinearTriangleStrain(self.points, self.faces, displacement)
            for weight in [-.7, 0., .2, 1., 2.]:
                with self.subTest(field=displacement.tolist(), weight=weight):
                    expected = independent_stretches(self.points, self.faces, displacement, weight)
                    actual = metric.stretches(weight)
                    np.testing.assert_allclose(actual, expected, atol=1e-13, rtol=1e-13)
                    self.assertAlmostEqual(metric.maximum_strain(weight), float(abs(expected - 1).max()), places=13)

    def test_uniform_scaling_has_the_same_two_stretches_for_every_face(self):
        metric = LinearTriangleStrain(self.points, self.faces, .2 * self.points)
        for weight in [-2., 0., .5, 2.]:
            with self.subTest(weight=weight):
                expected = abs(1 + .2 * weight)
                np.testing.assert_allclose(metric.stretches(weight), np.full((4, 2), expected), atol=1e-13)
                self.assertAlmostEqual(metric.maximum_strain(weight), abs(expected - 1), places=13)

    def test_uniform_translation_has_no_triangle_strain(self):
        shift = np.broadcast_to([.2, -.4, 1.7], self.points.shape).copy()
        metric = LinearTriangleStrain(self.points, self.faces, shift)
        for weight in [-3., 0., .5, 4.]:
            np.testing.assert_array_equal(metric.stretches(weight), np.ones((4, 2)))
            self.assertEqual(metric.maximum_strain(weight), 0.)

    def test_degenerate_source_triangles_are_rejected(self):
        face = np.array([[0, 1, 2]])
        for points in [np.zeros((3, 3)), np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]])]:
            with self.subTest(points=points.tolist()), self.assertRaisesRegex(ValueError, 'degenerate'):
                LinearTriangleStrain(points, face, np.zeros_like(points))

    def test_nonfinite_source_fields_and_weights_are_rejected(self):
        zeros = np.zeros_like(self.points)
        for bad in [float('nan'), float('inf'), -float('inf')]:
            for change_source in [True, False]:
                points, field = self.points.copy(), zeros.copy()
                (points if change_source else field)[1, 2] = bad
                with self.subTest(value=bad, source=change_source), self.assertRaises(ValueError):
                    LinearTriangleStrain(points, self.faces, field)
            metric = LinearTriangleStrain(self.points, self.faces, zeros)
            with self.subTest(weight=bad), self.assertRaises(ValueError):
                metric.stretches(bad)

    def test_source_and_field_inputs_remain_literal(self):
        points, faces = self.points.copy(), self.faces.copy()
        field = np.random.default_rng(7).normal(0, .01, points.shape)
        originals = [item.copy() for item in [points, faces, field]]
        metric = LinearTriangleStrain(points, faces, field)
        metric.maximum_strain(.7)
        for actual, original in zip([points, faces, field], originals):
            np.testing.assert_array_equal(actual, original)


if __name__ == '__main__':
    unittest.main()
