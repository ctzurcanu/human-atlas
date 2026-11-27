import unittest

import numpy as np

from geometry_tools.algorithms.axial_endpoint_field import endpoint_translation


class EndpointFieldTests(unittest.TestCase):
    def test_parallel_profiles_and_fixed_controls(self):
        section = np.array([[0., 1., 0.], [0., 0., 1.], [0., -1., 0.]])
        points = np.concatenate([section + [x, 0., 0.] for x in [-1., 0., .3, .8, 1., 2.]])
        normals = np.tile([0., 1., 0.], (len(points), 1))
        moved, transported, proof = endpoint_translation(points, normals, [1., 0., 0.], 0., 1., [.1, .2, 0.])
        for a in range(0, len(points), 3):
            np.testing.assert_allclose(moved[a:a+3] - moved[a], points[a:a+3] - points[a], atol=5e-16)
        np.testing.assert_array_equal(moved[:6], points[:6])
        np.testing.assert_array_equal(transported[:6], normals[:6])
        np.testing.assert_allclose(moved[-3:] - points[-3:], np.tile([.1, .2, 0.], (3, 1)))
        self.assertTrue(proof['fixed_controls_literal'])

    def test_normal_is_orthogonal_to_independently_differenced_tangent(self):
        point = np.array([[.42, .3, .2]])
        normal = np.array([[0., 1., 0.]])
        args = ([1., 0., 0.], 0., 1., [.1, .2, -.1])
        _, transported, _ = endpoint_translation(point, normal, *args)
        epsilon = 1e-6
        for tangent in [[1., 0., 0.], [0., 0., 1.]]:
            tangent = np.array(tangent)
            plus = endpoint_translation(point + epsilon*tangent, normal, *args)[0]
            minus = endpoint_translation(point - epsilon*tangent, normal, *args)[0]
            actual_tangent = (plus - minus)[0] / (2*epsilon)
            self.assertLess(abs(transported[0] @ actual_tangent), 1e-10)

    def test_folding_field_is_rejected_before_sampling(self):
        with self.assertRaisesRegex(ValueError, 'complete continuous'):
            endpoint_translation(np.array([[0., 0., 0.]]), np.array([[0., 1., 0.]]),
                                 [1., 0., 0.], 0., 1., [-1., 0., 0.])


if __name__ == '__main__':
    unittest.main()
