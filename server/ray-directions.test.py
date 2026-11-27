from pathlib import Path
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from ta98_ray_directions import basis, coordinates, directions


class RayDirections(unittest.TestCase):
    def setUp(self):
        self.rays = np.array([
            [.001, 0, 0],
            [0, -.0007, 0],
            [.0004, -.0008, .0012],
            [-.0013, .0009, -.0011],
        ])
        self.frame = basis(self.rays)
        self.parameters = np.array([[.4, -.3], [-.8, .2], [1.2, .9], [-.2, -.7]])

    def test_zero_thickness_is_rejected(self):
        for rays in (np.zeros((1, 3)), np.vstack((self.rays, [0, 0, 0]))):
            with self.subTest(rays=rays.tolist()):
                with self.assertRaisesRegex(ValueError, 'nonzero thickness'):
                    basis(rays)

    def test_zero_parameters_preserve_original_rays(self):
        rays, jacobian = directions(np.zeros_like(self.parameters), self.frame)
        np.testing.assert_allclose(rays, self.rays, rtol=2e-15, atol=1e-18)
        np.testing.assert_allclose(coordinates(self.rays, self.frame), 0, atol=2e-16)
        self.assertEqual(jacobian.shape, (len(self.rays), 3, 2))

    def test_nontrivial_rotations_preserve_length_and_round_trip(self):
        rays, jacobian = directions(self.parameters, self.frame)
        lengths, e, _, _ = self.frame
        np.testing.assert_allclose(np.linalg.norm(rays, axis=1), lengths, rtol=2e-15)
        np.testing.assert_allclose(coordinates(rays, self.frame), self.parameters,
                                   rtol=3e-15, atol=3e-15)
        expected_angle = 2 * np.arctan(np.linalg.norm(self.parameters, axis=1))
        actual_angle = np.arccos(np.clip(np.einsum('vi,vi->v', rays, e) / lengths, -1, 1))
        np.testing.assert_allclose(actual_angle, expected_angle, atol=2e-15)
        # Tangent columns must preserve length and have equal, orthogonal scales.
        np.testing.assert_allclose(np.einsum('vi,vij->vj', rays, jacobian), 0, atol=2e-21)
        scale = 2 * lengths / (1 + np.sum(self.parameters**2, axis=1))
        metric = np.einsum('via,vib->vab', jacobian, jacobian) / scale[:, None, None]**2
        np.testing.assert_allclose(metric, np.broadcast_to(np.eye(2), metric.shape),
                                   rtol=2e-15, atol=2e-15)

    def test_derivatives_match_central_differences(self):
        _, jacobian = directions(self.parameters, self.frame)
        step = 1e-6
        for axis in range(2):
            perturbation = np.zeros_like(self.parameters)
            perturbation[:, axis] = step
            numerical = (directions(self.parameters + perturbation, self.frame)[0]
                         - directions(self.parameters - perturbation, self.frame)[0]) / (2 * step)
            np.testing.assert_allclose(jacobian[..., axis], numerical, rtol=2e-8, atol=3e-13)

    def test_antipode_and_near_pole_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'chart pole'):
            coordinates(-self.rays, self.frame)
        lengths, e, t1, _ = self.frame
        almost_pole = -e + 1e-7 * t1
        almost_pole *= (lengths / np.linalg.norm(almost_pole, axis=1))[:, None]
        with self.assertRaisesRegex(ValueError, 'chart pole'):
            coordinates(almost_pole, self.frame)


if __name__ == '__main__':
    unittest.main()
