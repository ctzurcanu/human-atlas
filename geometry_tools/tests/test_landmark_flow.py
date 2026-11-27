import unittest

import numpy as np

from geometry_tools.algorithms.landmark_flow import landmark_flow


class LandmarkFlowTests(unittest.TestCase):
    def test_translation_and_derivative(self):
        p = np.array([[0., 0., 0.], [.1, 0., 0.], [0., .1, 0.], [0., 0., .1], [.1, .1, .1]])
        q, J, proof = landmark_flow(p, p, p+[.003, -.002, .001], steps=4)
        np.testing.assert_allclose(q, p+[.003, -.002, .001], atol=1e-12)
        np.testing.assert_allclose(J, np.broadcast_to(np.eye(3), J.shape), atol=1e-12)
        self.assertLess(proof['landmarkError_m'], 1e-12)

    def test_fixed_proximal_region(self):
        c = np.array([[0., 0., 0.], [.1, 0., 0.], [0., .1, 0.], [0., 0., .1], [.1, .1, .1]])
        p = np.r_[c, [[.02, 1.1, .03]]]
        q, J, _ = landmark_flow(p, c, c+[.001, 0., 0.], steps=4)
        np.testing.assert_array_equal(q[-1], p[-1])
        np.testing.assert_array_equal(J[-1], np.eye(3))

    def test_rejects_unbounded_and_repeated_controls(self):
        c = np.array([[0., 0., 0.], [.1, 0., 0.], [0., .1, 0.], [0., 0., .1]])
        for kwargs in [dict(steps=0), dict(wall_seconds=0), dict(steps=65)]:
            with self.assertRaises(ValueError):
                landmark_flow(c, c, c, **kwargs)
        with self.assertRaises(ValueError):
            landmark_flow(c, np.r_[c, c[:1]], np.r_[c, c[:1]])


if __name__ == '__main__':
    unittest.main()
