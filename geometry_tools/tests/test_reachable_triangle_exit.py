import unittest
import numpy as np
from geometry_tools.algorithms.reachable_triangle_exit import reachable_exit_axes


class ReachableTriangleExitTests(unittest.TestCase):
    def test_fixed_overlapping_triangles_cannot_have_an_exit(self):
        a = np.array([[[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.]]])
        _, available, proof = reachable_exit_axes(a, a, np.full((1, 3), -1), np.full((1, 3), -1),
                                                  maximum_motion_m=.0001, gap_m=.000001)
        self.assertFalse(available[0])
        self.assertEqual(proof['pairsUnreachableWithinDeclaredStations'], 1)

    def test_a_rigid_triangle_can_exit_through_its_normal(self):
        a = np.array([[[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.]]])
        axes, available, _ = reachable_exit_axes(a, a, np.zeros((1, 3), int), np.full((1, 3), -1),
                                               maximum_motion_m=.0001, gap_m=.000001)
        self.assertTrue(available[0])
        self.assertAlmostEqual(abs(axes[0, 2]), 1.)

    def test_a_common_station_cannot_separate_its_own_faces(self):
        a = np.array([[[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.]]])
        _, available, _ = reachable_exit_axes(a, a, np.zeros((1, 3), int), np.zeros((1, 3), int),
                                             maximum_motion_m=.001, gap_m=.000001)
        self.assertFalse(available[0])


if __name__ == '__main__':
    unittest.main()
