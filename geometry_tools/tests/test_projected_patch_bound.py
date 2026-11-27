import unittest

import numpy as np

from geometry_tools.algorithms.projected_patch_bound import projected_patch_bound


class ProjectedPatchBoundTests(unittest.TestCase):
    def fixture(self):
        return (np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,.01],[0.,1.,0.]])*.001,
                np.array([[0,1,2],[0,2,3]]), np.array([[0,1,3],[1,2,3]]))

    def test_complete_overlay_bounds_both_diagonal_surfaces(self):
        p, a, b = self.fixture()
        proof = projected_patch_bound(p, a, b, maximum_distance_m=6e-6)
        self.assertTrue(proof['completeOverlayCoverageChecked'])
        self.assertGreater(proof['maximumTwoWaySurfaceDistanceUpperBound_m'], 4.9e-6)
        self.assertLess(proof['maximumTwoWaySurfaceDistanceUpperBound_m'], 5.1e-6)
        reverse = projected_patch_bound(p, b, a, maximum_distance_m=6e-6)
        self.assertAlmostEqual(proof['maximumOverlayHeightDifference_m'],
                               reverse['maximumOverlayHeightDifference_m'], places=14)

    def test_distance_limit_cannot_be_silently_relaxed(self):
        p, a, b = self.fixture()
        with self.assertRaisesRegex(ValueError, 'exceeds'):
            projected_patch_bound(p, a, b, maximum_distance_m=4e-6)

    def test_incomplete_or_inconsistently_oriented_patch_rejected(self):
        p, a, b = self.fixture()
        with self.assertRaises(ValueError):
            projected_patch_bound(p, a, b[:1], maximum_distance_m=6e-6)
        b[1] = b[1, ::-1]
        with self.assertRaises(ValueError):
            projected_patch_bound(p, a, b, maximum_distance_m=6e-6)


if __name__ == '__main__':
    unittest.main()
