import unittest
import numpy as np
from geometry_tools.algorithms.surface_coverage_bound import symmetric_surface_bound

class SurfaceCoverage(unittest.TestCase):
    def test_parallel_surfaces_are_fully_bounded(self):
        a=np.array([[[0.,0,0],[1,0,0],[0,1,0]]]);b=a+np.array([0,0,.02])
        p=symmetric_surface_bound(a,b,maximum_distance_m=.06)
        self.assertGreaterEqual(p['symmetricDistanceUpperBound_m'],.02)
        self.assertLessEqual(p['symmetricDistanceUpperBound_m'],.06)
        self.assertGreater(p['directions'][0]['queriedPoints'],1)

    def test_displaced_surface_is_rejected(self):
        a=np.array([[[0.,0,0],[1,0,0],[0,1,0]]])
        with self.assertRaisesRegex(ValueError,'exceed'):symmetric_surface_bound(a,a+[0,0,.2],maximum_distance_m=.1)

    def test_budget_is_enforced_and_rigid_world_frame_is_supported(self):
        a=np.array([[[0.,0,0],[1,0,0],[0,1,0]]])
        with self.assertRaisesRegex(ValueError,'budget'):symmetric_surface_bound(a,a,maximum_distance_m=.0001,maximum_points=1)
        R=np.array([[0.,-1,0],[0,0,-1],[1,0,0]]);shift=np.array([2.,4.,-3.])
        p=symmetric_surface_bound(a@R.T+shift,(a+[0,0,.02])@R.T+shift,maximum_distance_m=.06)
        self.assertLessEqual(p['symmetricDistanceUpperBound_m'],.06)

if __name__=='__main__':unittest.main()
