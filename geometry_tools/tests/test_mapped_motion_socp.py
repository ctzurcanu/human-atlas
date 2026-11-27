import unittest
import numpy as np
from geometry_tools.algorithms.section_course_socp import solve_weighted_section_socp


class MappedMotionTests(unittest.TestCase):
    def solve(self,bound):
        return solve_weighted_section_socp(np.eye(2),np.array([[1.,1.]]),np.array([[1.,0,0]]),np.array([.002]),maximum_motion_m=bound,motion_weights=np.array([[1.,1.],[1.,-1.]]))

    def test_physical_station_bound_is_not_a_parameter_bound(self):
        # Each parameter can be only 1 mm while their source-station sum is 2.
        d,p=self.solve(.0015);self.assertIsNone(d);self.assertFalse(p['feasible'])

    def test_actual_source_motion_and_lower_bound_are_replayed(self):
        d,p=self.solve(.003);self.assertTrue(p['feasible']);self.assertAlmostEqual(p['maximumActualMotionMeters'],.002,places=8);self.assertAlmostEqual(p['necessaryRadiusLowerBoundMeters'],.002,places=8)

    def test_unbounded_parameter_direction_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'every movable'):
            solve_weighted_section_socp(np.eye(2),np.empty((0,2)),np.empty((0,3)),np.empty(0),maximum_motion_m=.002,motion_weights=np.array([[1.,1.]]))


if __name__=='__main__':unittest.main()
