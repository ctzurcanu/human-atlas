import unittest
import numpy as np
from geometry_tools.algorithms.triangle_ray_intervals import first_positive_triangle_hits
class TriangleRayIntervalsTest(unittest.TestCase):
 def test_nearest_actual_surface_and_scaled_direction(self):
  t=np.array([[[-1,-1,2],[1,-1,2],[0,1,2]],[[-1,-1,1],[1,-1,1],[0,1,1]]],float);r=first_positive_triangle_hits([[0,0,0]],[[0,0,7]],t,maximum_distance_m=3);self.assertEqual(r['face_ids'][0],1);self.assertAlmostEqual(r['distance_m'][0],1);np.testing.assert_allclose(r['points'][0],[0,0,1]);self.assertLess(r['ray_residual_m'][0],1e-14)
 def test_miss_and_finite_local_travel_are_explicit(self):
  t=np.array([[[-1,-1,2],[1,-1,2],[0,1,2]]],float);r=first_positive_triangle_hits([[0,0,0],[5,5,0]],[[0,0,1],[0,0,1]],t,maximum_distance_m=1);self.assertFalse(r['hit'].any());self.assertTrue(np.isnan(r['distance_m']).all());self.assertTrue(np.all(r['face_ids']==-1))
 def test_negative_and_parallel_surfaces_are_not_hits(self):
  t=np.array([[[-1,-1,-1],[1,-1,-1],[0,1,-1]]],float);r=first_positive_triangle_hits([[0,0,0],[0,0,0]],[[0,0,1],[1,0,0]],t,maximum_distance_m=3);self.assertFalse(r['hit'].any())
if __name__=='__main__':unittest.main()
