import unittest
import numpy as np
from geometry_tools.algorithms.planar_sliver_flip import flip_planar_slivers

class PlanarSliverFlipTests(unittest.TestCase):
 def test_planar_diagonal_improves_altitude_without_moving_controls(self):
  p=np.array([[0.,0,0],[1,0,0],[2,1e-8,0],[0,1,0]]);f=np.array([[0,1,2],[0,2,3]]);before=p.copy();out,report=flip_planar_slivers(p,f,minimum_altitude_m=1e-6)
  np.testing.assert_array_equal(p,before);self.assertEqual(len(report['flips']),1);self.assertEqual(report['remaining_below_threshold'],0);self.assertEqual(len(out),2)
 def test_declared_plane_bound_is_respected(self):
  p=np.array([[0.,0,0],[1,0,0],[2,1e-8,0],[0,1,.001]]);f=np.array([[0,1,2],[0,2,3]]);out,report=flip_planar_slivers(p,f,minimum_altitude_m=1e-6,plane_tolerance_m=1e-14)
  np.testing.assert_array_equal(out,f);self.assertEqual(len(report['flips']),0)
