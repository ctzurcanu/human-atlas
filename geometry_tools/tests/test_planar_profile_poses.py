import unittest
import numpy as np
import shapely
from geometry_tools.algorithms.planar_profile_poses import profile_pose_mask


class PlanarProfilePosesTest(unittest.TestCase):
 def test_rotation_can_fit_without_scaling_the_profile(self):
  source=np.array([[-.002,-.0004],[.002,-.0004],[.002,.0004],[-.002,.0004]])
  enclosure=shapely.box(-.001,-.004,.001,.004)
  result=profile_pose_mask(source,[[0,0,0],[0,0,np.pi/2]],container=enclosure)
  self.assertEqual(result['clear'].tolist(),[False,True]);self.assertLess(result['area_relative_error'].max(),1e-12)
 def test_open_sheet_crossing_rejected_even_with_clear_vertices(self):
  source=np.array([[-.001,-.001],[.001,-.001],[.001,.001],[-.001,.001]])
  barrier=shapely.linestrings([[-.002,0],[.002,0]])
  self.assertFalse(profile_pose_mask(source,[[0,0,0]],barrier=barrier)['clear'][0])
 def test_inside_a_solid_and_crossing_an_enclosure_are_distinct_failures(self):
  source=np.array([[-.001,-.001],[.001,-.001],[.001,.001],[-.001,.001]])
  self.assertFalse(profile_pose_mask(source,[[0,0,0]],solid=shapely.box(-.002,-.002,.002,.002))['clear'][0])
  self.assertFalse(profile_pose_mask(source,[[.002,0,0]],container=shapely.box(-.002,-.002,.002,.002))['clear'][0])
