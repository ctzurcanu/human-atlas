import unittest
import numpy as np
from geometry_tools.algorithms.triangle_cell_halfspaces import triangle_cell_halfspaces, moving_triangle_cell_halfspaces


class TriangleCellHalfspacesTest(unittest.TestCase):
 def test_plane_blocks_crossing_without_moving_fixed_triangle(self):
  first=np.array([[[0,0,0],[.002,0,0],[0,.002,0]]]);second=first+[0,0,.001];rows,protected,proof=triangle_cell_halfspaces(first,[[0,0,1]],second)
  self.assertFalse(protected);self.assertEqual(len(rows),2)
  for row in rows:self.assertGreaterEqual(np.dot(row['normalXZ'],[0,.0002]),row['minimum_m']);self.assertLess(np.dot(row['normalXZ'],[0,.0009]),row['minimum_m'])
  self.assertAlmostEqual(proof[0]['supportGap_m'],.001)
 def test_existing_intersection_is_protected_not_called_repaired(self):
  first=np.array([[[0,0,0],[.002,0,0],[0,.002,0]]]);rows,protected,proof=triangle_cell_halfspaces(first,[[1,2,3]],first)
  self.assertFalse(rows);self.assertEqual(protected,{1,2,3});self.assertTrue(proof[0]['protected'])
 def test_joint_motion_preserves_gap_but_relative_crossing_fails(self):
  first=np.array([[[0,0,0],[.002,0,0],[0,.002,0]]]);second=first+[0,0,.001]
  rows,unresolved,proof=moving_triangle_cell_halfspaces(first,[[0,1,2]],second,[[3,4,5]])
  self.assertFalse(unresolved);self.assertEqual(len(rows),9)
  def valid(delta):
   return all(sum(np.dot(value,delta[int(key)])for key,value in row['coefficients'].items())>=row['minimum_m']-1e-12 for row in rows)
  common=np.tile([.004,.003],(6,1));self.assertTrue(valid(common))
  common[:3,1]+=.002;self.assertFalse(valid(common))
  self.assertAlmostEqual(proof[0]['supportGap_m'],.001)
 def test_joint_existing_crossing_cannot_claim_a_separating_plane(self):
  triangle=np.array([[[0,0,0],[.002,0,0],[0,.002,0]]])
  rows,unresolved,proof=moving_triangle_cell_halfspaces(triangle,[[0,1,2]],triangle,[[3,4,5]])
  self.assertEqual(unresolved,[0]);self.assertFalse(rows);self.assertFalse(proof)


if __name__=='__main__':unittest.main()
