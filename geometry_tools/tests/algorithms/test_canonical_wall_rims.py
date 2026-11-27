import unittest
import numpy as np
from geometry_tools.algorithms.canonical_wall_rims import canonical_wall_rims

class CanonicalWallRimsTests(unittest.TestCase):
 def test_opposite_owners_have_literal_opposed_triangles_on_twisted_quad(self):
  for edge in [[0,1],[1,0],[2,0]]:
   a=canonical_wall_rims(np.array([edge]),3);b=canonical_wall_rims(np.array([edge[::-1]]),3)
   for face in a:self.assertTrue(any(any(np.array_equal(np.roll(face,k),other[[0,2,1]])for k in range(3))for other in b))
