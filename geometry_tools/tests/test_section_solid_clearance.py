import unittest
import numpy as np
from geometry_tools.algorithms.section_solid_clearance import translated_section_clearance

class SectionSolidClearanceTests(unittest.TestCase):
 def square(self,r):return np.array([[-r,-r],[r,-r],[r,r],[-r,r]])
 def test_small_nerve_enclosed_despite_clear_vessel_vertices_is_rejected(self):
  result=translated_section_clearance([self.square(.002)],[[0,0]],[self.square(.0003)],margin_m=.0002);self.assertFalse(result[0])
 def test_edges_cross_with_no_contained_vertices(self):
  horizontal=np.array([[-.003,-.0002],[.003,-.0002],[.003,.0002],[-.003,.0002]]);vertical=horizontal[:,::-1]
  self.assertFalse(translated_section_clearance([horizontal],[[0,0]],[vertical],margin_m=0)[0])
 def test_actual_gap_is_enforced_for_each_translated_bundle(self):
  result=translated_section_clearance([self.square(.001)],[[.0022,0],[.003,0]],[self.square(.001)],margin_m=.0003);np.testing.assert_array_equal(result,[False,True])
 def test_invalid_contour_is_never_repaired(self):
  bow=np.array([[-1,-1],[1,1],[-1,1],[1,-1]])
  with self.assertRaises(ValueError):translated_section_clearance([bow],[[0,0]],[self.square(.1)],margin_m=0)

if __name__=='__main__':unittest.main()
