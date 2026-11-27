import unittest
from geometry_tools.algorithms.layered_source_routes import layered_source_routes


class LayeredSourceRoutesTest(unittest.TestCase):
 def test_collateral_stays_separate_from_two_through_running_sources(self):
  points=[[[-.002,0],[.002,0]],[[-.002,0],[.002,0],[.02,0]],[[-.002,0],[.002,0]]];links=[[(0,0),(1,1),(1,2)],[(0,0),(1,1),(2,1)]];selected,proof=layered_source_routes(points,links)
  self.assertEqual(selected,[[0,1],[0,1],[0,1]]);self.assertTrue(proof['feasible']);self.assertEqual(proof['flowResidual'],0)
 def test_two_routes_cannot_share_one_section_profile(self):
  selected,proof=layered_source_routes([[[-1,0],[1,0]],[[0,0]],[[-1,0],[1,0]]],[[(0,0),(1,0)],[(0,0),(0,1)]])
  self.assertIsNone(selected);self.assertFalse(proof['feasible'])
