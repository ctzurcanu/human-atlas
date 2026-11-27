import unittest
import numpy as np
from geometry_tools.tests.test_declared_mouth_batch import DeclaredMouthBatchTests
from geometry_tools.algorithms.fixed_declared_mouth_area_budgets import reconcile_declared_mouths

class FixedMouthTests(unittest.TestCase):
 def solve(self,meshes,pairs,budget=25e-6):
  return reconcile_declared_mouths(meshes,pairs,maximum_gap_m=3e-5,maximum_motion_m=2e-5,maximum_fixed_rim_motion_m=budget,maximum_relative_area_change=.025,maximum_source_area_discrepancy=.05,maximum_facet_rotation_degrees=10)
 def test_fixed_first_side_is_literal_and_other_mouth_is_complete(self):
  meshes,pair=DeclaredMouthBatchTests().fixtures(gap=23e-6);pair['fixedOwner']='first';result,proof=self.solve(meshes,[pair])
  np.testing.assert_array_equal(result['first'].vertices,meshes['first'].vertices)
  np.testing.assert_array_equal(result['second'].vertices[:6],meshes['first'].vertices[6:])
  np.testing.assert_array_equal(result['second'].vertices[6:],meshes['second'].vertices[6:])
  for key in meshes:np.testing.assert_array_equal(result[key].faces,meshes[key].faces)
  self.assertEqual(proof['lumenDisksCreated'],0)
 def test_fixed_second_side_is_literal(self):
  meshes,pair=DeclaredMouthBatchTests().fixtures();pair['fixedOwner']='second';result,_=self.solve(meshes,[pair])
  np.testing.assert_array_equal(result['second'].vertices,meshes['second'].vertices)
  np.testing.assert_array_equal(result['first'].vertices[6:],meshes['second'].vertices[:6])
 def test_pin_budget_does_not_silently_override_other_constraints(self):
  meshes,pair=DeclaredMouthBatchTests().fixtures(gap=23e-6);pair['fixedOwner']='first'
  with self.assertRaisesRegex(ValueError,'movement'):self.solve(meshes,[pair],budget=20e-6)
 def test_unknown_fixed_owner_rejected(self):
  meshes,pair=DeclaredMouthBatchTests().fixtures();pair['fixedOwner']='unrelated'
  with self.assertRaisesRegex(ValueError,'fixed owner'):self.solve(meshes,[pair])
 def test_no_separate_budget_rejected(self):
  meshes,pair=DeclaredMouthBatchTests().fixtures();pair['fixedOwner']='first'
  with self.assertRaisesRegex(ValueError,'budget'):self.solve(meshes,[pair],budget=None)
if __name__=='__main__':unittest.main()
