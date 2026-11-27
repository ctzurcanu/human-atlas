import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.existing_section_side_axes import existing_section_side_axes
class ExistingSideTests(unittest.TestCase):
 def model(self,centers):
  ring=np.array([[.0001,0,0],[0,.0001,0],[-.0001,0,0],[0,-.0001,0]])
  return Mesh(np.concatenate([ring+c for c in centers]),np.array([[4*i,4*i+1,4*i+2] for i in range(len(centers))]),'test:metres')
 def test_opposite_existing_sides_are_preserved(self):
  a=self.model([[-.001,0,0],[.001,0,0]]);b=self.model([[0,0,0],[0,0,0]]);g=np.repeat([0,1],4);axes,p=existing_section_side_axes(a,b,g,g,[[0,0],[1,1]])
  self.assertEqual(p['resolvedPairs'],2);np.testing.assert_allclose(axes,[[-1,0,0],[1,0,0]],atol=1e-12)
 def test_coincident_profiles_are_explicitly_unresolved(self):
  a=self.model([[0,0,0]]);axes,p=existing_section_side_axes(a,a,np.zeros(4,int),np.zeros(4,int),[[0,0]])
  self.assertEqual(p['unresolvedPairs'],1);self.assertFalse(p['resolvedMask'][0]);np.testing.assert_array_equal(axes,[[0,0,0]])
if __name__=='__main__':unittest.main()
