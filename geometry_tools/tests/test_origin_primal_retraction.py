import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.origin_primal_retraction import retract_origin_feasible
class OriginRetractionTests(unittest.TestCase):
 def test_strict_original_cone_and_negative_rhs_preserved(self):
  d=np.array([[.050000005,0.,0.]]);A=sparse.csr_matrix([[-1.,0.,0.]]);C=sparse.vstack([sparse.eye(3),sparse.csr_matrix((3,3))])
  q,p=retract_origin_feasible(d,A,[-.06],maximum_motion_m=.1,distortion=C,distortion_bound=.05)
  self.assertLessEqual(np.linalg.norm(C@q.ravel()),.05);self.assertLessEqual(float(-q[0,0]),-.04999999);self.assertFalse(p['originalBoundsRelaxed']);self.assertTrue(p['applied'])
 def test_positive_origin_infeasibility_is_not_ignored(self):
  with self.assertRaisesRegex(ValueError,'zero displacement'):
   retract_origin_feasible(np.array([[.100000005,0.,0.]]),sparse.csr_matrix([[1.,0.,0.]]),[.1],maximum_motion_m=.1)
 def test_large_violation_cannot_be_hidden(self):
  with self.assertRaisesRegex(ValueError,'tiny retraction'):
   retract_origin_feasible(np.array([[.11,0.,0.]]),sparse.csr_matrix((0,3)),[],maximum_motion_m=.1)
 def test_already_feasible_point_stays_literal(self):
  d=np.array([[.01,.02,.03]]);q,p=retract_origin_feasible(d,sparse.csr_matrix((0,3)),[],maximum_motion_m=.1)
  np.testing.assert_array_equal(q,d);self.assertFalse(p['applied'])
if __name__=='__main__':unittest.main()
