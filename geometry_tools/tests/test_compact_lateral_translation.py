import unittest
import numpy as np
from geometry_tools.algorithms.compact_lateral_translation import CompactBilateralLateralTranslation

class CompactLateralTranslationTests(unittest.TestCase):
 def setUp(self):self.map=CompactBilateralLateralTranslation(.0025,(.087,.097,.123,.140),(.395,.419,.432,.455),(-.066,-.046,-.020,.002))
 def test_plateau_is_rigid_and_exterior_is_literal(self):
  p=np.asarray([[.11,.425,-.03],[-.11,.425,-.03],[.2,.425,-.03],[.11,.6,-.03]])
  q,J=self.map.apply(p);np.testing.assert_allclose(q-p,[[.0025,0,0],[-.0025,0,0],[0,0,0],[0,0,0]],atol=1e-16);np.testing.assert_array_equal(J,np.broadcast_to(np.eye(3),(4,3,3)))
 def test_analytic_jacobian_and_normal_transport_replay_finite_differences(self):
  p=np.asarray([[.092,.409,-.052],[-.092,.439,-.01],[.130,.425,-.03]]);q,J=self.map.apply(p);h=1e-7
  for j in range(3):
   step=np.zeros_like(p);step[:,j]=h;finite=(self.map.apply(p+step)[0]-self.map.apply(p-step)[0])/(2*h);np.testing.assert_allclose(finite,J[:,:,j],rtol=1e-7,atol=1e-8)
  self.assertGreater(self.map.determinant_lower_bound,0);self.assertTrue((np.linalg.det(J)>=self.map.determinant_lower_bound).all())
  n=np.tile([0.,0.,1.],(len(p),1));normal=self.map.normals(p,n);tangent=np.einsum('nij,nj->ni',J,np.tile([0.,1.,0.],(len(p),1)));np.testing.assert_allclose((normal*tangent).sum(1),0,atol=1e-14)
 def test_fold_cap_is_rejected(self):
  bad=CompactBilateralLateralTranslation(.02,self.map.absolute_x_support,self.map.y_support,self.map.z_support)
  with self.assertRaises(ValueError):bad.apply(np.asarray([[.11,.425,-.03]]))
if __name__=='__main__':unittest.main()
