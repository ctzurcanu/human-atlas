import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.algorithms.compact_similarity_flow import CompactSimilarityFlow

class SimilarityFlowTests(unittest.TestCase):
 def field(self):
  c=np.array([.03,1.2,-.02]);a=1.08*Rotation.from_rotvec([.03,-.04,.02]).as_matrix();m=np.eye(4);m[:3,:3]=a;m[:3,3]=c+np.array([.001,-.01,.004])-a@c
  return CompactSimilarityFlow(m,center=c,support_axes_m=[.08,.18,.12],core_fraction=.65,steps=64)

 def test_full_similarity_in_core_and_literal_identity_outside(self):
  f=self.field();p=f.center+np.array([[0,0,0],[.002,.003,-.002],[.2,0,0]]);q,j,proof=f.apply(p)
  np.testing.assert_allclose(q[:2],p[:2]@f.matrix[:3,:3].T+f.matrix[:3,3],atol=5e-15,rtol=0);np.testing.assert_allclose(j[:2],np.broadcast_to(f.matrix[:3,:3],(2,3,3)),atol=1e-13,rtol=0);np.testing.assert_array_equal(q[2],p[2]);np.testing.assert_array_equal(j[2],np.eye(3));self.assertGreater(proof['analytic_determinant_lower_bound'],0)

 def test_analytic_jacobian_and_positive_determinants_in_transition(self):
  f=self.field();rng=np.random.default_rng(73);u=rng.normal(size=(30,3));u/=np.linalg.norm(u,axis=1)[:,None];p=f.center+u*rng.uniform(.68,.98,(30,1))*f.axes;original=p.copy();q,j,proof=f.apply(p);eps=1e-7
  for axis in range(3):
   d=np.eye(3)[axis]*eps;finite=(f.apply(p+d)[0]-f.apply(p-d)[0])/(2*eps);np.testing.assert_allclose(finite,j[:,:,axis],atol=2e-7,rtol=2e-7)
  self.assertTrue(np.all(np.linalg.det(j)>=proof['analytic_determinant_lower_bound']));np.testing.assert_array_equal(p,original);self.assertTrue(np.all(np.linalg.norm((q-f.center)/f.axes,axis=1)<1))

 def test_repeated_shared_controls_and_invalid_transforms(self):
  f=self.field();p=f.center+np.array([[.04,.01,.01]]);a=f.apply(p)[0];b=f.apply(np.repeat(p,3,axis=0))[0];np.testing.assert_array_equal(b,np.repeat(a,3,axis=0))
  for diagonal in [[-1,1,1,1],[1,2,1,1]]:
   with self.assertRaises(ValueError):CompactSimilarityFlow(np.diag(diagonal),center=f.center,support_axes_m=f.axes)
  m=np.eye(4);m[:3,3]=[1,0,0]
  with self.assertRaises(ValueError):CompactSimilarityFlow(m,center=[0,0,0],support_axes_m=[.01,.01,.01],steps=1)

if __name__=='__main__':unittest.main()
