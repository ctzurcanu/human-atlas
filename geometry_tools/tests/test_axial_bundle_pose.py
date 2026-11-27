import unittest
import numpy as np
from geometry_tools.algorithms.axial_bundle_pose import AxialBundlePose

class AxialBundlePoseTests(unittest.TestCase):
 def field(self):
  controls=dict(y_m=[0.,1.],sourcePairMidpointXZ_m=[[0,0],[0,0]],targetPairMidpointXZ_m=[[.01,-.02],[.01,-.02]],angleRadians=[.5,.5])
  return AxialBundlePose(controls,[0,.1,.9,1],{'pair':[[.25,.26]],'joint':[[.5,.51]]},transition_m=.05)
 def test_separate_artery_and_veins_share_the_same_rigid_pose(self):
  f=self.field();R,t,w=f.pose(.6,['pair']);p=np.array([[0,.6,0],[.002,.6,0],[0,.6,.003],[.001,.6,.003]])
  q=p@R.T+t;np.testing.assert_allclose(np.linalg.norm(p[:,None]-p[None,:],axis=2),np.linalg.norm(q[:,None]-q[None,:],axis=2),atol=1e-15);self.assertAlmostEqual(np.linalg.det(R),1);self.assertEqual(w,1)
 def test_fixed_source_interval_has_identity_pose(self):
  R,t,w=self.field().pose(.255,['pair']);np.testing.assert_array_equal(R,np.eye(3));np.testing.assert_allclose(t,0);self.assertEqual(w,0)
 def test_joint_honours_both_explicit_bundle_anchors(self):
  R,t,w=self.field().pose(.505,['pair','joint']);np.testing.assert_array_equal(R,np.eye(3));np.testing.assert_allclose(t,0);self.assertEqual(w,0)
 def test_transition_is_proper_and_keeps_axial_coordinate(self):
  f=self.field();R,t,w=f.pose(.285,['pair']);self.assertGreater(w,0);self.assertLess(w,1);p=np.array([.02,.285,.03]);self.assertAlmostEqual((R@p+t)[1],p[1]);np.testing.assert_allclose(R.T@R,np.eye(3),atol=1e-15)
 def test_unknown_bundle_is_rejected(self):
  with self.assertRaises(ValueError):self.field().pose(.5,['unknown'])

if __name__=='__main__':unittest.main()
