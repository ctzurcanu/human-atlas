import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.algorithms.profile_pose_transfer import transfer_profile_poses


class ProfilePoseTransferTest(unittest.TestCase):
 def fixture(self):
  t=np.arange(6)*np.pi/3
  return np.c_[.001*np.cos(t),.001*np.sin(t),np.zeros(6)]

 def test_reflected_direction_map_keeps_recipient_shape_and_proper_rotation(self):
  p=self.fixture();R=Rotation.from_rotvec([.2,.1,.3]).as_matrix();d=np.array([.002,.001,.003]);M=np.diag([-1,1,1]);recipient=p@M.T+[.1,.2,.3]
  r=transfer_profile_poses(p,p@R.T+d,recipient,[np.arange(6)],M,maximum_source_residual_m=1e-10)
  self.assertAlmostEqual(np.linalg.det(r['rotations'][0]),1)
  q=(recipient-recipient.mean(0))@r['rotations'][0].T+r['target_centres'][0]
  np.testing.assert_allclose(q,(p@R.T+d)@M.T+[.1,.2,.3],atol=1e-15)
  np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None],axis=2),np.linalg.norm(recipient[:,None]-recipient[None],axis=2),atol=1e-15)

 def test_recipient_translation_and_distinct_shape_are_retained(self):
  p=self.fixture();c=p*np.array([1.2,.8,1])+[1,2,3]
  r=transfer_profile_poses(p,p+[.01,0,0],c,[np.arange(6)],np.eye(3),maximum_source_residual_m=1e-10)
  np.testing.assert_allclose(r['target_centres'][0],c.mean(0)+[.01,0,0]);np.testing.assert_allclose(r['rotations'][0],np.eye(3),atol=1e-14)

 def test_rejects_stretched_source_overlap_and_scaled_direction_map(self):
  p=self.fixture()
  with self.assertRaises(ValueError):transfer_profile_poses(p,2*p,p,[np.arange(6)],np.eye(3),maximum_source_residual_m=1e-8)
  with self.assertRaises(ValueError):transfer_profile_poses(p,p,p,[np.arange(6),np.arange(6)],np.eye(3),maximum_source_residual_m=1e-8)
  with self.assertRaises(ValueError):transfer_profile_poses(p,p,p,[np.arange(6)],2*np.eye(3),maximum_source_residual_m=1e-8)


if __name__=='__main__':unittest.main()
