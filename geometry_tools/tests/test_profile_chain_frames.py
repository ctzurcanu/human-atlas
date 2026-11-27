import unittest
import numpy as np
from geometry_tools.algorithms.profile_chain_frames import branch_anchored_frames


class ProfileChainFramesTest(unittest.TestCase):
    def test_unchanged_chain_is_exactly_identity(self):
        p=np.array([[0.,0,0],[1,0,0],[1.5,.2,0],[2,.8,.1]])
        np.testing.assert_array_equal(branch_anchored_frames(p,p,np.eye(3),np.eye(3)),np.repeat(np.eye(3)[None],4,axis=0))

    def test_common_rigid_pose_retains_oblique_section(self):
        p=np.array([[0.,0,0],[1,0,0],[1.5,.2,0],[2,.8,.1]])
        R=np.array([[0.,-1,0],[1.,0,0],[0,0,1]])
        rotations=branch_anchored_frames(p,p@R.T+[2,3,4],R,R)
        np.testing.assert_allclose(rotations,np.repeat(R[None],4,axis=0),atol=1e-12)
        oblique=np.array([.2,.1,.4])
        np.testing.assert_allclose(rotations@oblique,np.repeat((R@oblique)[None],4,axis=0),atol=1e-12)

    def test_distinct_prescribed_endpoint_roll_is_exact(self):
        p=np.array([[0.,0,0],[1,0,0],[2,0,0],[3,0,0]])
        end=np.array([[1.,0,0],[0,0,-1],[0,1,0]])
        r=branch_anchored_frames(p,p,np.eye(3),end)
        np.testing.assert_array_equal(r[0],np.eye(3))
        np.testing.assert_array_equal(r[-1],end)
        np.testing.assert_allclose(r.transpose(0,2,1)@r,np.repeat(np.eye(3)[None],4,axis=0),atol=1e-12)


if __name__=='__main__':unittest.main()
