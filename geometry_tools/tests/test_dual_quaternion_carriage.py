import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.algorithms.dual_quaternion_carriage import blend_proper_poses


class DualQuaternionCarriageTests(unittest.TestCase):
    def test_single_owner_replays_rigid_pose(self):
        p=np.array([[.1,.3,-.2],[.2,.7,-.8]])
        R=Rotation.from_euler('xyz',[[12,61,-90],[-70,23,81]],degrees=True).as_matrix()
        t=np.array([[.3,.4,.1],[-.7,.1,.8]])
        q,_=blend_proper_poses(p,np.eye(2),R,t)
        np.testing.assert_allclose(q,np.einsum('nij,nj->ni',R,p)+t,atol=3e-16)

    def test_opposing_twists_preserve_radius_and_translation(self):
        p=np.array([[1.,0,0],[0,1.,0]])
        R=Rotation.from_euler('z',[[-80],[80]],degrees=True).as_matrix()
        q,poses=blend_proper_poses(p,np.full((2,2),.5),R,np.tile([.3,.7,.4],(2,1)))
        np.testing.assert_allclose(q,p+[.3,.7,.4],atol=3e-16)
        np.testing.assert_allclose(np.linalg.det(poses),1,atol=3e-16)

    def test_permutation_does_not_change_blend(self):
        R=Rotation.from_euler('z',[[175],[-170],[12]],degrees=True).as_matrix()
        p=np.array([[1.,2.,3.],[-.4,.1,.6]])
        w=np.array([[.8,.15,.05],[.1,.7,.2]]);t=np.arange(9).reshape(3,3)*.1
        q,_=blend_proper_poses(p,w,R,t)
        perm=[2,0,1];z,_=blend_proper_poses(p,w[:,perm],R[perm],t[perm])
        np.testing.assert_allclose(q,z,atol=1e-15)

    def test_reflection_and_bad_weights_rejected(self):
        with self.assertRaises(ValueError):
            blend_proper_poses([[0,0,0]],[[1]],np.diag([-1,1,1])[None],[[0,0,0]])
        with self.assertRaises(ValueError):
            blend_proper_poses([[0,0,0]],[[.9]],np.eye(3)[None],[[0,0,0]])


if __name__=='__main__':unittest.main()
