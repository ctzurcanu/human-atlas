import unittest
import numpy as np
from geometry_tools.algorithms.observed_profile_carrier import ObservedProfileCarrier,rigid_profile_observations


class ObservedProfileCarrierTests(unittest.TestCase):
    def test_rigid_observations_recover_exact_pose_and_keep_profile_dimensions(self):
        ring=np.array([[1.,0,0],[0,1,0],[-1,0,0],[0,-1,0]])*.001
        source=np.array([ring,ring+[0,0,.01]])
        rotation=np.array([[0.,-1,0],[1,0,0],[0,0,1]])
        target=source@rotation.T+[.02,.03,0]
        c,R,t,error=rigid_profile_observations(source,target)
        self.assertLess(float(error.max()),1e-14)
        points=np.array([[.001,.002,.005],[.002,.001,.008]])
        q,poses,_=ObservedProfileCarrier(c,R,t,[5,5]).poses(points,[5,5])
        np.testing.assert_allclose(q,points@rotation.T+[.02,.03,0],atol=1e-14)
        np.testing.assert_allclose(poses,[rotation,rotation],atol=1e-14)

    def test_nearby_other_digit_never_contaminates_pose(self):
        centres=np.array([[0.,0,0],[.0001,0,0]])
        carrier=ObservedProfileCarrier(centres,np.repeat(np.eye(3)[None],2,axis=0),[[.001,0,0],[-.001,0,0]],[4,5])
        q,_,_=carrier.poses(np.array([[.00005,0,0],[.00005,0,0]]),[4,5])
        np.testing.assert_allclose(q[:,0],[.00105,-.00095],atol=1e-15)
        with self.assertRaises(ValueError):carrier.poses([[0.,0,0]],[3])
        with self.assertRaises(ValueError):carrier.poses([[.1,0,0]],[4])


if __name__=='__main__':unittest.main()
