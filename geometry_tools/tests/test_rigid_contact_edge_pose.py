import unittest
import numpy as np
from geometry_tools.algorithms.rigid_contact_edge_pose import contact_edge_pivot,pose_about_contact_edge


class ContactEdgePoseTests(unittest.TestCase):
    def profile(self):
        t=np.arange(6)*np.pi/3
        return np.c_[.001*np.cos(t),.001*np.sin(t),np.zeros(6)]+[.03,1.8,-.04]

    def test_contact_edge_stays_literal_and_calibre_is_preserved(self):
        p=self.profile();q,r,c,e=pose_about_contact_edge(p,[0,1],.2,maximum_angle_radians=.3)
        np.testing.assert_array_equal(p[[0,1]],q[[0,1]])
        np.testing.assert_allclose(np.linalg.norm(p[:,None]-p[None],axis=2),np.linalg.norm(q[:,None]-q[None],axis=2),atol=3e-16,rtol=0)
        np.testing.assert_allclose(r.T@r,np.eye(3),atol=1e-15,rtol=0)
        self.assertAlmostEqual(np.linalg.det(r),1)
        self.assertFalse(e['calibreChanged'])

    def test_three_noncollinear_contact_points_prohibit_rotation(self):
        with self.assertRaises(ValueError):contact_edge_pivot(self.profile(),[0,1,2])

    def test_literal_aliases_do_not_lock_a_line(self):
        p=self.profile();p=np.r_[p,p[[0]]]
        q,_,_,_=pose_about_contact_edge(p,[0,1,6],-.15,maximum_angle_radians=.2)
        np.testing.assert_array_equal(q[0],q[6])

    def test_out_of_bound_angle_is_not_clipped(self):
        with self.assertRaises(ValueError):pose_about_contact_edge(self.profile(),[0,1],.3,maximum_angle_radians=.2)
