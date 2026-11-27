import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.algorithms.common_pose_lift import common_pose_lift

class CommonPoseLift(unittest.TestCase):
    def test_derivative_matches_single_proper_motion_with_source_offsets(self):
        old=np.array([[.001,0,0],[-.001,.002,0]])
        new=np.array([[0.,0.,0.]])
        a=np.array([.006,.007]);b=np.array([.008]);mapping=np.array([0,0])
        T=common_pose_lift(old,a,new,b,mapping)
        p=np.array([[.0003,-.0001,.0002],[.0002,-.0003,.0001]])
        linear=(T@p.ravel()).reshape(4,3)
        epsilon=1e-5;R=Rotation.from_rotvec(p[1]/b[0]*epsilon).as_matrix()
        exact=p[0]*epsilon+(old-new)@(R-np.eye(3)).T
        np.testing.assert_allclose(linear[:2],exact/epsilon,atol=2e-11,rtol=0)
        np.testing.assert_allclose(linear[2:]/a[:,None],np.broadcast_to(p[1]/b[0],(2,3)))
    def test_identity_cells_leave_all_parameters_literal(self):
        c=np.array([[.1,.2,.3],[.2,.1,.4]]);a=np.array([.1,.2])
        np.testing.assert_array_equal(common_pose_lift(c,a,c,a,np.arange(2)).toarray(),np.eye(12))
