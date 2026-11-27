import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.algorithms.corotating_pair_intervals import corotating_vertex_pair_rows

class CorotatingPairIntervals(unittest.TestCase):
    def test_common_proper_pose_does_not_invent_interval_loss(self):
        a=np.array([[.001,.002,.003]]);b=np.array([[.002,.004,.001]])
        axis=np.array([[1.,0.,0.]])
        A,gap=corotating_vertex_pair_rows(a,b,[0],[0],[[0.,0.,0.]],[.02],axis,[[0]],[[1.]])
        np.testing.assert_allclose(A.toarray(),0,atol=1e-15)
        np.testing.assert_allclose(gap,[-.001])
    def test_local_row_matches_independent_finite_rotation_and_plane(self):
        a=np.array([[.001,.002,.003]]);b=np.array([[.002,.004,.001]])
        c=np.array([[0.,0.,0.],[.002,.003,0.]]);l=np.array([.02,.03]);n=np.array([[1.,0.,0.]])
        A,gap=corotating_vertex_pair_rows(a,b,[0],[1],c,l,n,[[0,1]],[[.5,.5]])
        p=np.array([[.0002,0.,.0003],[0.,-.0001,.0002],[.0001,.0002,0.],[.0001,0.,.0003]])
        eps=1e-5;omega=p[2:]/l[:,None]
        qa=c[0]+Rotation.from_rotvec(eps*omega[0]).apply(a[0]-c[0])+eps*p[0]
        qb=c[1]+Rotation.from_rotvec(eps*omega[1]).apply(b[0]-c[1])+eps*p[1]
        normal=Rotation.from_rotvec(eps*omega.mean(0)).apply(n[0])
        exact=((qa-qb)@normal-gap[0])/eps
        self.assertAlmostEqual(float((A@p.ravel())[0]),float(exact),places=10)
