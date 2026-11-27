import unittest
import numpy as np
from geometry_tools.algorithms.profile_support_separation import minimum_tangent_support_axis


class LocalSupportAxisTests(unittest.TestCase):
    def test_complete_cells_choose_short_separation_axis(self):
        a=np.array([[-2.,-.1,0.],[2.,-.1,0.],[2.,.1,0.],[-2.,.1,0.]])
        b=a+np.array([0.,.15,0.])
        u,p=minimum_tangent_support_axis(a,b,[0,0,1],[1,0,0])
        self.assertAlmostEqual(u[2],0)
        self.assertAlmostEqual(np.linalg.norm(u),1)
        self.assertAlmostEqual(p['requiredTranslation_m'],.05)
        gap=(a@u).min()-(b@u).max()
        self.assertAlmostEqual(gap,-p['requiredTranslation_m'])
        self.assertFalse(p['anatomicalAcceptance'])

    def test_rotated_plane_retains_all_points(self):
        a=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        b=a+np.array([0.,.1,.2])
        n=np.array([1.,2.,3.]);n/=np.linalg.norm(n)
        u,p=minimum_tangent_support_axis(a,b,n,[0,1,0])
        self.assertLess(abs(u@n),1e-14)
        self.assertAlmostEqual((b@u).max()-(a@u).min(),p['requiredTranslation_m'])
