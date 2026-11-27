import unittest
import numpy as np
from geometry_tools.algorithms.triangle_translation_clearance import separating_translation

class TriangleTranslation(unittest.TestCase):
    def setUp(self):
        self.a=np.array([[-1.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        self.b=np.array([[0.,.1,-1.],[0.,.1,1.],[0.,.8,0.]])

    def test_every_corner_has_the_declared_support_gap(self):
        d=np.array([1.,.3,.2]);d/=np.linalg.norm(d)
        motion,p=separating_translation(self.a,self.b,d,clearance_m=.02)
        moved=self.a+(motion+np.array([0.,.3,.7]))[:,None]*d
        n=np.array(p['separatingAxis']);gap=(moved@n).min()-(self.b@n).max()
        self.assertGreaterEqual(gap,.02-1e-12)
        self.assertGreater(motion,0)

    def test_rigid_world_frame_change_preserves_the_certificate(self):
        d=np.array([0.,1.,0.]);t,p=separating_translation(self.a,self.b,d,clearance_m=.01)
        rotation=np.array([[0.,-1,0],[0,0,-1],[1,0,0]]);origin=np.array([3.,-7.,12.])
        u,q=separating_translation(self.a@rotation.T+origin,self.b@rotation.T+origin,rotation@d,clearance_m=.01)
        self.assertAlmostEqual(t,u,places=12)
        self.assertAlmostEqual(p['certifiedProjectedGap_m'],q['certifiedProjectedGap_m'],places=12)

    def test_coplanar_support_and_invalid_inputs(self):
        t,p=separating_translation(self.a,self.a,[1,0,0],clearance_m=.01)
        self.assertAlmostEqual(t,2.01)
        self.assertGreaterEqual(p['certifiedProjectedGap_m'],.01-1e-12)
        for d in [[0,0,0],[np.nan,0,0]]:
            with self.assertRaises(ValueError):separating_translation(self.a,self.b,d)
        with self.assertRaises(ValueError):separating_translation(np.zeros((3,3)),self.b,[1,0,0])

if __name__=='__main__':unittest.main()
