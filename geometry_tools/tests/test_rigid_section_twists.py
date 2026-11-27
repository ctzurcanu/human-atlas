import unittest
import numpy as np
from scipy.spatial.distance import pdist
from geometry_tools.algorithms.rigid_section_twists import section_twist_lift,replay_section_twists


class SectionTwists(unittest.TestCase):
    def test_rigid_calibre_shared_gap_and_linearization(self):
        ring=np.c_[np.cos(np.arange(6)*np.pi/3),np.sin(np.arange(6)*np.pi/3),np.zeros(6)]*.001
        p=np.r_[ring,ring+[0,0,.00005]];g=np.zeros(len(p),int)
        L,c,l,proof=section_twist_lift(p,g,group_count=1,parameter_radius_m=.001,maximum_angle_radians=.15)
        x=np.array([[.0002,-.0001,.0001],[.0003,.0002,-.0001]])
        q,check=replay_section_twists(p,g,c,l,x)
        np.testing.assert_allclose(pdist(p),pdist(q),atol=1e-18,rtol=2e-15)
        eps=1e-5
        near,_=replay_section_twists(p,g,c,l,x*eps)
        other,_=replay_section_twists(p,g,c,l,-x*eps)
        np.testing.assert_allclose((near-other)/(2*eps),(L@x.ravel()).reshape(-1,3),atol=5e-12,rtol=0)
        self.assertLessEqual(check['maximumActualMotion_m'],proof['maximumExactControlMotionBound_m'])

    def test_fixed_controls_remain_literal(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])*.001
        L,c,l,_=section_twist_lift(p[:2],[0,0],group_count=1,parameter_radius_m=.001,maximum_angle_radians=.1)
        q,_=replay_section_twists(p,[0,0,-1],c,l,np.array([[.0002,0,0],[0,0,.0002]]))
        self.assertTrue(np.array_equal(q[2],p[2]))


if __name__=='__main__':unittest.main()
