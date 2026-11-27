import unittest
import numpy as np
from geometry_tools.algorithms.local_bone_cohort import local_bone_cohort


class LocalCohortTests(unittest.TestCase):
    def setUp(self):
        self.axis=np.array([[0.,0.,0.],[0.,0.,.03],[0.,0.,.06],[0.,0.,.08]])
        self.other=[self.axis+[.03,0,0],self.axis+[-.03,0,0]]
        self.ring=np.array([[.002*np.cos(t),.002*np.sin(t),.05]for t in np.linspace(0,2*np.pi,6,endpoint=False)])

    def test_identity_mapping_is_literal(self):
        q,R,proof=local_bone_cohort(self.ring,[np.arange(6)],self.axis,self.axis,[0,1,0],[0,1,0],self.other)
        np.testing.assert_allclose(q,self.ring,atol=1e-14)
        np.testing.assert_allclose(R,np.repeat(np.eye(3)[None],6,axis=0),atol=1e-14)

    def test_complete_section_translates_without_calibre_loss(self):
        q,R,proof=local_bone_cohort(self.ring,[np.arange(6)],self.axis,self.axis+[.005,0,0],[0,1,0],[0,1,0],self.other)
        np.testing.assert_allclose(q-self.ring,np.tile([.005,0,0],(6,1)),atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None],axis=2),np.linalg.norm(self.ring[:,None]-self.ring[None],axis=2),atol=1e-14)

    def test_neighbouring_digit_and_palm_root_stay_literal(self):
        points=np.r_[self.ring,self.ring+[.03,0,0],[[0,0,-.03]]]
        q,R,proof=local_bone_cohort(points,[np.arange(6),np.arange(6,12)],self.axis,self.axis+[.005,0,0],[0,1,0],[0,1,0],self.other)
        np.testing.assert_array_equal(q[6:],points[6:])

    def test_radial_placement_moves_centres_not_vessel_gauge(self):
        points=self.ring+[.005,0,0]
        q,_,proof=local_bone_cohort(points,[np.arange(6)],self.axis,self.axis,[0,1,0],[0,1,0],self.other,
            radial_scale=lambda t,d:np.full(len(t),1.5))
        np.testing.assert_allclose(q.mean(0)-points.mean(0),[.0025,0,0],atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None],axis=2),np.linalg.norm(points[:,None]-points[None],axis=2),atol=1e-14)

    def test_installed_palm_and_other_digit_are_not_reset_to_donor(self):
        p=np.r_[self.ring,self.ring+[.03,0,0],[[0,0,-.03]]]
        baseline=p+[.002,.001,0.]
        q,_,proof=local_bone_cohort(p,[np.arange(6),np.arange(6,12)],self.axis,self.axis+[.005,0,0],
            [0,1,0],[0,1,0],self.other,baseline_vertices=baseline)
        np.testing.assert_array_equal(q[6:],baseline[6:])
        np.testing.assert_allclose(q[:6]-p[:6],np.tile([.005,0,0],(6,1)),atol=1e-14)
        self.assertTrue(proof['installedBaselineOutsideSupportPreserved'])

    def test_radial_map_is_not_queried_outside_its_live_pose_domain(self):
        p=np.r_[self.ring,[[0,0,-.03]]]
        def field(station,delta):
            self.assertTrue(np.all(station>0))
            return np.ones(len(station))
        q,_,_=local_bone_cohort(p,[np.arange(6)],self.axis,self.axis,[0,1,0],[0,1,0],self.other,radial_scale=field)
        np.testing.assert_allclose(q,p,atol=1e-14)


if __name__=='__main__':unittest.main()
