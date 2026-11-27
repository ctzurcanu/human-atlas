import unittest
import numpy as np
from geometry_tools.algorithms.course_clearance_field import cubic_carrier


class CourseClearanceFieldTests(unittest.TestCase):
    def test_aliases_overlapping_sections_and_branch_anchor_share_translation(self):
        axis=np.array([[0,0,0],[0,0,1.]])
        p=np.array([[.1,0,.1],[-.1,0,.1],[0,.1,.2],[.1,0,.1],[.2,0,.3],[0,0,.9]])
        w,s,_,proof=cubic_carrier(p,axis,[[0,1],[1,2,4]],anchors={1:np.array([0,0,.15])})
        np.testing.assert_array_equal(w[0],w[3])
        np.testing.assert_array_equal(w[[0,1,2,4]],np.repeat(w[:1],4,axis=0))
        np.testing.assert_allclose(s[[0,1,2,3,4]],.15)
        self.assertTrue(proof['exactSourceAliasesCoMove'])

    def test_control_bound_bounds_every_point(self):
        axis=np.array([[0,0,0],[0,0,1.]])
        p=np.c_[np.zeros(101),np.zeros(101),np.linspace(0,1,101)]
        w,_,metric,_=cubic_carrier(p,axis,[])
        rng=np.random.default_rng(7);d=rng.normal(size=(12,3));d*=.002/np.linalg.norm(d,axis=1)[:,None]
        self.assertLessEqual(np.linalg.norm(w@d,axis=1).max(),.002+1e-15)
        np.linalg.cholesky(metric)

    def test_does_not_merge_nearby_distinct_controls(self):
        axis=np.array([[0,0,0],[0,0,1.]])
        p=np.array([[0,0,.1],[0,0,.10000001]])
        w,_,_,_=cubic_carrier(p,axis,[])
        self.assertFalse(np.array_equal(w[0],w[1]))


if __name__=='__main__':unittest.main()
