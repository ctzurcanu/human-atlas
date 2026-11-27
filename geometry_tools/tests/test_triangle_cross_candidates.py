import unittest
import numpy as np
from geometry_tools.algorithms.triangle_cross_candidates import TriangleCrossIndex


class ConservativeCrossCandidates(unittest.TestCase):
    def test_matches_every_literal_box_overlap_across_scales(self):
        rng=np.random.default_rng(431)
        a=rng.normal(size=(91,3,3));b=rng.normal(size=(73,3,3))
        b[:20]*=1e-5;b[20:25]=0;b[25:30]*=100
        expected=set(map(tuple,np.argwhere(np.all(a.max(1)[:,None]>=b.min(1)[None],axis=2)&np.all(b.max(1)[None]>=a.min(1)[:,None],axis=2))))
        actual=set(map(tuple,np.concatenate(list(TriangleCrossIndex(b).batches(a,batch_size=13)))))
        self.assertEqual(expected,actual)

    def test_budget_fails_instead_of_returning_partial_success(self):
        t=np.zeros((4,3,3))
        with self.assertRaises(RuntimeError):list(TriangleCrossIndex(t).batches(t,maximum_proposals=2))


if __name__=='__main__':unittest.main()
