import unittest
import numpy as np
from geometry_tools.algorithms.coupled_triangle_separation import separating_pair_proposals


class CoupledTriangleSeparationTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[-1.,0,0],[1,0,0],[0,1,0],[0,-.2,-1],[0,-.2,1],[0,.8,0]])

    def test_all_emitted_support_certificates_and_fixed_corner(self):
        rows=separating_pair_proposals(self.p,[0,1,2],[3,4,5],locked=[2],gap_m=.01)
        self.assertTrue(rows)
        for r in rows:
            q=self.p.copy();q[r['ids']]=r['positions']
            np.testing.assert_array_equal(q[2],self.p[2])
            self.assertGreaterEqual((q[:3]@r['axis']).min()-(q[3:]@r['axis']).max(),.01-1e-12)

    def test_rigid_frame_and_scale_covariance(self):
        first=separating_pair_proposals(self.p,[0,1,2],[3,4,5],gap_m=.01)[0]
        R=np.array([[0,1.,0],[0,0,1],[1,0,0]])
        other=separating_pair_proposals(self.p@R.T*.003+[3,4,5],[0,1,2],[3,4,5],gap_m=.00003)[0]
        self.assertAlmostEqual(other['squaredMotion'],first['squaredMotion']*.003**2,places=14)

    def test_protected_crossing_and_shared_corner_not_released(self):
        self.assertEqual(separating_pair_proposals(self.p,[0,1,2],[3,4,5],locked=range(6)),[])
        with self.assertRaises(ValueError):separating_pair_proposals(self.p,[0,1,2],[2,4,5])
