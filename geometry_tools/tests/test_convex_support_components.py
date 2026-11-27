import unittest
import numpy as np
from geometry_tools.algorithms.convex_support_components import separating_components


class ConvexSupportComponentsTest(unittest.TestCase):
    def test_disjoint_concave_profiles_with_overlapping_convex_hulls(self):
        a=np.array([[0,0],[4,0],[4,4],[3,4],[3,1],[1,1],[1,4],[0,4]])*.001
        b=np.array([[1.5,2],[2.5,2],[2.5,3],[1.5,3]])*.001
        pieces=separating_components(a,b,clearance_m=.0001)
        self.assertGreater(len(pieces),1)
        for i,j,n,low,high in pieces:
            self.assertGreaterEqual(float((a[i]@n).min()-(b[j]@n).max()),.0001-1e-12)
            self.assertTrue(set(i)<=set(range(len(a))));self.assertTrue(set(j)<=set(range(len(b))))

    def test_actual_overlap_is_rejected(self):
        a=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        with self.assertRaises(ValueError):separating_components(a,a+.1,clearance_m=.001)

    def test_open_segment_remains_a_constraint_segment(self):
        a=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        r=separating_components(a,[[2.,0.],[2.,1.]],clearance_m=.1)
        self.assertEqual(len(r),1);self.assertEqual(len(r[0][1]),2)


if __name__=='__main__':unittest.main()
