import unittest
import numpy as np
from geometry_tools.distance import exact_intersection_segment


class ExactIntersectionSegmentTests(unittest.TestCase):
    def test_complete_interval_not_a_single_closest_witness(self):
        a=np.array([[0,0,0],[2,0,0],[0,2,0]],float)
        b=np.array([[.5,-1,-1],[.5,2,-1],[.5,.5,1]],float)
        r=exact_intersection_segment(a,b)
        self.assertEqual(r['relation'],'crossing')
        np.testing.assert_array_equal(r['endpoints'],[[.5,0,0],[.5,1.25,0]])
        np.testing.assert_array_equal(r['coordinate_error_upper_m'],0)

    def test_translation_and_input_order_keep_the_same_segment(self):
        a=np.array([[0,0,0],[2,0,0],[0,2,0]],float)
        b=np.array([[.5,-1,-1],[.5,2,-1],[.5,.5,1]],float)
        shift=np.array([1024.,-2048.,4096.])
        for x,y in [(a,b),(b,a),(a[::-1],b[::-1])]:
            r=exact_intersection_segment(x+shift,y+shift)
            np.testing.assert_array_equal(r['endpoints']-shift,[[.5,0,0],[.5,1.25,0]])

    def test_non_binary_endpoint_has_an_explicit_rounding_bound(self):
        a=np.array([[0,0,0],[2,0,0],[0,2,0]],float)
        b=np.array([[0,0,-1],[1,0,2],[0,1,2]],float)
        r=exact_intersection_segment(a,b)
        self.assertEqual(r['relation'],'crossing')
        np.testing.assert_allclose(r['endpoints'],[[0,1/3,0],[1/3,0,0]],rtol=0,atol=1e-16)
        self.assertTrue(np.any(r['coordinate_error_upper_m']>0))
        self.assertLess(r['coordinate_error_upper_m'].max(),1e-16)

    def test_shared_edge_is_not_an_interior_crossing(self):
        a=[[0,0,0],[1,0,0],[0,1,0]];b=[[0,0,0],[1,0,0],[0,0,1]]
        r=exact_intersection_segment(a,b)
        self.assertEqual(r['relation'],'touching')
        np.testing.assert_array_equal(r['endpoints'],[[0,0,0],[1,0,0]])

    def test_coplanar_and_parallel_disjoint_are_not_invented_segments(self):
        a=np.array([[0,0,0],[1,0,0],[0,1,0]],float)
        for delta,relation in [(0,'coplanar'),(1,'disjoint')]:
            r=exact_intersection_segment(a,a+[0,0,delta])
            self.assertEqual(r['relation'],relation);self.assertEqual(r['endpoints'].shape,(0,3))


if __name__=='__main__':unittest.main()
