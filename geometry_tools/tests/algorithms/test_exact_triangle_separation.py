import unittest
import numpy as np
from geometry_tools.algorithms.exact_triangle_separation import exact_triangle_separation,exact_segment_triangle_intersection


class ExactTriangleSeparationTests(unittest.TestCase):
    def setUp(self):self.a=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    def test_plane_gap_smaller_than_a_generic_distance_tolerance_is_retained(self):
        b=self.a.copy();b[:,2]=2**-50
        r=exact_triangle_separation(self.a,b);self.assertTrue(r['disjoint'])
        self.assertEqual(r['projected_separation_m'],2**-50)
    def test_perpendicular_crossing_and_boundary_touch_are_not_separated(self):
        b=np.array([[.25,.25,-1.],[.25,.25,1.],[.75,.25,0.]])
        self.assertFalse(exact_triangle_separation(self.a,b)['disjoint'])
        b=self.a+[1.,0.,0.];self.assertFalse(exact_triangle_separation(self.a,b)['disjoint'])
    def test_coplanar_triangles_need_in_plane_axes(self):
        b=self.a+[.8,.8,0.];self.assertTrue(exact_triangle_separation(self.a,b)['disjoint'])
        b=self.a+[.1,.1,0.];self.assertFalse(exact_triangle_separation(self.a,b)['disjoint'])
    def test_large_translation_does_not_destroy_a_representable_small_gap(self):
        a=self.a*2**-20+[1000.,-500.,250.];b=a.copy();b[:,2]=np.nextafter(b[:,2],np.inf)
        self.assertTrue(exact_triangle_separation(a,b)['disjoint'])
    def test_degenerate_inputs_are_explicit(self):
        with self.assertRaisesRegex(ValueError,'Degenerate'):
            exact_triangle_separation(self.a,np.zeros((3,3)))

    def test_segment_crossing_and_declared_endpoint_are_distinct(self):
        r=exact_segment_triangle_intersection([.2,.2,-1],[.2,.2,1],self.a)
        self.assertTrue(r['intersects']);self.assertFalse(r['at_first_endpoint']);self.assertFalse(r['at_last_endpoint'])
        r=exact_segment_triangle_intersection([0,0,0],[.2,.2,1],self.a)
        self.assertTrue(r['at_first_endpoint'])
        self.assertFalse(exact_segment_triangle_intersection([2,2,-1],[2,2,1],self.a)['intersects'])

    def test_coplanar_segment_resolves_overlap_and_separation_exactly(self):
        self.assertTrue(exact_segment_triangle_intersection([-.2,.2,0],[.5,.2,0],self.a)['intersects'])
        self.assertFalse(exact_segment_triangle_intersection([1.1,0,0],[2,0,0],self.a)['intersects'])
        self.assertFalse(exact_segment_triangle_intersection([0,0,2**-50],[1,0,2**-50],self.a)['intersects'])


if __name__=='__main__':unittest.main()
