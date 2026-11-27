import unittest
import numpy as np
from geometry_tools.algorithms.rigid_translation_intervals import triangle_overlap_interval,clear_translation_intervals


class RigidTranslationIntervals(unittest.TestCase):
    def test_transverse_triangle_has_analytic_contact_interval(self):
        a=np.array([[-2.,-.5,-1.],[-2.,.5,-1.],[-2.,0.,1.]])
        b=np.array([[-1.,-1.,0.],[1.,-1.,0.],[0.,1.,0.]])
        row=triangle_overlap_interval(a,b,[1.,0.,0.],maximum_distance_m=4.)
        lo,hi=row['forbidden_interval_m']
        self.assertLessEqual(lo,1.375);self.assertGreaterEqual(hi,2.625)
        np.testing.assert_allclose([lo,hi],[1.375,2.625],atol=1e-10,rtol=0)
        self.assertEqual(row['forbidden_interval_m'],triangle_overlap_interval(a,b,[7.,0.,0.],maximum_distance_m=4.)['forbidden_interval_m'])

    def test_coplanar_slide_uses_in_plane_axes(self):
        a=np.array([[-3.,0.,0.],[-2.,0.,0.],[-3.,1.,0.]])
        b=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        row=triangle_overlap_interval(a,b,[1.,0.,0.],maximum_distance_m=5.)
        np.testing.assert_allclose(row['forbidden_interval_m'],[2.,4.],atol=1e-10,rtol=0)

    def test_separate_plane_and_union(self):
        a=np.array([[0.,0.,1.],[1.,0.,1.],[0.,1.,1.]])
        b=a-[0.,0.,1.]
        self.assertIsNone(triangle_overlap_interval(a,b,[1.,0.,0.],maximum_distance_m=2.)['forbidden_interval_m'])
        row=clear_translation_intervals([[.3,.5],[.1,.4],[.8,1.2],None],maximum_distance_m=1.)
        self.assertEqual(row['blocked_closed_intervals_m'],[[.1,.5],[.8,1.]])
        self.assertEqual(row['clear_open_interiors_m'],[[0.,.1],[.5,.8]])

    def test_degenerate_input_is_not_an_aperture(self):
        a=np.zeros((3,3));b=np.eye(3)
        with self.assertRaises(ValueError):triangle_overlap_interval(a,b,[1.,0.,0.],maximum_distance_m=1.)


if __name__=='__main__':unittest.main()
