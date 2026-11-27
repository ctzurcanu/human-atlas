import unittest
import numpy as np
from geometry_tools.algorithms.surface_deformation_guards import surface_deformation_guards,evaluate_surface_guards
from geometry_tools.algorithms.bounded_cell_fairing import fair_rigid_cells


class SourceSurfaceGradientTest(unittest.TestCase):
    def test_rigid_translation_has_zero_gradient_even_on_thin_triangle(self):
        v=np.array([[0.,0.,0.],[.01,0.,0.],[.003,.000003,0.]])
        g=surface_deformation_guards(v,[[0,1,2]],np.arange(3),range(3));d=np.tile([.01,0.,.005],(3,1));r=evaluate_surface_guards(d,g)
        self.assertLess(r['maximumNorm'],1e-10);self.assertEqual(r['violatingFaces'],0)

    def test_micron_face_stretch_is_detected_from_actual_surface_not_cell_centres(self):
        v=np.array([[0.,0.,0.],[.01,0.,0.],[.003,.000003,0.]])
        g=surface_deformation_guards(v,[[0,1,2]],np.arange(3),range(3));d=np.zeros((3,3));d[2,2]=.0002;r=evaluate_surface_guards(d,g)
        self.assertAlmostEqual(r['maximumNorm'],.0002/.000003);self.assertEqual(r['violatingFaces'],1)

    def test_repeated_rigid_profile_cells_cancel_exactly(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        g=surface_deformation_guards(v,[[0,1,2]],[0,0,1],[0,1]);d=np.array([[.1,0.,0.],[.1,0.,.2]])
        self.assertAlmostEqual(evaluate_surface_guards(d,g)['maximumNorm'],.2)

    def test_fitter_prevents_the_actual_thin_face_stretch_without_moving_its_fixed_edge(self):
        v=np.array([[0.,0.,0.],[.01,0.,0.],[.003,.000003,0.]])
        g=surface_deformation_guards(v,[[0,1,2]],np.arange(3),[2]);d,p=fair_rigid_cells(v,[[0,2],[1,2]],[2],{0,1},{2:dict(translation_m=[0.,.0002],weight=100)},surface_guards=g)
        self.assertTrue(p['feasible']);self.assertTrue(np.array_equal(d[:2],np.zeros((2,3))))
        self.assertLessEqual(evaluate_surface_guards(d,g)['maximumNorm'],.9+1e-8);self.assertLessEqual(abs(d[2,2]),.9*.000003+1e-10)


if __name__=='__main__':unittest.main()
