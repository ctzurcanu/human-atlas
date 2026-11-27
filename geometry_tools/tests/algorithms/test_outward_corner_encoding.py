import unittest
import numpy as np
from geometry_tools.algorithms.surface_normals import quantize_outward_corner_words


class OutwardCornerEncodingTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        self.f=np.array([[0,1,2]])

    def test_bounded_projection_retains_tangent_direction_and_geometry(self):
        original=self.p.copy()
        normals=np.array([[[1.,0.,-.03],[0.,1.,.6],[0.,0.,1.]]])
        words,proof=quantize_outward_corner_words(self.p,self.f,normals,
                                                maximum_adjustment_degrees=2.)
        self.assertTrue(np.all(words[:,:,2]>0))
        self.assertEqual(proof['adjusted_corners'],1)
        self.assertLess(proof['maximum_adjustment_degrees'],2.)
        self.assertEqual(words[0,0,1],0)
        np.testing.assert_array_equal(self.p,original)

    def test_exceeded_bound_is_rejected(self):
        normals=np.tile([1.,0.,-.2],(1,3,1))
        with self.assertRaises(ValueError):
            quantize_outward_corner_words(self.p,self.f,normals,maximum_adjustment_degrees=1.)

    def test_reciprocal_winding_and_field_encode_exact_opposites(self):
        normals=np.array([[[1.,0.,-.03],[0.,1.,.6],[0.,0.,1.]]])
        a,_=quantize_outward_corner_words(self.p,self.f,normals,maximum_adjustment_degrees=2.)
        b,_=quantize_outward_corner_words(self.p,self.f[:,[0,2,1]],-normals[:,[0,2,1]],
                                        maximum_adjustment_degrees=2.)
        np.testing.assert_array_equal(a,-b[:,[0,2,1]])

if __name__=='__main__':unittest.main()
