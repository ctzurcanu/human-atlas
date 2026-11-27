import unittest
import numpy as np
from geometry_tools.algorithms.layer_tangent_separation import layer_cross_course_axes


class LayerTangentSeparation(unittest.TestCase):
    def test_nearly_parallel_courses_separate_in_layer_not_in_depth(self):
        normal=np.array([[0.,0.,1.]])
        a=np.array([[0.,1.,0.]]);b=np.array([[.05,1.,0.]])
        old=np.cross(a,b);old/=np.linalg.norm(old,axis=1)[:,None]
        self.assertAlmostEqual(abs((old*normal).sum()),1.)
        axis,proof=layer_cross_course_axes(normal,a,b,side_reference=[1.,0.,0.])
        self.assertAlmostEqual((axis*normal).sum(),0.,places=14)
        self.assertGreater(axis[0,0],0.)
        self.assertLess(proof['maximumCourseComponent'],1e-14)
        flipped,_=layer_cross_course_axes(normal,-a,-b,side_reference=[1.,0.,0.])
        np.testing.assert_allclose(axis,flipped,atol=1e-14)

    def test_perpendicular_entry_and_ambiguous_side_are_rejected(self):
        with self.assertRaises(ValueError):
            layer_cross_course_axes([[0.,0.,1.]],[[0.,0.,1.]],[[0.,0.,1.]],side_reference=[1.,0.,0.])
        with self.assertRaises(ValueError):
            layer_cross_course_axes([[0.,0.,1.]],[[0.,1.,0.]],[[0.,1.,0.]],side_reference=[0.,0.,1.])


if __name__=='__main__':unittest.main()
