import unittest
import numpy as np
from geometry_tools.algorithms.surface_normals import transport_corner_words


class FacetNormalTransportTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        angle=np.radians(30);r=np.array([[1.,0.,0.],[0.,np.cos(angle),-np.sin(angle)],[0.,np.sin(angle),np.cos(angle)]])
        return v,v@r.T,np.array([[0,1,2]]),np.array([[[0,0,32767]]*3],np.int16),r

    def test_corner_fields_follow_the_proper_facet_rotation(self):
        v,p,f,w,r=self.fixture();out,e=transport_corner_words(v,p,f,w,maximum_facet_rotation_degrees=31)
        expected=np.array([0.,0.,1.])@r.T;actual=out[0,0].astype(float);actual/=np.linalg.norm(actual)
        np.testing.assert_allclose(actual,expected,atol=2e-5);self.assertAlmostEqual(e['maximum_facet_rotation_degrees'],30)

    def test_rotation_bound_is_not_expanded(self):
        v,p,f,w,_=self.fixture()
        with self.assertRaises(ValueError):transport_corner_words(v,p,f,w,maximum_facet_rotation_degrees=29)

    def test_a_reversed_source_field_is_rejected(self):
        v,p,f,w,_=self.fixture()
        with self.assertRaises(ValueError):transport_corner_words(v,p,f,-w,maximum_facet_rotation_degrees=31)


if __name__=='__main__':unittest.main()
