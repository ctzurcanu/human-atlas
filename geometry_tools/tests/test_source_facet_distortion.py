import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_facet_distortion import facet_tangent_rows

class SourceFacetDistortionTests(unittest.TestCase):
    def test_rigid_translation_has_zero_tangent_change(self):
        m=Mesh(np.array([[0,0,0],[2,0,0],[0,1,0]],float),np.array([[0,1,2]]),'fixture:metres')
        C,p=facet_tangent_rows(m,np.zeros(3,int),1)
        np.testing.assert_allclose(C@np.array([.1,.2,.3]),0,atol=1e-15)

    def test_affine_change_matches_literal_tangent_frame(self):
        m=Mesh(np.array([[0,0,0],[2,0,0],[0,1,0]],float),np.array([[0,1,2]]),'fixture:metres')
        C,p=facet_tangent_rows(m,np.array([0,1,2]),3)
        d=np.array([[.01,.02,.03],[.04,-.01,.1],[.1,.05,.08]])
        E=np.c_[m.vertices[1]-m.vertices[0],m.vertices[2]-m.vertices[0]]
        u,s,v=np.linalg.svd(E,full_matrices=False)
        R=v.T/s[None,:]
        expected=(np.c_[d[1]-d[0],d[2]-d[0]]@R).T.ravel()
        np.testing.assert_allclose(C@d.ravel(),expected,atol=1e-14)

    def test_degenerate_source_is_rejected(self):
        m=Mesh(np.array([[0,0,0],[1,0,0],[2,0,0]],float),np.array([[0,1,2]]),'fixture:metres')
        with self.assertRaises(AssertionError):facet_tangent_rows(m,np.array([0,1,2]),3)

if __name__=='__main__':unittest.main()
