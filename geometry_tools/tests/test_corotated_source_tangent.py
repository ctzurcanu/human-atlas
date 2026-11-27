import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.corotated_source_tangent import corotated_tangent_rows

class CorotatedSourceTest(unittest.TestCase):
    def setUp(self):
        self.v=np.array([[0.,0.,0.],[.001,0.,0.],[0.,.001,0.]])
        self.a=Mesh(self.v,np.array([[0,1,2]]),'test')
        angle=.13;self.rotation=np.array([[1.,0.,0.],[0.,np.cos(angle),-np.sin(angle)],[0.,np.sin(angle),np.cos(angle)]])

    def test_proper_bending_does_not_count_as_stretch(self):
        q=self.v@self.rotation.T;b=Mesh(q,self.a.faces,'test')
        C,offset,p=corotated_tangent_rows(self.a,b,np.arange(3),3)
        residual=(C@(q-self.v).ravel()).reshape(-1,6)+offset
        self.assertLess(np.linalg.norm(residual),1e-12)
        self.assertLess(p['maximumReviewedOriginalSingularStrain'],1e-12)
        self.assertTrue(p['properUnscaledBendFrames'])

    def test_stretch_remains_relative_to_original_after_changing_frame(self):
        q=(self.v@np.diag([1.04,1.,1.]))@self.rotation.T;b=Mesh(q,self.a.faces,'test')
        C,offset,_=corotated_tangent_rows(self.a,b,np.arange(3),3)
        residual=(C@(q-self.v).ravel()).reshape(-1,6)+offset
        self.assertAlmostEqual(np.linalg.norm(residual),.04,places=12)
        # A second 4% stretch is 8.16% cumulative, not a fresh 4% allowance.
        twice=(self.v@np.diag([1.04**2,1.,1.]))@self.rotation.T
        self.assertGreater(np.linalg.norm((C@(twice-self.v).ravel()).reshape(-1,6)+offset),.08)

if __name__=='__main__':unittest.main()
