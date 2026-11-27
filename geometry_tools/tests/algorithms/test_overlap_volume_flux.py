import unittest
import numpy as np
import manifold3d as md
from geometry_tools.model import Mesh
from geometry_tools.algorithms.occupancy_classifier import close_literal_rims
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.algorithms.overlap_volume_flux import overlap_volume_flux


class OverlapFluxTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
        f=np.array([[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        m,_=close_literal_rims(Mesh(v,f,'test:consistent-length'))
        fixed,_=make_solid(m.vertices+[.4,.1,.2],m.faces)
        return m,fixed

    def test_rigid_and_nonuniform_source_gradients_match_finite_differences(self):
        m,fixed=self.fixture();operator=np.c_[np.ones(len(m.vertices)),m.vertices[:,0]];oid=md.Manifold.reserve_ids(1)
        volume,gradient,proof=overlap_volume_flux(m.vertices,m.faces,operator,fixed,original_id=oid)
        self.assertAlmostEqual(volume,.6*.9*.8);self.assertTrue(proof['originalFaceProvenanceVerified'])
        for control in range(2):
            for dim,axis in enumerate([0,2]):
                delta=np.zeros_like(m.vertices);delta[:,axis]=operator[:,control]*1e-6
                values=[overlap_volume_flux(m.vertices+s*delta,m.faces,operator,fixed,original_id=oid)[0]for s in [-1,1]]
                self.assertAlmostEqual(gradient[control,dim],(values[1]-values[0])/2e-6,places=6)

    def test_disjoint_shapes_have_zero_overlap_and_gradient(self):
        m,fixed=self.fixture();volume,gradient,_=overlap_volume_flux(m.vertices+[10,0,0],m.faces,np.ones((len(m.vertices),1)),fixed,original_id=md.Manifold.reserve_ids(1))
        self.assertEqual(volume,0);np.testing.assert_array_equal(gradient,0)
