import unittest
import numpy as np
from geometry_tools.algorithms.coupled_surface_quality import improve_triangles


class CoupledSurfaceQualityTests(unittest.TestCase):
    def fixture(self):
        return np.array([[0.,0.,0.],[1.,0.,0.],[.5,1e-8,0.],[.5,-1.,0.]]),np.array([[0,1,2],[1,0,3]])

    def test_removes_a_sliver_without_moving_any_control_or_boundary(self):
        p,f=self.fixture();before=p.copy();result,evidence=improve_triangles(p,p,f)
        np.testing.assert_array_equal(p,before)
        self.assertEqual(len(evidence['flips']),1)
        self.assertTrue(evidence['patch_boundary_edges_literal'])
        t=p[result];self.assertTrue(np.all(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2]>0))
        self.assertGreater(min(x['minimum_quality_after']for x in evidence['flips']),.4)

    def test_preserves_a_nonplanar_patch_above_its_deviation_budget(self):
        p,f=self.fixture();p[2,1]=.02;p[3,2]=.1
        result,evidence=improve_triangles(p,p,f,plane_tolerance_m=1e-6,poor_quality=.5)
        np.testing.assert_array_equal(result,f)
        self.assertEqual(evidence['flips'],[])

    def test_opposed_tissue_owner_gets_the_same_reversed_patch(self):
        p,f=self.fixture();a,_=improve_triangles(p,p,f);b,_=improve_triangles(p,p,f[:,[0,2,1]])
        self.assertEqual({tuple(sorted(t))for t in a},{tuple(sorted(t))for t in b})
        t=p[b];self.assertTrue(np.all(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2]<0))

    def test_preserves_a_declared_attribute_seam(self):
        p,f=self.fixture()
        result,evidence=improve_triangles(p,p,f,protected_edges=[(0,1)])
        np.testing.assert_array_equal(result,f)
        self.assertEqual(evidence['flips'],[])
        self.assertTrue(evidence['protected_edges_literal'])
        self.assertEqual(evidence['protected_feature_edges'],1)

    def test_rejects_an_unknown_protected_edge(self):
        p,f=self.fixture()
        with self.assertRaisesRegex(ValueError,'Protected edges'):
            improve_triangles(p,p,f,protected_edges=[(2,3)])


if __name__=='__main__':unittest.main()
