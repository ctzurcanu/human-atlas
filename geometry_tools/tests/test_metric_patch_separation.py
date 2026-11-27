import unittest
import numpy as np
from geometry_tools.algorithms.metric_patch_separation import separate_patch
from geometry_tools.algorithms.hinge_isometry import intrinsic_facet_strains
from geometry_tools.algorithms.exact_triangle_separation import exact_triangle_separation


class MetricPatchTests(unittest.TestCase):
    def source(self):
        # Two parallel, overlapping triangles; only the first may move.
        p=np.array([[0,0,0],[1e-3,0,0],[0,1e-3,0],
                    [0,0,0],[1e-3,0,0],[0,1e-3,0]],float)
        return p,np.array([[0,1,2],[3,4,5]])

    def test_rigid_translation_and_literal_fixed_triangle(self):
        p,f=self.source()
        q,proof=separate_patch(p,f,f[0],f[1],[0,1,2],np.array([0.,0.,1.]),
            gap_m=1e-6,maximum_motion_m=2e-6,maximum_strain=.075,wall_seconds=1)
        np.testing.assert_array_equal(q[3:],p[3:])
        np.testing.assert_allclose(q[:3]-p[:3],np.tile([0,0,1e-6],(3,1)),atol=1e-12)
        self.assertLess(np.abs(intrinsic_facet_strains(p,q,f)).max(),1e-9)
        self.assertTrue(exact_triangle_separation(q[f[0]],q[f[1]])['disjoint'])
        self.assertTrue(proof['success'])

    def test_impossible_fixed_corner_is_rejected(self):
        p,f=self.source()
        with self.assertRaises(ValueError):
            separate_patch(p,f,f[0],f[1],[1,2],np.array([0.,0.,1.]),
                gap_m=1e-6,maximum_motion_m=2e-6,maximum_strain=.075,wall_seconds=1)

    def test_motion_domain_cannot_be_extended(self):
        p,f=self.source()
        with self.assertRaises(ValueError):
            separate_patch(p,f,f[0],f[1],[0,1,2],np.array([0.,0.,1.]),
                gap_m=2e-6,maximum_motion_m=1e-6,maximum_strain=.075,wall_seconds=1)


if __name__=='__main__':unittest.main()
