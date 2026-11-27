import unittest
import numpy as np
from geometry_tools.algorithms.bounded_float32_rounding import round_diagnosed_crossings
from geometry_tools.distance import triangle_pairs


class CutRoundingBounds(unittest.TestCase):
    def fixture(self):
        p=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[.1,.900000001,0]],np.float32).astype(float)
        f=np.array([[0,1,2],[1,3,4]])
        self.assertEqual(triangle_pairs(p[f[:1]],p[f[1:]])['relation'][0],'coplanar-overlap')
        return p,f

    def test_local_cap_cannot_be_overridden_by_global_bound(self):
        p,f=self.fixture();limits=np.full(len(p),1e-7);limits[4]=0.
        r=round_diagnosed_crossings(p,f,[[0,1]],origin=[0,0,0],protected=range(4),
            maximum_displacement_m=1e-7,maximum_displacement_by_vertex=limits,seconds=2.)
        self.assertFalse(r['evidence']['passed']);np.testing.assert_array_equal(r['positions'],p)

    def test_one_mutable_point_preserves_all_other_controls(self):
        p,f=self.fixture()
        r=round_diagnosed_crossings(p,f,[[0,1]],origin=[0,0,0],protected=range(4),
            maximum_displacement_m=1e-7,maximum_displacement_by_vertex=np.full(len(p),1e-7),seconds=2.)
        self.assertTrue(r['evidence']['passed']);np.testing.assert_array_equal(r['positions'][:4],p[:4])
        self.assertLessEqual(np.linalg.norm(r['positions'][4]-p[4]),1e-7)


if __name__=='__main__':unittest.main()
