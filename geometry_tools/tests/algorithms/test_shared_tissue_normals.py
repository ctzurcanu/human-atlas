import unittest
import numpy as np
from geometry_tools.algorithms.shared_tissue_normals import shared_tissue_corner_normals


class SharedTissueNormalsTests(unittest.TestCase):
    def fixture(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.]])
        return p,{'upper':np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]]),
                  'lower':np.array([[0,1,2],[0,4,1],[1,4,2],[2,4,0]])}

    def test_shared_cap_stays_sharp_and_reciprocal_while_outer_surface_is_smooth(self):
        p,cells=self.fixture();original=p.copy();n,proof=shared_tissue_corner_normals(p,cells)
        np.testing.assert_array_equal(p,original);self.assertEqual(proof['shared_interface_triangles'],1)
        self.assertEqual(proof['normal_charts'],2)
        np.testing.assert_array_equal(n['upper'][0],np.tile([0.,0.,-1.],(3,1)))
        np.testing.assert_array_equal(n['lower'][0],np.tile([0.,0.,1.],(3,1)))
        # Shared exterior vertices retain the same chart across both cells;
        # the dissected planar cap is not averaged into that exterior.
        np.testing.assert_array_equal(n['upper'][1,0],n['lower'][1,0])
        self.assertFalse(np.array_equal(n['upper'][0,0],n['upper'][1,0]))

    def test_duplicate_same_facing_and_third_owner_faces_are_rejected(self):
        p,cells=self.fixture()
        with self.assertRaises(ValueError):shared_tissue_corner_normals(p,{'a':cells['upper'],'b':cells['upper']})
        with self.assertRaises(ValueError):shared_tissue_corner_normals(p,{**cells,'extra':cells['lower'][:1]})


if __name__=='__main__':unittest.main()
