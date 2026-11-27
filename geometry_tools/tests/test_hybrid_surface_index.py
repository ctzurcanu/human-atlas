import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.hybrid_surface_index import HybridSurfaceIndex


class HybridSurfaceIndexTests(unittest.TestCase):
    def test_tiny_triangle_keeps_exact_witness_and_blocks_segment_after_native_collapse(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],
                    [.75,.75,.5],[.75+1e-10,.75,.5],[.75,.75+1e-10,.5]])
        m=Mesh(v,np.array([[0,1,2],[3,4,5]]),'atlas-world');before=v.copy();index=HybridSurfaceIndex(m)
        self.assertEqual(len(index.fallback_ids),1)
        q=np.array([[.75+2e-11,.75+2e-11,.5+1e-9]])
        hit,proof=index.closest_points(q);self.assertEqual(hit['face_ids'][0],1)
        self.assertAlmostEqual(hit['distance_m'][0],1e-9,delta=1e-15)
        self.assertFalse(index.segment_clear(q,q-[0,0,2e-9])[0])
        self.assertTrue(proof['allSourceFacesRetained']);np.testing.assert_array_equal(m.vertices,before)


if __name__=='__main__':unittest.main()
