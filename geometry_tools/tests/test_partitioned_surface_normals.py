import unittest
import numpy as np
from geometry_tools.algorithms.partitioned_surface_normals import transport_partitioned_vertex_words

class PartitionedNormalTests(unittest.TestCase):
    def test_collapsed_inner_corner_cannot_change_stationary_outer_normals(self):
        v=np.array([[0.,0.,0.],[.001,0.,0.],[0.,.001,0.],
                    [0.,0.,0.],[.001,0.,-.00001],[0.,.001,-.00001]])
        f=np.array([[0,1,2],[3,5,4]])
        words=np.array([[0,0,32767]]*3+[[0,0,-32767]]*3,dtype=np.int16)
        q=v.copy();q[3,2]=-.00001
        result,proof=transport_partitioned_vertex_words(v,q,f,words,[[0],[1]],maximum_facet_rotation_degrees=2)
        np.testing.assert_array_equal(result[:3],words[:3])
        self.assertTrue(proof['noCrossSheetNormalStarPooling'])
        self.assertTrue(proof['partitions'][0]['allPartitionCoordinatesFixed'])

    def test_undeclared_or_shared_face_ownership_is_rejected(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        f=np.array([[0,1,2],[0,3,1]]);words=np.tile(np.array([0,0,32767],np.int16),(4,1))
        with self.assertRaises(AssertionError):
            transport_partitioned_vertex_words(v,v,f,words,[[0],[1]],maximum_facet_rotation_degrees=2)
        with self.assertRaises(AssertionError):
            transport_partitioned_vertex_words(v,v,f,words,[[0]],maximum_facet_rotation_degrees=2)

if __name__=='__main__':unittest.main()
