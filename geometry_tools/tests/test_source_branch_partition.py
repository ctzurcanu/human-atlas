import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_branch_partition import partition_source_faces

class PartitionTests(unittest.TestCase):
    def test_faces_and_packed_attributes_have_single_ownership(self):
        v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[2,0,0],[2,1,0]],float)
        f=np.array([[0,1,2],[0,2,3],[1,4,5],[1,5,2]])
        attrs=dict(normals=np.tile(np.array([0,0,32767],np.int16),(len(v),1)),normalSigns=np.ones(len(v),np.int8))
        m=Mesh(v,f,'test:source',attrs)
        pieces=partition_source_faces(m,[[2,3]])
        self.assertEqual(sum(len(p[0].faces)for p in pieces),len(f))
        for child,ids,raw in pieces:
            np.testing.assert_array_equal(child.triangles,m.triangles[ids])
            for key in attrs:np.testing.assert_array_equal(child.attributes[key],attrs[key][raw])
        # Shared border uses identical original coordinates, without a cap.
        self.assertEqual(set(map(tuple,pieces[0][0].vertices))&set(map(tuple,pieces[1][0].vertices)),{(1.,0.,0.),(1.,1.,0.)})

    def test_overlap_is_rejected(self):
        m=Mesh(np.array([[0,0,0],[1,0,0],[0,1,0]],float),np.array([[0,1,2]]),'test:source')
        with self.assertRaises(ValueError):partition_source_faces(m,[[0],[0]])

if __name__=='__main__':unittest.main()
