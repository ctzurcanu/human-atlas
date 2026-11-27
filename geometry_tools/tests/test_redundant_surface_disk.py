import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.redundant_surface_disk import remove_redundant_surface_disk

def tetra(extra=True):
    v=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.01,0.],[0.,0.,.01],[.002,.002,.001]])
    f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]+([[0,1,4],[1,2,4],[2,0,4]] if extra else []))
    return Mesh(v,f,'test:metres')

class RedundantSurfaceDiskTests(unittest.TestCase):
    def test_extra_disk_only_is_removed(self):
        mesh=tetra();keep,p=remove_redundant_surface_disk(mesh,[4,5,6])
        np.testing.assert_array_equal(keep,[0,1,2,3]);self.assertEqual(p['facesAdded'],0)
        self.assertTrue(p['sourceOpenBoundaryEdgesLiteral']);np.testing.assert_array_equal(mesh.vertices,tetra().vertices)

    def test_true_outer_face_not_replaced_by_inconsistent_extra_disk(self):
        with self.assertRaisesRegex(ValueError,'winding'):remove_redundant_surface_disk(tetra(),[0])

    def test_regular_outer_face_is_not_an_extra_sheet(self):
        with self.assertRaisesRegex(ValueError,'opening|triple'):remove_redundant_surface_disk(tetra(False),[0])

if __name__=='__main__':unittest.main()
