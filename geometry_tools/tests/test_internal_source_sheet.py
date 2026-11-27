import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.internal_source_sheet import remove_internal_source_sheet

def fixture(outside=False):
    v=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.01,0.],[0.,0.,.01],
                [.002,-.002 if outside else .002,.002]])
    f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3],[0,1,4]])
    return Mesh(v,f,'test:metres',{'colors':np.arange(15,dtype=np.uint8).reshape(5,3)})

class InteriorSourceSheetTests(unittest.TestCase):
    def test_only_declared_internal_partition_removed(self):
        m=fixture();keep,p=remove_internal_source_sheet(m,[4])
        np.testing.assert_array_equal(keep,[0,1,2,3])
        self.assertEqual(p['newBoundaryEdges'],0)
        self.assertEqual(p['retainedTopology']['boundaryLoops'],0)
        self.assertEqual(p['computationalCapsEmitted'],0)
        np.testing.assert_array_equal(m.vertices,fixture().vertices)
        np.testing.assert_array_equal(m.attributes['colors'],fixture().attributes['colors'])

    def test_outside_flap_not_classified_as_internal(self):
        with self.assertRaisesRegex(ValueError,'strictly inside'):
            remove_internal_source_sheet(fixture(True),[4])

    def test_exterior_wall_must_not_be_deleted(self):
        with self.assertRaises(ValueError):remove_internal_source_sheet(fixture(),[0])

if __name__=='__main__':unittest.main()
