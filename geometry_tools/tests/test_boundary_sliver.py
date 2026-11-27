import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.boundary_sliver import remove_dangling_sliver_sheet
from geometry_tools.algorithms.tube_continuity import boundary_loops


class BoundarySliverTests(unittest.TestCase):
    def fixture(self,width=.000001):
        v=np.array([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0],[-1,-1,2],[1,-1,2],[1,1,2],[-1,1,2]],float)*.001
        v=np.r_[v,v[[0,5]]+[0,0,width]]
        f=np.array([[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7],[0,8,5],[5,8,9]])
        return Mesh(v,f,'fixture',{'uvs':np.zeros((len(v),2))})

    def test_only_declared_microscopic_appendage_is_removed(self):
        old=self.fixture();new,p=remove_dangling_sliver_sheet(old,[8,9],maximum_width_m=.000005)
        np.testing.assert_array_equal(new.vertices,old.vertices)
        np.testing.assert_array_equal(new.faces,old.faces[:8])
        np.testing.assert_array_equal(new.attributes['uvs'],old.attributes['uvs'])
        self.assertEqual(sorted(map(len,boundary_loops(new.faces))),[4,4])
        self.assertEqual(p['newCapsAdded'],0);self.assertEqual(p['newFacesAdded'],0)
        self.assertLessEqual(p['maximumEndpointWidth_m'],.000005)

    def test_a_broad_sheet_is_not_silently_deleted(self):
        with self.assertRaises(ValueError):remove_dangling_sliver_sheet(self.fixture(.0001),[8,9],maximum_width_m=.000005)

    def test_real_tube_faces_are_not_classified_as_dangling(self):
        with self.assertRaises(ValueError):remove_dangling_sliver_sheet(self.fixture(),[0,1],maximum_width_m=.000005)

    def test_invalid_face_selection_is_rejected(self):
        for ids,width in [([8,8],.000005),([100],.000005),([8,9],.001)]:
            with self.subTest(ids=ids,width=width),self.assertRaises(ValueError):remove_dangling_sliver_sheet(self.fixture(),ids,maximum_width_m=width)
