import unittest
import numpy as np
from geometry_tools.urinary_delivery import stored_source_face_lineage


class StoredSourceLineageTests(unittest.TestCase):
    def test_collapsed_reference_fragment_does_not_shift_the_remaining_parent(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,0.],[2.,2.,0.]])
        source=np.array([[0,3,1],[0,1,2]])
        parents=stored_source_face_lineage(v,np.array([[2,1,0]]),source,np.array([7,9]))
        np.testing.assert_array_equal(parents,[9])
        with self.assertRaises(ValueError):stored_source_face_lineage(v,np.array([[0,1,4]]),source,np.array([7,9]))

    def test_internal_retriangulation_cannot_inherit_a_false_anatomical_source_face(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.]])
        original=np.array([[0,1,2],[0,2,3]]);new=np.array([[0,1,3],[1,2,3]])
        np.testing.assert_array_equal(stored_source_face_lineage(v,new,original,np.array([-1,-1])),[-1,-1])
        with self.assertRaises(ValueError):stored_source_face_lineage(v,new,original,np.array([7,9]))


if __name__=='__main__':unittest.main()
