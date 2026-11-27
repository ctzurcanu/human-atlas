import unittest
import numpy as np
from geometry_tools.algorithms.source_face_ancestry import retained_source_face_ancestry


class AncestryTests(unittest.TestCase):
    def test_removed_faces_and_compacted_controls_are_replayed(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[2,0,0],[2,1,0]],float)
        f=np.array([[0,1,2],[1,3,2],[1,4,3],[4,5,3]])
        raw=np.array([1,3,4,5]);g=np.array([[0,2,1],[2,3,1]])
        kept,pairs,proof=retained_source_face_ancestry(v,f,v[raw],g)
        np.testing.assert_array_equal(kept,[2,3]);np.testing.assert_array_equal(pairs,[[0,1]])
        self.assertEqual(proof['removedOriginalFaceIds'],[0,1])

    def test_wrong_source_frame_is_rejected(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0]],float);f=np.array([[0,1,2],[2,1,0]])
        with self.assertRaises(ValueError):retained_source_face_ancestry(v,f,v+.001,f)

    def test_reversed_face_is_not_an_original_face(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0]],float);f=np.array([[0,1,2],[1,3,2]])
        with self.assertRaises(ValueError):retained_source_face_ancestry(v,f,v,f[:,::-1])


if __name__=='__main__':unittest.main()
