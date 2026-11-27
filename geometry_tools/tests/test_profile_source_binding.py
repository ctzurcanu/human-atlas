import unittest
import numpy as np
from geometry_tools.algorithms.profile_source_binding import verify_profile_topology_replay


class ProfileSourceBindingTests(unittest.TestCase):
    def test_repaired_diagonal_keeps_the_same_owned_boundary(self):
        source=np.array([[0,1,2],[0,2,3],[4,5,6]])
        repaired=np.array([[1,2,3],[1,3,0],[4,5,6]])
        result=verify_profile_topology_replay(source,repaired)
        self.assertEqual(result['changedFaces'],2)
        self.assertEqual(result['changedPatches'][0]['faceIds'],[0,1])

    def test_reversed_boundary_or_different_source_control_is_rejected(self):
        source=np.array([[0,1,2],[0,2,3]])
        for other in [np.array([[3,2,1],[0,3,1]]),np.array([[1,2,4],[1,4,0]])]:
            with self.assertRaises(ValueError):verify_profile_topology_replay(source,other)

    def test_face_removal_is_not_a_cache_replay(self):
        with self.assertRaises(ValueError):verify_profile_topology_replay(np.array([[0,1,2],[0,2,3]]),np.array([[0,1,2]]))
