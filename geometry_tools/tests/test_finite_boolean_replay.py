import unittest
import numpy as np
from geometry_tools.algorithms.finite_tissue_boolean import same_oriented_triangles


class OrientedReplayTests(unittest.TestCase):
    def test_order_and_cyclic_rotation(self):
        self.assertTrue(same_oriented_triangles(
            [[0,1,2],[2,1,3]], [[3,2,1],[1,2,0]]))

    def test_reverse_face_rejected(self):
        self.assertFalse(same_oriented_triangles([[0,1,2]], [[0,2,1]]))

    def test_duplicate_coverage_rejected(self):
        self.assertFalse(same_oriented_triangles(
            [[0,1,2],[2,1,3]], [[0,1,2],[0,1,2]]))

    def test_same_outer_edges_do_not_prove_same_triangles(self):
        self.assertFalse(same_oriented_triangles(
            [[0,1,2],[0,2,3]], [[0,1,3],[1,2,3]]))


if __name__=='__main__':unittest.main()
