import unittest
import numpy as np
from geometry_tools.algorithms.shared_simplex_contact import only_shared_simplex


class SharedSimplexTests(unittest.TestCase):
    a=np.array([[0.,0,0],[1,0,0],[0,1,0]])
    def test_coplanar_point_and_hidden_overlap(self):
        self.assertTrue(only_shared_simplex(self.a,[[0,0,0],[-1,0,0],[0,-1,0]])['passed'])
        self.assertFalse(only_shared_simplex(self.a,[[0,0,0],[.5,.2,0],[.2,.5,0]])['passed'])
    def test_edge_only_and_hidden_coplanar_overlap(self):
        self.assertTrue(only_shared_simplex(self.a,[[0,0,0],[1,0,0],[0,-1,0]])['passed'])
        self.assertFalse(only_shared_simplex(self.a,[[0,0,0],[1,0,0],[.2,.2,0]])['passed'])
    def test_noncoplanar_point_and_hidden_segment(self):
        self.assertTrue(only_shared_simplex(self.a,[[0,0,0],[0,0,1],[-1,-1,1]])['passed'])
        self.assertFalse(only_shared_simplex(self.a,[[0,0,0],[.2,.2,1],[.2,.2,-1]])['passed'])


if __name__=='__main__':unittest.main()
