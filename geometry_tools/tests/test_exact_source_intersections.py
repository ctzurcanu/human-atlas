from fractions import Fraction as Q
import unittest
import numpy as np
from geometry_tools.algorithms.exact_source_intersections import edge_plane,contains
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank


class ExactSourceTests(unittest.TestCase):
    def test_edge_plane_and_finite_triangle(self):
        tri=((Q(1,2),Q(-1),Q(-1)),(Q(1,2),Q(1),Q(-1)),(Q(1,2),Q(0),Q(1)))
        hit=edge_plane((Q(0),Q(0),Q(0)),(Q(1),Q(0),Q(0)),tri)
        self.assertEqual(hit,(Q(1,2),Q(0),Q(0)))
        self.assertTrue(contains(hit,tri))
        self.assertFalse(contains((Q(1,2),Q(2),Q(0)),tri))

    def test_distinct_cut_features_remain_distinct(self):
        bank=CommonPositionBank(np.array([[0.,0.,0.]]),0.)
        a=bank.insert([.5,0.,0.]);b=bank.insert([.5+2**-48,0.,0.])
        self.assertNotEqual(a,b)
        self.assertEqual(bank.insert([.5,0.,0.]),a)
        self.assertEqual(bank.maximum_snap,0.)


if __name__=='__main__':unittest.main()
