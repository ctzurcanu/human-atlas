import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.minimax_section_clearance import solve_minimax_sections

class MinimaxSections(unittest.TestCase):
    def test_unattainable_full_clearance_has_explicit_best_deficit(self):
        A=sparse.csc_matrix([[-1.,0,0]])
        C=sparse.csc_matrix([[1.,0,0]])
        d,p=solve_minimax_sections(sparse.csc_matrix((0,1)),A,np.array([-.0001]),C,np.array([.0004]),
            clearance_rows=C,clearance_rhs=np.array([.0004]),maximum_motion_m=.001,
            condition_hard_rows=True,soft_weight=10000.)
        self.assertIsNotNone(d)
        self.assertLessEqual(d[0,0],.0001+1e-10)
        self.assertAlmostEqual(p['worstClearanceDeficit_m'],.0003,places=8)
        self.assertLessEqual(p['maximumEpigraphResidual_m'],1e-10)

    def test_fixed_context_conflict_is_still_rejected(self):
        A=sparse.csc_matrix([[1.,0,0],[-1.,0,0]])
        C=sparse.csc_matrix([[1.,0,0]])
        d,p=solve_minimax_sections(sparse.csc_matrix((0,1)),A,np.array([.0002,-.0001]),C,np.array([.0004]),
            clearance_rows=C,clearance_rhs=np.array([.0004]),maximum_motion_m=.001,condition_hard_rows=True)
        self.assertIsNone(d)

    def test_clearance_priority_cannot_be_exchanged_for_smoothing(self):
        A=sparse.csc_matrix([[-1.,0,0]]);C=sparse.csc_matrix([[1.,0,0]])
        d,p=solve_minimax_sections(sparse.csc_matrix([[1000.]]),A,np.array([-.0001]),C,np.array([.0004]),
            clearance_rows=C,clearance_rhs=np.array([.0004]),maximum_motion_m=.001,
            condition_hard_rows=True,soft_weight=10000.,lexicographic=True)
        self.assertIsNotNone(d)
        self.assertTrue(p['firstPriority']['optimalityCertified'])
        self.assertGreaterEqual(d[0,0],.0001-1.1e-8)
        self.assertLessEqual(p['worstClearanceDeficit_m'],.0003+1.1e-8)

if __name__=='__main__':unittest.main()
