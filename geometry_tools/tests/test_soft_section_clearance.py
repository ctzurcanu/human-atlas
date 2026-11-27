import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.soft_section_clearance import solve_soft_sections,prune_motion_ball_rows


class SoftSectionClearanceTests(unittest.TestCase):
    def test_motion_ball_pruning_retains_every_reachable_constraint(self):
        A=sparse.csr_matrix([[1,0,0],[0,1,0],[0,0,0]])
        B,b,p=prune_motion_ball_rows(A,np.array([-.002,-.0005,-.01]),.001)
        self.assertEqual(p['droppedGuaranteedRows'],2)
        np.testing.assert_array_equal(B.toarray(),[[0,1,0]])
        np.testing.assert_array_equal(b,[-.0005])

    def test_soft_goal_cannot_override_a_fixed_context(self):
        hard=sparse.csc_matrix([[-1.,0,0]])
        goal=sparse.csc_matrix([[1.,0,0]])
        d,p=solve_soft_sections(sparse.csc_matrix((0,1)),hard,np.array([-.00002]),goal,np.array([.0001]),maximum_motion_m=.001)
        self.assertIsNotNone(d)
        self.assertLessEqual(d[0,0],.0000200001)
        self.assertGreater(p['maximumRemainingSlack_m'],.000079)
        self.assertTrue(p['softGoalsAreNotClearanceAcceptance'])

    def test_reachable_goal_and_motion_ball(self):
        d,p=solve_soft_sections(sparse.csc_matrix((0,1)),sparse.csc_matrix((0,3)),np.empty(0),sparse.csc_matrix([[1.,0,0]]),np.array([.0001]),maximum_motion_m=.0002)
        self.assertTrue(p['feasible'])
        self.assertAlmostEqual(d[0,0],.0001,delta=2e-8)
        self.assertLessEqual(np.linalg.norm(d[0]),.0002)

    def test_infeasible_hard_context_emits_no_candidate(self):
        d,p=solve_soft_sections(sparse.csc_matrix((0,1)),sparse.csc_matrix([[1.,0,0],[-1.,0,0]]),np.array([.001,.001]),sparse.csc_matrix([[1.,0,0]]),np.array([.0001]),maximum_motion_m=.001)
        self.assertIsNone(d)
        self.assertFalse(p['feasible'])

if __name__=='__main__':unittest.main()
