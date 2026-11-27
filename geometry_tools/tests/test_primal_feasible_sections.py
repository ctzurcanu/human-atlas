import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.primal_feasible_section_clearance import solve_primal_feasible_sections
from geometry_tools.algorithms.constraint_generation_primal_sections import solve_sections_generated


class FeasibilityFirstSections(unittest.TestCase):
    def test_actual_point_and_metric_constraints_still_hold(self):
        A=sparse.csr_matrix([[-1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        b=np.array([-.0001,-.01,-.01]);C=sparse.csr_matrix(np.r_[np.eye(3),np.zeros((3,3))]*100)
        d,proof=solve_primal_feasible_sections(sparse.eye(1),A,b,sparse.eye(3),np.array([.0008,0.,0.]),maximum_motion_m=.001,
            distortion=C,distortion_bound=.075,time_limit_seconds=2.)
        self.assertIsNotNone(d);self.assertTrue(proof['feasible'])
        self.assertLessEqual(float(np.max(b-A@d.ravel())),1e-10)
        self.assertLessEqual(float(np.linalg.norm(C@d.ravel())),.075+1e-9)
        self.assertIn(proof['status'],['Solved','AlmostSolved','CallbackTerminated'])
        if proof['status']=='CallbackTerminated':self.assertTrue(proof['intentionalPrimalFeasibilityTermination'])

    def test_infeasible_problem_emits_no_point(self):
        A=sparse.csr_matrix([[1.,0.,0.],[-1.,0.,0.]])
        d,proof=solve_primal_feasible_sections(sparse.eye(1),A,np.array([.002,.002]),sparse.eye(3),np.zeros(3),maximum_motion_m=.001,time_limit_seconds=2.)
        self.assertIsNone(d);self.assertFalse(proof['feasible'])

    def test_full_original_constraint_generation_still_checks_all_rows(self):
        A=sparse.csr_matrix([[-1.,0.,0.]])
        d,proof=solve_sections_generated(sparse.eye(1),A,np.array([-.0001]),sparse.eye(3),np.array([.0008,0.,0.]),maximum_motion_m=.001,
            maximum_rounds=3,stage_seconds=2.,wall_seconds=8.)
        self.assertIsNotNone(d);self.assertTrue(proof['allOriginalHalfspacesChecked']);self.assertGreater(len(proof['rounds']),1)
        self.assertLessEqual(float(d[0,0]),.0001+1e-10)


if __name__=='__main__':unittest.main()
