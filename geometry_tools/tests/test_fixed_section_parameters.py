import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.fixed_section_parameters import solve_with_fixed_sections
from geometry_tools.algorithms.corotated_soft_section_clearance import solve_soft_sections

class FixedSectionParameters(unittest.TestCase):
    def test_fixed_variables_are_literal_zero_and_coupling_is_retained(self):
        S=sparse.csc_matrix([[2.,-3.]])
        A=sparse.csc_matrix([[0,0,0,1.,0,0]])
        goal=sparse.csc_matrix([[1.,0,0,1.,0,0]])
        d,p=solve_with_fixed_sections(solve_soft_sections,S,A,np.array([.0001]),goal,np.array([.0002]),
            fixed=[True,False],maximum_motion_m=.001,condition_hard_rows=True)
        self.assertIsNotNone(d)
        np.testing.assert_array_equal(d[0],np.zeros(3))
        self.assertGreaterEqual(d[1,0],.0001-1e-10)
        self.assertEqual(p['declaredFixedSectionVectors'],1)

    def test_impossible_fixed_row_is_not_dropped(self):
        d,p=solve_with_fixed_sections(solve_soft_sections,sparse.csc_matrix((0,2)),
            sparse.csc_matrix([[1.,0,0,0,0,0]]),np.array([.0001]),
            sparse.csc_matrix([[0,0,0,1.,0,0]]),np.array([.0002]),
            fixed=[True,False],maximum_motion_m=.001,condition_hard_rows=True)
        self.assertIsNone(d)

if __name__=='__main__':unittest.main()
