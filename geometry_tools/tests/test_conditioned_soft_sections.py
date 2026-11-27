import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.corotated_soft_section_clearance import solve_soft_sections

class ConditionedSoftSections(unittest.TestCase):
    def test_equivalent_scaled_rows_preserve_original_feasibility(self):
        A=np.array([[1.,0,0],[-1.,0,0],[0,1.,0],[0,-1.,0]])
        b=np.array([.0002,-.00025,-.00001,-.00001])
        scale=np.array([1e-5,1e6,1e-3,1e4])
        d,p=solve_soft_sections(sparse.csc_matrix((0,1)),sparse.csc_matrix(A*scale[:,None]),b*scale,
            sparse.csc_matrix([[1.,0,0]]),np.array([.0004]),maximum_motion_m=.001,
            condition_hard_rows=True)
        self.assertIsNotNone(d)
        self.assertTrue(p['originalUnscaledResidualsRequired'])
        self.assertLessEqual(np.maximum(b*scale-(A*scale[:,None])@d.ravel(),0).max(),1e-10)
        self.assertGreaterEqual(d[0,0],.0002-1e-10)
        self.assertLessEqual(d[0,0],.00025+1e-10)

    def test_inconsistent_original_rows_remain_rejected(self):
        d,p=solve_soft_sections(sparse.csc_matrix((0,1)),sparse.csc_matrix([[1e5,0,0],[-1e-4,0,0]]),
            np.array([50.,-4e-8]),sparse.csc_matrix([[1.,0,0]]),np.array([.0004]),
            maximum_motion_m=.001,condition_hard_rows=True)
        self.assertIsNone(d)
        self.assertFalse(p['feasible'])

if __name__=='__main__':unittest.main()
