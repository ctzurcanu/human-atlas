import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.constraint_generation_sections import solve_sections_generated


class GeneratedConstraints(unittest.TestCase):
    def test_unseeded_constraint_is_recovered_and_all_rows_checked(self):
        A=sparse.csr_matrix([[-1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        b=np.array([-.0001,-.01,-.01]);B=sparse.eye(3)
        d,proof=solve_sections_generated(sparse.eye(1),A,b,B,np.array([.0008,0.,0.]),maximum_motion_m=.001,
            maximum_rounds=3,stage_seconds=2,wall_seconds=8)
        self.assertIsNotNone(d)
        self.assertGreater(len(proof['rounds']),1)
        self.assertLessEqual(float(d[0,0]),.0001+1e-10)
        self.assertTrue(proof['allOriginalHalfspacesChecked'])
        self.assertLessEqual(float(np.max(b-A@d.ravel())),1e-10)


if __name__=='__main__':unittest.main()
