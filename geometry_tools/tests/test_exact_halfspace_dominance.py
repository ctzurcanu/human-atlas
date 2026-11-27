import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.exact_halfspace_dominance import exact_halfspace_dominance


class ExactDominance(unittest.TestCase):
    def test_only_identical_rows_and_strongest_original_bound(self):
        A=sparse.csr_matrix([[1.,2.,0.],[1.,2.,0.],[1.,2.+1e-12,0.],[0.,0.,0.],[0.,0.,0.]])
        b=np.array([1.,2.,1.5,-2.,-1.]);C,d,m,proof=exact_halfspace_dominance(A,b)
        self.assertEqual(len(d),3)
        np.testing.assert_array_equal(C.toarray()[m],A.toarray())
        self.assertTrue(np.all(d[m]>=b))
        self.assertFalse(proof['coefficientRoundingUsed'])


if __name__=='__main__':unittest.main()
