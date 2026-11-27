import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.nonregressing_singular_metric_rows import nonregressing_singular_metric_rows
class SingularBudgetTests(unittest.TestCase):
    def test_two_axis_strain_cannot_spend_a_larger_frobenius_allowance(self):
        f=np.array([[0,1,2]]);a=Mesh(np.array([[0.,0,0],[1,0,0],[0,1,0]]),f,'fixture:metres')
        b=Mesh(np.array([[0.,0,0],[1.2,0,0],[0,1.1,0]]),f,'fixture:metres')
        C,o,p=nonregressing_singular_metric_rows(a,b,np.arange(3),3)
        self.assertAlmostEqual(p['sourceFacetBudgets'][0],.2)
        self.assertGreater(np.linalg.norm(o),.075)
        d=np.zeros((3,3));d[2,1]=-.1
        self.assertAlmostEqual(np.linalg.norm((C@d.ravel()).reshape(-1,6)+o),.075)
        self.assertTrue(p['sourceReferenceNeverReset'])
if __name__=='__main__':unittest.main()
