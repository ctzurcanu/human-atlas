import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.conic_clearance_diagnostic import diagnose_direction


class DirectionalModelDiagnostic(unittest.TestCase):
    def test_explicit_fixed_body_halfspace_limits_inward_motion(self):
        A=sparse.csr_matrix([[-1.,0,0]])
        d,r=diagnose_direction(A,[-.0002],1,0,[1,0,0],maximum_motion_m=.001)
        self.assertTrue(r['optimalityCertified']);self.assertAlmostEqual(r['attainedDirectionalMotion_m'],.0002,places=9)
        self.assertGreater(r['hardRowDualWeights'][0],.9)

    def test_motion_ball_is_a_real_limit_without_a_body_barrier(self):
        d,r=diagnose_direction(sparse.csr_matrix((0,3)),[],1,0,[0,1,0],maximum_motion_m=.0004)
        self.assertAlmostEqual(r['attainedDirectionalMotion_m'],.0004,places=9)
        self.assertFalse(r['globalRoutingImpossibilityEstablished'])

if __name__=='__main__':unittest.main()
