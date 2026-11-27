import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.hybrid_section_clearance import solve_hybrid_sections


class HybridSectionClearance(unittest.TestCase):
    def test_reserve_does_not_demand_unpromised_inherited_improvement(self):
        offset=np.array([[.12,0,0,0,-.04,0.]])
        d,proof=solve_hybrid_sections(sparse.eye(1),sparse.csr_matrix((0,3)),np.array([]),sparse.eye(3),np.zeros(3),
            maximum_motion_m=.001,distortion=sparse.csr_matrix((6,3)),distortion_offset=offset,
            distortion_bound=.12,numerical_cone_reserve=1e-7,time_limit_seconds=2.)
        self.assertIsNotNone(d)
        self.assertTrue(proof['feasible'])
        self.assertLessEqual(proof['maximumOriginalCorotatedResidual'],.12+1e-12)

    def test_preserves_permitted_inherited_state_rejected_by_frobenius(self):
        offset=np.array([[.12,0,0,0,-.04,0.]])
        C=sparse.csr_matrix((6,3))
        d,proof=solve_hybrid_sections(sparse.eye(1),sparse.csr_matrix((0,3)),np.array([]),sparse.eye(3),np.array([.0002,0,0]),
            maximum_motion_m=.001,distortion=C,distortion_offset=offset,distortion_bound=.12,time_limit_seconds=2.)
        self.assertIsNotNone(d);self.assertTrue(proof['feasible'])
        self.assertEqual(proof['mixedTangentCones']['exactSpectralFacets'],1)
        self.assertLessEqual(proof['maximumOriginalCorotatedResidual'],.12+1e-9)
        self.assertGreater(proof['maximumRawFrobeniusResidual'],.12)

    def test_rejects_actual_spectral_budget_violation(self):
        offset=np.array([[.13,0,0,0,-.04,0.]])
        d,proof=solve_hybrid_sections(sparse.eye(1),sparse.csr_matrix((0,3)),np.array([]),sparse.eye(3),np.zeros(3),
            maximum_motion_m=.001,distortion=sparse.csr_matrix((6,3)),distortion_offset=offset,distortion_bound=.12,time_limit_seconds=2.)
        self.assertIsNone(d);self.assertFalse(proof['feasible'])


if __name__=='__main__':unittest.main()
