import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.lazy_source_shape_sections import solve_lazy_source_shape


class LazyOriginalSourceShape(unittest.TestCase):
    def test_full_shape_replay_adds_original_facet_budget(self):
        C = sparse.csr_matrix([[750., 0., 0.]] + [[0., 0., 0.]] * 5)
        def check(d):
            strain = abs(750. * d[0, 0])
            return np.array([0] if strain > .075 + 1e-8 else [], dtype=int), dict(strain=float(strain))
        d, proof = solve_lazy_source_shape(sparse.eye(1), sparse.csr_matrix((0, 3)), np.array([]),
            sparse.eye(3), np.array([.0008, 0., 0.]), all_distortion=C, all_offsets=np.zeros((1, 6)),
            exact_shape_violations=check, maximum_motion_m=.001, maximum_stages=3, wall_seconds=8.)
        self.assertIsNotNone(d)
        self.assertTrue(proof['fullOriginalSourceShapeChecked'])
        self.assertEqual(proof['stages'][0]['activeOriginalFacetCones'], 0)
        self.assertEqual(proof['stages'][1]['activeOriginalFacetCones'], 1)
        self.assertLessEqual(abs(750. * d[0, 0]), .075 + 1e-8)

    def test_exact_pose_failure_cannot_be_accepted_by_linear_cone(self):
        C = sparse.csr_matrix((6, 3))
        def still_invalid(d):
            return np.array([0]), dict(exactPoseBudgetViolated=True)
        d, proof = solve_lazy_source_shape(sparse.eye(1), sparse.csr_matrix((0, 3)), np.array([]),
            sparse.eye(3), np.array([.0008, 0., 0.]), all_distortion=C, all_offsets=np.zeros((1, 6)),
            exact_shape_violations=still_invalid, maximum_motion_m=.001, maximum_stages=3, wall_seconds=8.)
        self.assertIsNone(d)
        self.assertFalse(proof['feasible'])
        self.assertIn('relinearization required', proof['reason'])


if __name__ == '__main__':
    unittest.main()
