import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.support_exchange_sections import solve_sections_generated


class SupportExchange(unittest.TestCase):
    def test_full_replay_catches_other_inequality_on_same_motion_support(self):
        # Both halfspaces act on the same 3D section; satisfying the initially
        # worse row does not imply satisfying the other, differently oriented row.
        A = sparse.csr_matrix([[-1., -.2, 0.], [-.2, -1., 0.]])
        b = np.array([-.00015, -.0002])
        d, proof = solve_sections_generated(sparse.eye(1), A, b, sparse.eye(3),
            np.array([.0008, .0007, 0.]), maximum_motion_m=.001,
            maximum_rounds=5, stage_seconds=1., wall_seconds=8.)
        self.assertIsNotNone(d)
        self.assertTrue(proof['allOriginalHalfspacesChecked'])
        self.assertEqual(proof['fullOriginalRows'], 2)
        self.assertGreaterEqual(len(proof['rounds']), 3)
        self.assertLessEqual(float((b - A @ d.ravel()).max()), 1e-10)

    def test_exhausted_exchange_never_emits_subset_solution(self):
        A = sparse.csr_matrix([[-1., 0., 0.]])
        d, proof = solve_sections_generated(sparse.eye(1), A, np.array([-.0001]),
            sparse.eye(3), np.array([.0008, 0., 0.]), maximum_motion_m=.001,
            maximum_rounds=1, stage_seconds=1., wall_seconds=3.)
        self.assertIsNone(d)
        self.assertFalse(proof['feasible'])
        self.assertEqual(proof['rounds'][0]['violatedOriginalRows'], 1)


if __name__ == '__main__':
    unittest.main()
