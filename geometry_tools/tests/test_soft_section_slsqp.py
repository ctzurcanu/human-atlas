import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.soft_section_slsqp import solve_small_compartment


class SmallCompartmentTests(unittest.TestCase):
    def test_hard_recipient_interval_stops_a_soft_target(self):
        d, proof = solve_small_compartment(sparse.eye(1), sparse.csr_matrix([[-1., 0, 0]]),
                       np.array([-.00002]), sparse.csr_matrix([[1., 0, 0]]), np.array([.0001]),
                       sparse.csr_matrix((0, 3)), maximum_motion_m=.0002, distortion_bound=.05)
        self.assertTrue(proof['feasible'])
        self.assertLessEqual(d[0, 0], .0000200001)
        self.assertGreater(proof['maximumRemainingSoftSlack_m'], .000079)

    def test_facet_bound_is_enforced_even_when_motion_limit_allows_more(self):
        C = sparse.csr_matrix(np.array([[1000., 0, 0], [0, 0, 0], [0, 0, 0],
                                        [0, 0, 0], [0, 0, 0], [0, 0, 0]]))
        d, proof = solve_small_compartment(sparse.eye(1), sparse.csr_matrix((0, 3)), np.empty(0),
                       sparse.csr_matrix([[1., 0, 0]]), np.array([.0001]), C,
                       maximum_motion_m=.0002, distortion_bound=.05)
        self.assertTrue(proof['feasible'])
        self.assertLessEqual(d[0, 0], .0000500001)


if __name__ == '__main__':
    unittest.main()
