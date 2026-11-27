import unittest
import numpy as np

from geometry_tools.algorithms.quantized_surface_topology import contract_declared_intersection_edges


class DeclaredIntersectionContractionTests(unittest.TestCase):
    def fixture(self):
        p = np.array([[0., 0., 0.], [.001, 0., 0.], [.001, 1e-10, 0.], [.002, .001, 0.], [.001, .002, 0.]])
        q = p.copy(); q[2] = p[1]
        f = np.array([[1, 2, 3], [2, 4, 3]])
        return p, q, {'first': f, 'second': f[:, [0, 2, 1]]}

    def test_identification_is_reciprocal_and_keeps_parent_faces(self):
        p, q, surfaces = self.fixture()
        literal = p.copy()
        output, ancestry, evidence = contract_declared_intersection_edges(p, surfaces, q, [[1, 2]], protected_prefix=1)
        np.testing.assert_array_equal(output['first'], output['second'][:, [0, 2, 1]])
        np.testing.assert_array_equal(ancestry['first'], [1])
        np.testing.assert_array_equal(p, literal)
        self.assertFalse(evidence['literal_source_controls_moved'])

    def test_source_controls_and_close_nonedge_controls_are_protected(self):
        p, q, surfaces = self.fixture()
        with self.assertRaisesRegex(ValueError, 'new intersection'):
            contract_declared_intersection_edges(p, surfaces, q, [[1, 2]], protected_prefix=2)
        with self.assertRaisesRegex(ValueError, 'construction edge'):
            contract_declared_intersection_edges(p, surfaces, q, [[1, 4]], protected_prefix=1)

    def test_short_edge_alone_does_not_authorize_contraction(self):
        p, _, surfaces = self.fixture()
        with self.assertRaisesRegex(ValueError, 'combined storage error'):
            contract_declared_intersection_edges(p, surfaces, p.copy(), [[1, 2]], protected_prefix=1)


if __name__ == '__main__':
    unittest.main()
