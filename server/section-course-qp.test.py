"""Independent analytic energy, contradictory-side and fixed-port checks."""
from pathlib import Path
import sys, unittest
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.algorithms.section_course_clearance import smooth_course_metric
from geometry_tools.algorithms.section_course_qp import project_section_halfspaces, project_weighted_section_halfspaces
from geometry_tools.algorithms.triangle_contact_witness import intersection_witness


class SectionQPTests(unittest.TestCase):
    def test_single_constraint_matches_closed_form_metric_minimum(self):
        K = smooth_course_metric(17, smoothing=64., fixed=[0, 1, 15, 16])
        n = np.array([1., 2., -1.]); n /= np.linalg.norm(n)
        d, proof = project_section_halfspaces(K, [8], [n], [.0017])
        expected = K[:, 8, None] * n * (.0017 / K[8, 8])
        np.testing.assert_allclose(d, expected, rtol=0., atol=2e-10)
        np.testing.assert_array_equal(d[[0, 1, 15, 16]], np.zeros((4, 3)))
        self.assertLessEqual(proof['maximumResidualMeters'], 1e-10)

    def test_obstacle_side_conflict_is_distinguished_from_slow_convergence(self):
        K = smooth_course_metric(9)
        with self.assertRaisesRegex(ValueError, 'Incompatible local obstacle halfspaces'):
            project_section_halfspaces(K, [4, 4], [[1., 0., 0.], [-1., 0., 0.]], [.001, .001])

    def test_cannot_release_a_fixed_attachment(self):
        K = smooth_course_metric(9, fixed=[0, 8])
        with self.assertRaisesRegex(ValueError, 'protected source station'):
            project_section_halfspaces(K, [0], [[1., 0., 0.]], [.001])

    def test_contact_between_two_sections_has_analytic_metric_solution(self):
        K = smooth_course_metric(13, fixed=[0, 12])
        w = np.zeros(13); w[5] = .25; w[6] = .75
        n = np.array([0., 1., 0.])
        d, proof = project_weighted_section_halfspaces(K, [w], [n], [.002])
        expected = (K @ w)[:, None] * n * (.002 / (w @ K @ w))
        np.testing.assert_allclose(d, expected, rtol=0., atol=2e-10)
        self.assertLessEqual(proof['maximumResidualMeters'], 1e-10)

    def test_proper_contact_witness_stays_on_both_actual_material_triangles(self):
        a = np.array([[-2., -2., 0.], [2., -2., 0.], [0., 2., 0.]])
        b = np.array([[0., -1., -1.], [0., -1., 1.], [0., 1., 0.]])
        p, wa, wb = intersection_witness(a, b)
        np.testing.assert_allclose(p, [0., 0., 0.], rtol=0., atol=1e-14)
        np.testing.assert_allclose(wa @ a, p, rtol=0., atol=1e-14)
        np.testing.assert_allclose(wb @ b, p, rtol=0., atol=1e-14)


if __name__ == '__main__':
    unittest.main()
