import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.original_singular_pose_rows import OriginalSingularPoseRows


class OriginalSingularPoses(unittest.TestCase):
    def model(self, common=False):
        p = np.array([[0., 0., 0.], [.01, 0., 0.], [0., .008, .001]])
        q = p*np.array([1.03, .98, 1.])
        f = np.array([[0, 1, 2]])
        centers = np.array([[0., -.003, 0.]]) if common else p+np.array([.001, -.002, .001])
        groups = np.zeros(3, int) if common else np.arange(3)
        return OriginalSingularPoseRows(Mesh(p, f, 'test'), Mesh(q, f, 'test'), groups,
                                        centers, np.full(len(centers), .03))

    def test_analytic_derivative_matches_exact_rotated_source(self):
        m = self.model()
        x = np.array([[.00004, -.00001, .00002], [-.00003, .00002, -.00001],
                      [.00001, -.00003, .00004], [.0002, -.0001, .0003],
                      [-.0001, .0003, .0002], [.0002, .0001, -.0002]])
        values, gradient = m.singular_gradients(x, [0])
        step = 1e-8
        finite = np.zeros((2, x.size))
        for column in range(x.size):
            plus, minus = x.copy(), x.copy()
            plus.ravel()[column] += step; minus.ravel()[column] -= step
            finite[:, column] = ((m.values(plus)-m.values(minus))/(2*step)).ravel()
        np.testing.assert_allclose(gradient.toarray(), finite, rtol=2e-6, atol=2e-7)

    def test_common_proper_pose_preserves_source_shape_and_derivative(self):
        m = self.model(common=True)
        x = np.array([[.0002, -.0001, .0003], [.001, -.0004, .0002]])
        np.testing.assert_allclose(m.values(x), m.initial_singular, atol=1e-14)
        _, gradient = m.singular_gradients(x, [0])
        self.assertLess(abs(gradient.toarray()).max(), 1e-11)

    def test_rows_reproduce_actual_budget_at_linearization_point(self):
        m = self.model()
        x = np.zeros((6, 3)); x[1, 0] = .001
        values = m.values(x).ravel()
        rows, b, _ = m.local_rows(x, [0])
        residual = b-rows@x.ravel()
        self.assertEqual(bool(np.any(residual > 0)), bool(np.any(abs(values-1) > m.budget[0])))
        bad, proof = m.violations(x)
        self.assertEqual(bad.tolist(), [0])
        self.assertTrue(proof['originalReferenceAndBudgetsUnchanged'])

    def test_small_facet_row_units_retain_stricter_shape_precision(self):
        m = self.model()
        x = np.zeros((6, 3)); x[1, 0] = .001
        _, gradient = m.singular_gradients(x, [0])
        rows, _, proof = m.local_rows(x, [0])
        np.testing.assert_allclose(rows.toarray(),
            np.r_[-gradient.toarray(), gradient.toarray()]/10., atol=1e-12)
        self.assertTrue(proof['positiveRowScalingDoesNotChangeInequalities'])


if __name__ == '__main__':
    unittest.main()
