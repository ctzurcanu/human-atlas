import unittest

import numpy as np

from geometry_tools.algorithms.tube_c2_course import C2Course


class C2AttachmentJets(unittest.TestCase):
    def test_owned_attachment_and_nonzero_root_curvature(self):
        s = np.linspace(0., 40., 12)
        reference = np.column_stack([s, .001 * s**2, np.zeros_like(s)])
        course = C2Course(reference, 40., len(s), 10., boundary_condition='not-a-knot')
        parameter = 17.3  # Deliberately not a uniform control station.
        jets = [(0., 1, [1., 0., 0.], 'root tangent'),
                (0., 2, [0., .002, 0.], 'root curvature'),
                (parameter, 0, [parameter, .001 * parameter**2, 0.], 'child attachment'),
                (parameter, 1, [1., .002 * parameter, 0.], 'child tangent')]
        controls, report = course.solve([], equalities=jets, time_limit=2.)
        self.assertIsNotNone(controls, report)
        self.assertTrue(report['certificate']['continuous_curvature_bound_pass'])
        for p, order, value, _ in jets:
            np.testing.assert_allclose(course.evaluate(controls, p, order), value, atol=1e-7, rtol=0)

    def test_incompatible_shared_attachment_is_rejected(self):
        s = np.linspace(0., 20., 8)
        reference = np.column_stack([s, np.zeros_like(s), np.zeros_like(s)])
        course = C2Course(reference, 20., len(s), 5.)
        controls, report = course.solve([], equalities=[
            (10., 0, [10., 0., 0.], 'parent attachment'),
            (10., 0, [10., 1., 0.], 'contradictory child attachment')], time_limit=2.)
        self.assertIsNone(controls)
        self.assertIn('Infeasible', report['status'])

    def test_legacy_endpoint_and_validation(self):
        s = np.linspace(0., 20., 8)
        reference = np.column_stack([s, np.zeros_like(s), np.zeros_like(s)])
        course = C2Course(reference, 20., len(s), 5.)
        controls, report = course.solve([], time_limit=2.)
        self.assertIsNotNone(controls, report)
        np.testing.assert_array_equal(controls[[0, -1]], reference[[0, -1]])
        self.assertTrue(report['certificate']['continuous_curvature_bound_pass'])
        with self.assertRaises(ValueError):
            course.solve([], equalities=[(21., 0, [0., 0., 0.])])


if __name__ == '__main__':
    unittest.main()
