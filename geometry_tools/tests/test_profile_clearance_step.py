import unittest
import numpy as np
from geometry_tools.algorithms.profile_clearance_step import clearance_translation_step


class ProfileClearanceStepTest(unittest.TestCase):
    def test_preserves_fixed_sections_and_meets_two_surface_constraints(self):
        r = clearance_translation_step(2, [[0, 1]], [1, 1],
                                      [[1, 0, 0], [0, -1, 0]], [.001, .0005],
                                      fixed_profiles=[0])
        np.testing.assert_array_equal(r['translations_m'][0], [0, 0, 0])
        self.assertLess(r['maximum_residual_m'], 1e-12)
        self.assertGreaterEqual(r['translations_m'][1, 0], .001-1e-12)
        self.assertLessEqual(r['translations_m'][1, 1], -.0005+1e-12)

    def test_conflicting_corridor_returns_violation_without_resizing(self):
        r = clearance_translation_step(1, [], [0, 0],
                                      [[1, 0, 0], [-1, 0, 0]], [.001, .001])
        self.assertGreaterEqual(sum(r['inequality_residual_m']), .002-1e-12)
        self.assertFalse(r['profile_radius_changed'])
        self.assertFalse(r['surface_containment_proved'])

    def test_trust_boundary_is_not_silently_extended(self):
        r = clearance_translation_step(1, [], [0], [[0, 0, 1]], [.004],
                                      maximum_step_m=.001)
        self.assertAlmostEqual(r['translations_m'][0, 2], .001)
        self.assertAlmostEqual(r['maximum_residual_m'], .003)

    def test_invalid_witness_normal_is_rejected(self):
        with self.assertRaises(ValueError):
            clearance_translation_step(1, [], [0], [[2, 0, 0]], [0])


if __name__ == '__main__':
    unittest.main()
