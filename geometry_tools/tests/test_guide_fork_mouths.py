import unittest

import numpy as np

from geometry_tools.commands.pose_hand_vessel_courses import separate_guide_mouths, trim_course, retain_proximal_course


class GuideForkMouthTests(unittest.TestCase):
    def match(self, sign, component, radius=.0005):
        return dict(source=dict(owner=dict(physicalComponent=component)),
                    target=dict(axis=np.array([[0., 0., 0.], [sign*.002, -.004, 0.],
                                               [sign*.003, -.01, 0.]])),
                    mouthRadius_m=radius)

    def test_shared_guide_point_becomes_two_separated_finite_mouths(self):
        rows = [self.match(-1, 0), self.match(1, 1)]
        distal = [r['target']['axis'][-1].copy() for r in rows]
        report = separate_guide_mouths(rows)
        self.assertEqual(len(report), 1)
        self.assertGreaterEqual(report[0]['mouthCentreDistance_m'], .0012)
        self.assertLessEqual(report[0]['guideAdvance_m'], .004)
        for row, end in zip(rows, distal):
            np.testing.assert_array_equal(row['target']['axis'][-1], end)
            self.assertEqual(row['mouthRadius_m'], .0005)

    def test_does_not_move_distinct_source_ports(self):
        rows = [self.match(-1, 0), self.match(1, 1)]
        rows[1]['target']['axis'][:, 2] += .002
        original = rows[0]['target']['axis'].copy()
        self.assertEqual(separate_guide_mouths(rows), [])
        np.testing.assert_array_equal(rows[0]['target']['axis'], original)

    def test_nearby_distinct_points_still_require_finite_mouth_clearance(self):
        rows = [self.match(-1, 0), self.match(1, 1)]
        rows[1]['target']['axis'][:, 2] += .00015
        report = separate_guide_mouths(rows)
        self.assertEqual(len(report), 1)
        self.assertGreaterEqual(report[0]['mouthCentreDistance_m'], .0012)
        self.assertGreater(report[0]['guideAdvance_m'], 0.)

    def test_rejects_insufficient_fork_length(self):
        rows = [self.match(-1, 0, radius=.003), self.match(1, 1, radius=.003)]
        with self.assertRaisesRegex(ValueError, 'trim limit'):
            separate_guide_mouths(rows, maximum=.004)

    def test_trim_uses_arclength_and_keeps_distal_points(self):
        axis = np.array([[0., 0., 0.], [0., 3., 0.], [4., 3., 0.]])
        np.testing.assert_allclose(trim_course(axis, 5), [[2., 3., 0.], [4., 3., 0.]])
        with self.assertRaises(ValueError):
            trim_course(axis, 7)

    def test_source_fork_pose_is_exact_at_proximal_boundary(self):
        old = np.c_[np.zeros(100), np.linspace(0, 1, 100), np.zeros(100)]
        new = old+np.array([.1, .02, 0])
        tangent = np.tile([0., 1., 0.], (100, 1))
        posed, directions = retain_proximal_course(old, new, tangent, tangent)
        np.testing.assert_array_equal(posed[0], old[0])
        np.testing.assert_array_equal(directions[0], tangent[0])
        np.testing.assert_allclose(posed[40:], new[40:])
        self.assertTrue(np.isfinite(directions).all())


if __name__ == '__main__':
    unittest.main()
