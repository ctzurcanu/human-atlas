import unittest
import numpy as np
from geometry_tools.algorithms.source_profile_alignment import proper_profile_rotation


class SourceProfileAlignmentTests(unittest.TestCase):
    def test_proper_alignment_retains_an_elliptical_source_section(self):
        t = np.arange(6)*np.pi/3
        source = np.c_[.002*np.cos(t), .001*np.sin(t), np.zeros(6)]
        angle = .45
        expected = np.array([[1., 0., 0.], [0., np.cos(angle), -np.sin(angle)], [0., np.sin(angle), np.cos(angle)]])
        target = source@expected.T+[.01, -.02, .03]
        rotation = proper_profile_rotation(source, target)
        np.testing.assert_allclose(rotation, expected, atol=1e-14)
        self.assertAlmostEqual(np.linalg.det(rotation), 1.)
        fitted = (source-source.mean(0))@rotation.T+source.mean(0)
        np.testing.assert_allclose(np.linalg.norm(fitted[:, None]-fitted[None], axis=2),
                                   np.linalg.norm(source[:, None]-source[None], axis=2), atol=1e-17)
        np.testing.assert_allclose(fitted.mean(0), source.mean(0), atol=1e-17)

    def test_collinear_sections_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'plane'):
            proper_profile_rotation(np.c_[np.arange(6), np.zeros((6, 2))],
                                    np.c_[np.arange(6), np.zeros((6, 2))])


if __name__ == '__main__':
    unittest.main()
