import unittest

import numpy as np

from geometry_tools.algorithms.compact_rbf_translation import CompactRBFTranslation


class CompactRBFTranslationTest(unittest.TestCase):
    def field(self):
        return CompactRBFTranslation(
            [[0, 0, 0], [.006, 0, 0]],
            [[0, .00002, 0], [.00001, -.00001, 0]],
            .01, 'test-metres')

    def test_interpolates_controls_preserves_shared_points_and_fixed_exterior(self):
        field = self.field()
        p = np.array([[0, 0, 0], [.006, 0, 0], [.003, .001, 0],
                      [.003, .001, 0], [.1, 0, 0]])
        q = field.map(p)
        np.testing.assert_allclose(q[:2], [[0, .00002, 0], [.00601, -.00001, 0]],
                                   rtol=0, atol=1e-14)
        np.testing.assert_array_equal(q[2], q[3])
        np.testing.assert_array_equal(q[4], p[4])
        self.assertLess(field.evidence()['sumDisplacementLipschitzBound'], .1)

    def test_jacobian_and_transported_normal_match_mapped_tangent_plane(self):
        field = self.field()
        p = np.array([[.003, .002, .001]])
        h = 1e-7
        finite = np.column_stack([(field.map(p + np.eye(3)[i]*h) -
                                   field.map(p - np.eye(3)[i]*h))[0]/(2*h)
                                  for i in range(3)])
        np.testing.assert_allclose(field.jacobians(p)[0], finite, atol=1e-9)
        normal = field.transform_normals(p, [[0, 0, 1]])[0]
        self.assertAlmostEqual(np.linalg.norm(normal), 1)
        self.assertLess(abs(normal @ finite[:, 0]), 1e-9)
        self.assertLess(abs(normal @ finite[:, 1]), 1e-9)

    def test_rejects_ambiguous_controls_and_excess_distortion(self):
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            CompactRBFTranslation([[0, 0, 0]]*2, [[0, 0, 0]]*2, .01, 'test')
        with self.assertRaisesRegex(ValueError, 'distortion'):
            CompactRBFTranslation([[0, 0, 0]], [[.001, 0, 0]], .01, 'test',
                                  maximum_lipschitz=.1)


if __name__ == '__main__':
    unittest.main()
