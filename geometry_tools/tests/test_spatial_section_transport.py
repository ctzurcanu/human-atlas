import unittest
import numpy as np
from geometry_tools.algorithms.tube_sections import transport_spatial_sections, TubeGeometryError


class SpatialSectionTransportTests(unittest.TestCase):
    def source(self):
        angles = np.arange(12) * np.pi / 6
        contour = np.c_[.01 * np.cos(angles), .008 * np.sin(angles), .002 * np.cos(2 * angles)]
        return np.array([[contour + [0, 0, z] for z in [0, .02, .04, .06]]])

    def test_nonplanar_shape_and_proper_orientation_are_retained(self):
        p = self.source()
        target = np.array([[0, 0, 0], [.003, 0, .02], [.008, .002, .04], [.015, .006, .06]])
        q, rotations, report = transport_spatial_sections(p, target)
        self.assertTrue(report['sectionwise_isometry_pass'])
        self.assertFalse(report['profiles_flattened'])
        self.assertGreater(report['maximum_source_axial_extent_m'], .003)
        np.testing.assert_allclose(q[0].mean(1), target, atol=1e-16)
        np.testing.assert_allclose(np.linalg.det(rotations), 1., atol=1e-14)
        before = np.linalg.norm(p[..., :, None, :] - p[..., None, :, :], axis=-1)
        after = np.linalg.norm(q[..., :, None, :] - q[..., None, :, :], axis=-1)
        np.testing.assert_allclose(after, before, atol=1e-16)

    def test_invalid_target_shape_is_rejected(self):
        with self.assertRaises(TubeGeometryError):
            transport_spatial_sections(self.source(), np.zeros((3, 3)))

    def test_duplicate_stations_are_rejected(self):
        with self.assertRaises(TubeGeometryError):
            transport_spatial_sections(self.source(), np.zeros((4, 3)))


if __name__ == '__main__':
    unittest.main()
