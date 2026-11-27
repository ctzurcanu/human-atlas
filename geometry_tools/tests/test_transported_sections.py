import unittest
import numpy as np
from geometry_tools.commands.audit_transported_sections import section_polygons


def tube(ring):
    vertices = np.concatenate([np.c_[ring[:, 0], np.full(len(ring), y), ring[:, 1]]
                               for y in [-1., 1.]])
    n = len(ring)
    faces = []
    for a in range(n):
        b = (a + 1) % n
        faces.extend([[a, b, b + n], [a, b + n, a + n]])
    return vertices, np.array(faces)


class TransportedSectionTests(unittest.TestCase):
    def test_recut_after_rotation_uses_actual_plane(self):
        v, f = tube(np.array([[-.1, -.1], [.1, -.1], [.1, .1], [-.1, .1]]))
        a = .3
        R = np.array([[np.cos(a), -np.sin(a), 0], [np.sin(a), np.cos(a), 0], [0, 0, 1]])
        q = v @ R.T
        polygons, evidence = section_polygons(q, f, .23)
        self.assertTrue(evidence['closed'])
        self.assertFalse(evidence['invalidLoops'])
        self.assertAlmostEqual(polygons[0].area, .04 / np.cos(a), places=12)

    def test_crossed_surface_is_rejected(self):
        v, f = tube(np.array([[-1, -1], [1, 1], [-1, 1], [1, -1]]))
        _, evidence = section_polygons(v, f, .23)
        self.assertEqual(evidence['invalidLoops'], [0])

    def test_open_wall_is_not_silently_filled(self):
        v, f = tube(np.array([[-1, -1], [1, -1], [1, 1], [-1, 1]]))
        polygons, evidence = section_polygons(v, f[:-2], .23)
        self.assertFalse(evidence['closed'])
        self.assertEqual(len(polygons), 0)


if __name__ == '__main__':
    unittest.main()
