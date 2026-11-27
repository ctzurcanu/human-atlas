import unittest

import numpy as np

from geometry_tools.algorithms.source_profile_axis import source_profile_graph, primary_axis


def tube(offset=0.):
    angle = np.arange(6)*np.pi/3
    y = np.linspace(1, 0, 7)
    v = np.array([[offset+.05*np.cos(a)+.01*s*s, s, .04*np.sin(a)]
                  for s in y for a in angle])
    faces = []
    for ring in range(6):
        for j in range(6):
            a, b = ring*6+j, ring*6+(j+1)%6
            c, d = a+6, b+6
            faces.extend([[a, c, d], [a, d, b]])
    return v, np.array(faces)


class SourceProfileAxisTests(unittest.TestCase):
    def test_explicit_exporter_pairs_survive_reordered_triangles(self):
        v, f = tube()
        order = np.roll(np.arange(len(f)), 1)
        reverse = np.argsort(order)
        pairing = reverse[np.arange(len(f)).reshape(-1, 2)]
        moved_faces = f[order]
        before = v.copy()
        guide = source_profile_graph(v, moved_faces, triangle_pairing=pairing)
        self.assertEqual(len(guide['profiles']), 7)
        axis, proof = primary_axis(guide, moved_faces)
        np.testing.assert_allclose(axis[0], v[:6].mean(0), atol=1e-12)
        np.testing.assert_allclose(axis[-1], v[-6:].mean(0), atol=1e-12)
        self.assertGreater(proof['courseLength_m'], 1)
        np.testing.assert_array_equal(v, before)

    def test_nearby_disconnected_tube_is_never_an_attachment(self):
        a, fa = tube()
        b, fb = tube(.101)
        vertices = np.r_[a, b]
        faces = np.r_[fa, fb+len(a)]
        guide = source_profile_graph(vertices, faces,
                                     triangle_pairing=np.arange(len(faces)).reshape(-1, 2))
        axis, _ = primary_axis(guide, fa)
        self.assertLess(axis[:, 0].max(), .011)
        self.assertEqual(len(guide['boundary_nodes']), 4)

    def test_reused_exporter_face_is_rejected(self):
        v, f = tube()
        pairs = np.arange(len(f)).reshape(-1, 2)
        pairs[-1] = pairs[0]
        with self.assertRaisesRegex(ValueError, 'multiple exporter quads'):
            source_profile_graph(v, f, triangle_pairing=pairs)


if __name__ == '__main__':
    unittest.main()
