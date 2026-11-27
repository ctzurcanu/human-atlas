import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.terminal_profile_repair import reconstructed_terminal_profile
from geometry_tools.open3d_backend import diagnostics


class TerminalProfileReconstruction(unittest.TestCase):
    def source(self, minor=1e-7):
        angle = np.arange(6) * np.pi / 3
        first = np.c_[.0003 * np.cos(angle), .0003 * np.sin(angle), np.zeros(6)]
        second = np.c_[.00008 * np.cos(angle), minor * np.sin(angle), np.full(6, .003)]
        faces = []
        for i in range(6):
            j = (i + 1) % 6
            faces += [[i, j, 6+j], [i, 6+j, 6+i]]
        return Mesh(np.r_[first, second], np.array(faces), 'fixture:metres')

    def test_collapsed_tip_uses_neighbor_shape_at_same_centre_and_diameter(self):
        source = self.source()
        vertices, proof = reconstructed_terminal_profile(source, list(range(6, 12)))
        np.testing.assert_array_equal(vertices[:6], source.vertices[:6])
        np.testing.assert_allclose(vertices[6:].mean(0), source.vertices[6:].mean(0), atol=5e-16, rtol=0)
        diameter = lambda p: np.linalg.norm(p[:, None] - p[None], axis=2).max()
        self.assertAlmostEqual(diameter(vertices[6:]), diameter(source.vertices[6:]), places=15)
        self.assertTrue(proof['rebuiltProfile']['convex'])
        self.assertFalse(proof['originalInvalidAreaPreserved'])
        self.assertEqual(diagnostics(Mesh(vertices, source.faces, source.frame))['self_intersection_pairs'], 0)

    def test_valid_round_source_profile_cannot_be_reconstructed(self):
        with self.assertRaises(ValueError):
            reconstructed_terminal_profile(self.source(.00008), list(range(6, 12)))


if __name__ == '__main__':
    unittest.main()
