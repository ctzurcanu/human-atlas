"""Source-triangle displacement energy, including thin and obtuse facets."""
import unittest

import numpy as np

from geometry_tools.algorithms.coupled_source_displacement import _triangle_fem_blocks
from geometry_tools.model import Mesh
from geometry_tools.tests.test_coupled_source_displacement import solve, pair, declaration


def fan(center=(.35, .4), *, rotate=False):
    p = np.array([[0., 0., 0.], [2., 0., 0.], [2., 1., 0.], [0., 1., 0.],
                  [center[0], center[1], 0.]])
    # Prescribe a true affine vector displacement, never absolute smoothed q.
    field = p @ np.array([[.04, -.02, .02], [.03, .01, -.02], [0., 0., 0.]])
    if rotate:
        angle = .47
        rotation = np.array([[np.cos(angle), 0., np.sin(angle)], [0., 1., 0.],
                             [-np.sin(angle), 0., np.cos(angle)]])
        p = p @ rotation.T + [.03, 1.2, -.04]
        field = field @ rotation.T
    mesh = Mesh(p, np.array([[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]]), 'test:metres')
    return mesh, field


class CoupledSourceFEMTests(unittest.TestCase):
    def check_affine_fan(self, center=(.35, .4), *, rotate=False):
        mesh, field = fan(center, rotate=rotate)
        prescribed = {'patch': dict(rawIds=[1, 2, 3], targetPoints_m=mesh.vertices[[1, 2, 3]] + field[[1, 2, 3]])}
        actual, proof = solve({'patch': mesh}, [], prescribed=prescribed,
                               fixed={'patch': [0]}, energy='source-triangle-fem')
        np.testing.assert_allclose(actual['patch'], mesh.vertices + field, atol=2e-12, rtol=0)
        np.testing.assert_array_equal(actual['patch'][0], mesh.vertices[0])
        self.assertEqual(proof['energy'], 'source-triangle-fem')
        self.assertTrue(proof['sourceEdgeAdjacencyRetainedForComponents'])
        self.assertFalse(proof['triangleQualityAndCalibreVerified'])
        self.assertFalse(proof['anatomicalAcceptance'])
        return mesh, field, prescribed

    def test_irregular_flat_patch_has_affine_displacement_exactness(self):
        mesh, field, prescribed = self.check_affine_fan()
        graph, _ = solve({'patch': mesh}, [], prescribed=prescribed, fixed={'patch': [0]})
        self.assertGreater(np.linalg.norm(graph['patch'][4] - (mesh.vertices + field)[4]), .001)

    def test_thin_original_triangle_keeps_the_same_affine_displacement_solution(self):
        mesh, _, _ = self.check_affine_fan((.35, 1e-8))
        _, proof = _triangle_fem_blocks(mesh.vertices, mesh.faces)
        self.assertLess(proof['minimumSourceArea_m2'], 2e-8)

    def test_full_world_float64_facets_work_in_an_oblique_plane(self):
        self.check_affine_fan(rotate=True)

    def test_obtuse_triangle_local_operator_is_psd_and_matches_actual_gradient_energy(self):
        p = np.array([[0., 0., 0.], [2., 0., 0.], [.2, .1, 0.]])
        block, proof = _triangle_fem_blocks(p, np.array([[0, 1, 2]]))
        matrix = block[0]
        np.testing.assert_allclose(matrix, matrix.T, atol=1e-14, rtol=0)
        np.testing.assert_allclose(matrix.sum(1), 0., atol=2e-14, rtol=0)
        self.assertGreaterEqual(np.linalg.eigvalsh(matrix).min(), -1e-12)
        self.assertGreater(matrix[0, 1], 0)  # Negative cotangent edge weight is legal.
        self.assertGreater(proof['positiveOffDiagonalEntries'], 0)
        values = np.array([.2, -.3, .7])
        edges = p[1:, :2] - p[0, :2]
        gradient = np.linalg.solve(edges, values[1:] - values[0])
        independently_integrated_energy = abs(np.linalg.det(edges)) * .5 * (gradient @ gradient)
        self.assertAlmostEqual(float(values @ matrix @ values), independently_integrated_energy, places=13)
        self.assertTrue(proof['positiveSemidefiniteByAreaWeightedGradientGram'])

    def test_fem_keeps_exact_aliases_paired_offsets_and_fixed_boundaries(self):
        meshes, interfaces = pair(gap=2e-7, alias=True)
        shift = np.array([.0001, -.00002, .00003])
        prescribed = {'a': declaration(meshes['a'], np.arange(18, 24), shift)}
        before = {pid: mesh.vertices.copy() for pid, mesh in meshes.items()}
        actual, proof = solve(meshes, interfaces, prescribed=prescribed,
            fixed={'a': list(range(6)), 'b': list(range(18, 24))}, energy='source-triangle-fem')
        np.testing.assert_array_equal(actual['a'][18], actual['a'][24])
        np.testing.assert_array_equal(actual['a'][:6], before['a'][:6])
        np.testing.assert_array_equal(actual['b'][18:24], before['b'][18:24])
        np.testing.assert_allclose(actual['b'][:6] - actual['a'][18:24],
                                   before['b'][:6] - before['a'][18:24], atol=2e-18, rtol=0)
        for pid, mesh in meshes.items():
            np.testing.assert_array_equal(mesh.vertices, before[pid])
        self.assertTrue(proof['fixedSourceCoordinatesLiteral'])

    def test_fem_zero_displacement_keeps_original_unsmoothed_coordinates(self):
        meshes, interfaces = pair(gap=2e-7, alias=True)
        actual, _ = solve(meshes, interfaces, energy='source-triangle-fem')
        for pid, mesh in meshes.items():
            np.testing.assert_array_equal(actual[pid], mesh.vertices)

    def test_fem_does_not_soften_conflicting_hard_constraints(self):
        meshes, interfaces = pair(gap=2e-7)
        with self.assertRaisesRegex(ValueError, 'hard fixed'):
            solve(meshes, interfaces, fixed={'a': [18]},
                  prescribed={'b': declaration(meshes['b'], [0], np.array([1e-10, 0., 0.]))},
                  tolerance=1e-6, energy='source-triangle-fem')

    def test_degenerate_source_triangles_and_unknown_energy_are_rejected(self):
        mesh = Mesh(np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]),
                    np.array([[0, 1, 2]]), 'test:metres')
        with self.assertRaisesRegex(ValueError, 'nondegenerate original facets'):
            solve({'line': mesh}, [], energy='source-triangle-fem')
        meshes, interfaces = pair()
        with self.assertRaises(ValueError):
            solve(meshes, interfaces, energy='unspecified-smoothing')


if __name__ == '__main__':
    unittest.main()
