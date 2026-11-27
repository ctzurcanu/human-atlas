"""Source boundary ownership and locality of the shared-mouth field."""
import unittest

import numpy as np

from geometry_tools.algorithms.coupled_mouth_field import mouth_displacement_basis
from geometry_tools.model import Mesh


def open_tube(*, rings=6, sides=6):
    angle = np.arange(sides) * (2 * np.pi / sides)
    points = np.array([
        [.001 * np.cos(a), .001 * np.sin(a), ring * .001]
        for ring in range(rings) for a in angle
    ])
    faces = []
    for ring in range(rings - 1):
        for side in range(sides):
            a = ring * sides + side
            b = ring * sides + (side + 1) % sides
            faces.extend([[a, b, b + sides], [a, b + sides, a + sides]])
    return Mesh(points, np.array(faces), 'test:metres')


def mouth(mesh, raw_ids, displacement):
    return dict(sourceRawIds=list(raw_ids),
                targetPoints_m=mesh.vertices[raw_ids] + displacement)


class CoupledMouthFieldTests(unittest.TestCase):
    def test_actual_opening_moves_and_other_opening_and_far_surface_stay_fixed(self):
        mesh = open_tube()
        source = mesh.vertices.copy()
        shift = np.array([.0001, -.00003, .00002])
        field, proof = mouth_displacement_basis(
            mesh, [mouth(mesh, np.arange(6), shift)], support_m=.0025)

        np.testing.assert_allclose(field[:6], np.broadcast_to(shift, (6, 3)), atol=1e-18)
        np.testing.assert_array_equal(field[18:], np.zeros((18, 3)))
        self.assertTrue(np.any(np.linalg.norm(field[6:18], axis=1) > 0))
        # Positive source-edge weights interpolate without displacement overshoot.
        factors = field[6:18, 0] / shift[0]
        self.assertTrue(np.all((factors > 0) & (factors < 1)))
        np.testing.assert_allclose(field[6:18], factors[:, None] * shift, atol=1e-18)
        np.testing.assert_array_equal(mesh.vertices, source)
        self.assertTrue(proof['sourceGraphWeightsFixed'])
        self.assertFalse(proof['anatomicalAcceptance'])

    def test_exact_coordinate_aliases_share_one_displacement_across_a_seam(self):
        source = open_tube()
        points = np.vstack([source.vertices, source.vertices[[0, 7]]])
        faces = source.faces.copy()
        # Split selected raw corners without changing the physical surface.
        faces[0, 0] = len(source.vertices)
        location = np.argwhere(faces == 7)[0]
        faces[tuple(location)] = len(source.vertices) + 1
        mesh = Mesh(points, faces, source.frame)
        root = np.array([len(source.vertices), 1, 2, 3, 4, 5])
        field, _ = mouth_displacement_basis(
            mesh, [mouth(mesh, root, np.array([.0001, 0, 0]))], support_m=.0025)
        np.testing.assert_array_equal(field[0], field[-2])
        np.testing.assert_array_equal(field[7], field[-1])

    def test_fields_are_additive_for_fixed_source_weights(self):
        mesh = open_tube()
        root = np.arange(6)
        a = np.array([.0001, -.00002, 0])
        b = np.array([-.00003, 0, .00004])
        fa, _ = mouth_displacement_basis(mesh, [mouth(mesh, root, a)], support_m=.0025)
        fb, _ = mouth_displacement_basis(mesh, [mouth(mesh, root, b)], support_m=.0025)
        combined, _ = mouth_displacement_basis(mesh, [mouth(mesh, root, a + b)], support_m=.0025)
        np.testing.assert_allclose(combined, fa + fb, atol=1e-18, rtol=1e-12)

    def test_multiple_disjoint_mouth_declarations_move_only_their_components(self):
        first = open_tube()
        second = open_tube()
        points = np.vstack([first.vertices, second.vertices + [.01, 0, 0]])
        faces = np.vstack([first.faces, second.faces + len(first.vertices)])
        mesh = Mesh(points, faces, first.frame)
        root_a = np.arange(6)
        root_b = root_a + len(first.vertices)
        a = np.array([.0001, 0, 0])
        b = np.array([0, .00005, 0])
        field, _ = mouth_displacement_basis(
            mesh, [mouth(mesh, root_a, a), mouth(mesh, root_b, b)], support_m=.0025)
        np.testing.assert_allclose(field[root_a], np.broadcast_to(a, (6, 3)), atol=1e-18)
        np.testing.assert_allclose(field[root_b], np.broadcast_to(b, (6, 3)), atol=1e-18)
        np.testing.assert_array_equal(field[30:36], np.zeros((6, 3)))
        np.testing.assert_array_equal(field[66:72], np.zeros((6, 3)))
        np.testing.assert_array_equal(field[:36, 1:], np.zeros((36, 2)))
        np.testing.assert_array_equal(field[36:, [0, 2]], np.zeros((36, 2)))

    def test_incomplete_opening_and_internal_ring_are_not_anatomical_mouths(self):
        mesh = open_tube()
        for root in [np.arange(5), np.arange(6, 12)]:
            with self.subTest(root=root.tolist()), self.assertRaisesRegex(ValueError, 'complete actual source mouth'):
                mouth_displacement_basis(
                    mesh, [mouth(mesh, root, np.array([.0001, 0, 0]))], support_m=.0025)

    def test_two_declarations_cannot_own_the_same_physical_mouth(self):
        mesh = open_tube()
        declaration = mouth(mesh, np.arange(6), np.array([.0001, 0, 0]))
        with self.assertRaisesRegex(ValueError, 'share physical ownership'):
            mouth_displacement_basis(mesh, [declaration, declaration], support_m=.0025)

    def test_all_actual_openings_may_be_explicitly_owned(self):
        mesh = open_tube()
        start, end = np.arange(6), np.arange(30, 36)
        a, b = np.array([.0001, 0, 0]), np.array([0, .00005, 0])
        field, proof = mouth_displacement_basis(
            mesh, [mouth(mesh, start, a), mouth(mesh, end, b)], support_m=.0025)
        np.testing.assert_allclose(field[start], np.broadcast_to(a, (6, 3)), atol=1e-18)
        np.testing.assert_allclose(field[end], np.broadcast_to(b, (6, 3)), atol=1e-18)
        self.assertEqual(proof['fixedOtherBoundaryPhysicalIds'], [])


if __name__ == '__main__':
    unittest.main()
