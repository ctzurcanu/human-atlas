import unittest
import numpy as np
import manifold3d as md

from geometry_tools.algorithms.finite_cell_partition import partition_closed_cell
from geometry_tools.audit import audit
from geometry_tools.model import Mesh


def arrays(solid):
    mesh = solid.to_mesh64()
    return np.asarray(mesh.vert_properties)[:, :3], np.asarray(mesh.tri_verts, int)


class FinitePartitionTests(unittest.TestCase):
    def test_partial_source_faces_and_reciprocal_cut_preserve_whole_volume(self):
        vertices, faces = arrays(md.Manifold.cube([.02, .03, .04], center=True))
        mask_vertices, mask_faces = arrays(md.Manifold.cube([.03, .04, .05], center=True).translate([.015, .002, 0]))
        result = partition_closed_cell(vertices, faces, mask_vertices, mask_faces)
        np.testing.assert_array_equal(result['positions'][:len(vertices)], vertices)
        evidence = result['evidence']
        self.assertAlmostEqual(sum(evidence['cell_volumes_m3'].values()), .02*.03*.04, places=18)
        for face_array in result['cell_faces'].values():
            check = audit(Mesh(result['positions'], face_array, 'fixture'), require_closed=True,
                          expected_components=1, expected_genus=[0])
            self.assertTrue(check['passed'], check['failures'])

    def test_partition_of_wall_does_not_cap_its_lumen(self):
        outer = md.Manifold.cube([.04, .04, .04], center=True)
        lumen = md.Manifold.cube([.02, .02, .06], center=True)
        vertices, faces = arrays(outer-lumen)
        mask_vertices, mask_faces = arrays(md.Manifold.cube([.03, .05, .05], center=True).translate([.015, 0, 0]))
        result = partition_closed_cell(vertices, faces, mask_vertices, mask_faces)
        cut = result['positions'][result['cut_faces']]
        self.assertTrue(np.all(abs(cut.mean(1)[:, 1]) >= .01-1e-14))
        self.assertAlmostEqual(sum(result['evidence']['cell_volumes_m3'].values()), .04**3-.02**2*.04, places=18)
        for face_array in result['cell_faces'].values():
            check = audit(Mesh(result['positions'], face_array, 'fixture'), require_closed=True,
                          expected_components=1, expected_genus=[0])
            self.assertTrue(check['passed'], check['failures'])

    def test_classifier_must_create_two_nonempty_cells(self):
        vertices, faces = arrays(md.Manifold.cube([1., 1., 1.]))
        mask_vertices, mask_faces = arrays(md.Manifold.cube([1., 1., 1.]).translate([3., 0, 0]))
        with self.assertRaisesRegex(ValueError, 'two finite nonempty cells'):
            partition_closed_cell(vertices, faces, mask_vertices, mask_faces)


if __name__ == '__main__':
    unittest.main()
