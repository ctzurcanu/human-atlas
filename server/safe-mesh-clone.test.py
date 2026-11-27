"""A posed mesh must never retain the source's already-computed witnesses."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from geometry_tools.algorithms.safe_mesh_clone import clone_with_geometry
from geometry_tools.model import Mesh


def fixture():
    return Mesh(
        np.array([[0., 0., 0.], [.01, 0., 0.], [0., .01, 0.], [0., 0., .01]]),
        np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]]),
        'fixture:metres',
        {'normals': np.tile(np.array([0, 0, 32767], dtype='<i2'), (4, 1)),
         'normalSigns': np.array([1, -1, 1, -1], dtype='i1'),
         'uvs': np.array([[0., 0.], [1., 0.], [0., 1.], [1., 1.]], dtype='<f4')},
        {'record': {'id': 'source-fixture'}, 'source_vertex_ids': [11, 13, 17, 19]})


class SafeMeshClone(unittest.TestCase):
    def test_primed_source_caches_are_not_reused_after_geometry_and_normal_changes(self):
        source = fixture()
        old_geometry, old_attributes, old_surface = source.geometry_sha, source.attributes_sha, source.surface_sha
        old_triangles, old_bounds = source.triangles.copy(), source.bounds.copy()
        points = source.vertices.copy()
        points[3] += [.003, .002, .004]
        normals = np.tile(np.array([0, 32767, 0], dtype='<i2'), (4, 1))

        posed = clone_with_geometry(source, points, attribute_updates={'normals': normals})

        self.assertTrue({'geometry_sha', 'attributes_sha', 'surface_sha', 'triangles', 'bounds'}.isdisjoint(posed.__dict__))
        self.assertNotEqual(posed.geometry_sha, old_geometry)
        self.assertNotEqual(posed.attributes_sha, old_attributes)
        self.assertNotEqual(posed.surface_sha, old_surface)
        np.testing.assert_array_equal(posed.triangles, points[source.faces])
        np.testing.assert_array_equal(posed.bounds, [points.min(0), points.max(0)])
        self.assertEqual(source.geometry_sha, old_geometry)
        self.assertEqual(source.attributes_sha, old_attributes)
        np.testing.assert_array_equal(source.triangles, old_triangles)
        np.testing.assert_array_equal(source.bounds, old_bounds)
        np.testing.assert_array_equal(source.attributes['normals'], [[0, 0, 32767]]*4)

    def test_input_banks_and_provenance_are_independent_and_stored_arrays_are_frozen(self):
        source = fixture()
        points = source.vertices.copy()
        normals = source.attributes['normals'].copy()
        posed = clone_with_geometry(source, points, attribute_updates={'normals': normals})
        points[:] = 9.
        normals[:] = 0

        np.testing.assert_array_equal(posed.vertices, source.vertices)
        np.testing.assert_array_equal(posed.attributes['normals'], source.attributes['normals'])
        for original, copied in [(source.vertices, posed.vertices), (source.faces, posed.faces),
                                 *[(source.attributes[k], posed.attributes[k]) for k in source.attributes]]:
            self.assertFalse(np.shares_memory(original, copied))
            self.assertFalse(copied.flags.writeable)
            with self.assertRaises(ValueError):
                copied.flat[0] = 0
        posed.source['record']['id'] = 'clone-only'
        posed.source['source_vertex_ids'][0] = 999
        self.assertEqual(source.source['record']['id'], 'source-fixture')
        self.assertEqual(source.source['source_vertex_ids'][0], 11)

    def test_attribute_only_change_and_fresh_npz_roundtrip_have_correct_hashes(self):
        source = fixture()
        original_geometry, original_attributes = source.geometry_sha, source.attributes_sha
        normals = np.tile(np.array([32767, 0, 0], dtype='<i2'), (4, 1))
        posed = clone_with_geometry(source, attribute_updates={'normals': normals})
        self.assertEqual(posed.geometry_sha, original_geometry)
        self.assertNotEqual(posed.attributes_sha, original_attributes)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'posed.npz'
            posed.save(path)
            loaded = Mesh.load(path)
        self.assertEqual(loaded.geometry_sha, posed.geometry_sha)
        self.assertEqual(loaded.attributes_sha, posed.attributes_sha)
        self.assertEqual(loaded.surface_sha, posed.surface_sha)
        np.testing.assert_array_equal(loaded.vertices, posed.vertices)
        np.testing.assert_array_equal(loaded.faces, source.faces)
        for key in posed.attributes:
            np.testing.assert_array_equal(loaded.attributes[key], posed.attributes[key])
            self.assertFalse(loaded.attributes[key].flags.writeable)
        self.assertEqual(loaded.source['source_vertex_ids'], source.source['source_vertex_ids'])

    def test_vertex_count_changes_require_explicit_new_ancestry(self):
        source = fixture()
        with self.assertRaisesRegex(ValueError, 'source vertex-bank shape'):
            clone_with_geometry(source, source.vertices[:3])


if __name__ == '__main__':
    unittest.main()
