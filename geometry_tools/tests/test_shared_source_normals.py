"""Normal fields follow complete declared original-bank stars across slices."""
import unittest
from unittest.mock import patch

import numpy as np

from geometry_tools.algorithms.shared_source_normals import transport_shared_source_normals
from geometry_tools.algorithms.surface_normals import transport_vertex_words
from geometry_tools.model import Mesh


def source_prism():
    angles = np.arange(6) * np.pi / 3
    p = np.array([[np.cos(a), np.sin(a), float(z)] for z in range(3) for a in angles])
    faces = []
    for ring in range(2):
        for side in range(6):
            a, b = ring * 6 + side, ring * 6 + (side + 1) % 6
            faces.extend([[a, b, b + 6], [a, b + 6, a + 6]])
    normals = np.rint(np.tile(np.column_stack([np.cos(angles), np.sin(angles), np.zeros(6)]), (3, 1)) * 32767).astype('<i2')
    uv = np.array([[side / 6, ring / 2] for ring in range(3) for side in range(6)], dtype='<f4')
    return p, np.array(faces), dict(normals=normals, normalSigns=np.ones(18, dtype='i1'), uvs=uv)


def split_sources(*, duplicate_face=False):
    p, f, attrs = source_prism()
    record = dict(chunk=0, positions=0, normals=1000, normalSigns=1500, uvs=1600,
                  vertexCount=100, geometryOrigin=[0., 0., 0.], geometryScale=1.)
    ids_a, ids_b = np.arange(12), np.arange(6, 18)
    face_a, face_b = f[:12], f[12:]
    if duplicate_face:
        face_b = np.vstack([f[12:], np.roll(f[0], 1)])
        ids_b = np.unique(face_b)
    meshes = {}
    for pid, ids, faces in [('a', ids_a, face_a), ('b', ids_b, face_b)]:
        local = np.full(18, -1, dtype=int)
        local[ids] = np.arange(len(ids))
        meshes[pid] = Mesh(p[ids], local[faces], 'test:metres',
                           {key: value[ids] for key, value in attrs.items()},
                           dict(source_vertex_ids=(10 + 3 * ids).tolist(), record=record,
                                bank_path='/explicit/test/source.bank', bank_sha256='fixture-source'))
    return p, f, attrs, meshes


def sliced_positions(q, meshes):
    return {pid: q[(np.array(mesh.source['source_vertex_ids']) - 10) // 3] for pid, mesh in meshes.items()}


def canonical_faces(faces):
    begin = np.argmin(faces, axis=1)
    return np.unique(np.take_along_axis(faces, (begin[:, None] + np.arange(3)) % 3, axis=1), axis=0)


class SharedSourceNormalTests(unittest.TestCase):
    def test_split_curved_prism_matches_complete_union_and_repairs_shared_word_conflict(self):
        p, f, attrs, meshes = split_sources()
        q = p.copy()
        q[:, 0] *= 1 + .1 * p[:, 2]
        q[:, 1] *= 1 - .06 * p[:, 2]
        q[:, 2] += .05 * p[:, 0] * p[:, 1]
        positions = sliced_positions(q, meshes)
        expected, _ = transport_vertex_words(p, q, canonical_faces(f), attrs['normals'], maximum_facet_rotation_degrees=65)
        independent = {pid: transport_vertex_words(mesh.vertices, positions[pid], mesh.faces,
                       mesh.attributes['normals'], maximum_facet_rotation_degrees=65)[0]
                       for pid, mesh in meshes.items()}
        self.assertTrue(np.any(independent['a'][6:12] != independent['b'][:6]))
        actual, proof = transport_shared_source_normals(meshes, positions, [['a', 'b']])
        for pid, mesh in meshes.items():
            ids = (np.array(mesh.source['source_vertex_ids']) - 10) // 3
            np.testing.assert_array_equal(actual[pid], expected[ids])
        np.testing.assert_array_equal(actual['a'][6:12], actual['b'][:6])
        self.assertEqual(proof['families'][0]['normalComputationFaces'], 24)
        self.assertFalse(proof['anatomicalAcceptance'])

    def test_oriented_duplicate_face_is_counted_once_only_for_normal_computation(self):
        p, f, attrs, meshes = split_sources(duplicate_face=True)
        q = p.copy()
        q[12:, 0] += .2
        expected, _ = transport_vertex_words(p, q, canonical_faces(f), attrs['normals'], maximum_facet_rotation_degrees=65)
        before = {pid: mesh.faces.copy() for pid, mesh in meshes.items()}
        actual, proof = transport_shared_source_normals(meshes, sliced_positions(q, meshes), [['a', 'b']])
        for pid, mesh in meshes.items():
            ids = (np.array(mesh.source['source_vertex_ids']) - 10) // 3
            np.testing.assert_array_equal(actual[pid], expected[ids])
            np.testing.assert_array_equal(mesh.faces, before[pid])
        self.assertEqual(proof['families'][0]['orientedDuplicateFacesCountedOnce'], 1)

    def test_unchanged_complete_source_stars_retain_literal_words(self):
        p, _, attrs, meshes = split_sources()
        q = p.copy()
        q[12:, 0] += .2
        before = {pid: {name: value.copy() for name, value in mesh.attributes.items()}
                  for pid, mesh in meshes.items()}
        actual, proof = transport_shared_source_normals(meshes, sliced_positions(q, meshes), [['a', 'b']])
        np.testing.assert_array_equal(actual['a'][:6], attrs['normals'][:6])
        for pid, mesh in meshes.items():
            for name, values in before[pid].items():
                np.testing.assert_array_equal(mesh.attributes[name], values)
        self.assertTrue(proof['unchangedCompleteSourceStarsLiteral'])

    def test_shared_pose_or_source_words_signs_and_uv_conflicts_are_rejected(self):
        p, _, _, meshes = split_sources()
        positions = sliced_positions(p, meshes)
        positions['b'] = positions['b'].copy()
        positions['b'][0, 0] += 1e-12
        with self.assertRaisesRegex(ValueError, 'new coordinates'):
            transport_shared_source_normals(meshes, positions, [['a', 'b']])
        for key in ['normals', 'normalSigns', 'uvs']:
            _, _, _, meshes = split_sources()
            b = meshes['b']
            modified = {name: value.copy() for name, value in b.attributes.items()}
            if key == 'normalSigns':
                modified[key][0] = -1
            elif key == 'normals':
                modified[key][0, 1] += 1
            else:
                modified[key][0, 0] += .01
            meshes['b'] = Mesh(b.vertices, b.faces, b.frame, modified, b.source)
            with self.subTest(key=key), self.assertRaisesRegex(ValueError, 'conflicting'):
                transport_shared_source_normals(meshes, sliced_positions(p, meshes), [['a', 'b']])

    def test_unrelated_coincident_families_are_transported_separately(self):
        p, f, attrs = source_prism()
        meshes = {'a': Mesh(p, f, 'test:metres', attrs), 'b': Mesh(p, f, 'test:metres', attrs)}
        qa, qb = p.copy(), p.copy()
        qa[12:, 0] += .2
        qb[12:, 1] -= .15
        actual, proof = transport_shared_source_normals(meshes, {'a': qa, 'b': qb}, [['a'], ['b']])
        for pid, q in [('a', qa), ('b', qb)]:
            expected, _ = transport_vertex_words(p, q, canonical_faces(f), attrs['normals'], maximum_facet_rotation_degrees=65)
            np.testing.assert_array_equal(actual[pid], expected)
        self.assertTrue(np.any(actual['a'] != actual['b']))
        self.assertEqual(len(proof['families']), 2)

    def test_mismatched_layout_missing_shared_identity_and_invalid_family_scope_reject(self):
        p, _, _, meshes = split_sources()
        b = meshes['b']
        source = dict(b.source, record=dict(b.source['record'], normals=1001))
        meshes['b'] = Mesh(b.vertices, b.faces, b.frame, b.attributes, source)
        with self.assertRaisesRegex(ValueError, 'bank/layout'):
            transport_shared_source_normals(meshes, sliced_positions(p, meshes), [['a', 'b']])
        _, _, _, meshes = split_sources()
        source = dict(meshes['b'].source)
        del source['source_vertex_ids']
        b = meshes['b']
        meshes['b'] = Mesh(b.vertices, b.faces, b.frame, b.attributes, source)
        positions = {'a': meshes['a'].vertices.copy(), 'b': meshes['b'].vertices.copy()}
        with self.assertRaisesRegex(ValueError, 'source_vertex_ids'):
            transport_shared_source_normals(meshes, positions, [['a', 'b']])
        with self.assertRaisesRegex(ValueError, 'singly owned'):
            transport_shared_source_normals(meshes, positions, [['a'], ['a']])

    def test_native_deadline_never_emits_partial_normals(self):
        _, _, _, meshes = split_sources()
        positions = {pid: mesh.vertices.copy() for pid, mesh in meshes.items()}
        with patch('geometry_tools.algorithms.shared_source_normals.time.monotonic', side_effect=[0., 10.]):
            with self.assertRaises(TimeoutError):
                transport_shared_source_normals(meshes, positions, [['a', 'b']])


if __name__ == '__main__':
    unittest.main()
