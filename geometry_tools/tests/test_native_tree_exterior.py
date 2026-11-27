"""Independent geometric obligations for an explicitly anchored branching cast."""
import unittest

import numpy as np

from geometry_tools.algorithms.native_tree_exterior import (
    NativeTreeExteriorError,
    _cone_fields,
    construct,
)
from geometry_tools.model import Mesh


def branching_fixture():
    vertices, normals, uvs, faces, groups, mouths = [], [], [], [], [], []

    def tube(start, end, radius, internal_first=False):
        start, end = np.asarray(start, float), np.asarray(end, float)
        direction = end - start
        direction /= np.linalg.norm(direction)
        helper = np.array([0., 1., 0.])
        first = np.cross(helper, direction)
        first /= np.linalg.norm(first)
        second = np.cross(direction, first)
        offset, first_face = len(vertices), len(faces)
        count = 12
        for station, center in enumerate([start, end]):
            controls = []
            for index in range(count):
                angle = 2 * np.pi * index / count
                normal = np.cos(angle) * first + np.sin(angle) * second
                controls.append(len(vertices))
                vertices.append(center + radius * normal)
                normals.append(np.rint(normal * 32767).astype(np.int16))
                uvs.append([index / count, station])
            mouths.append(dict(kind='internal-source-junction' if station == 0 and internal_first
                               else 'original-exterior-mouth',
                               positions_m=[vertices[index].tolist() for index in controls],
                               source_control_ids=controls))
        for index in range(count):
            following = (index + 1) % count
            a, b, c, d = offset + index, offset + following, offset + count + index, offset + count + following
            faces.extend([[a, b, d], [a, d, c]])
        groups.append(np.arange(first_face, len(faces)))

    tube([0, 0, 0], [0, 0, .02], .004)
    # The source branch begins outside the parent.  A distinct, exposed sleeve
    # must connect it; merely deleting its start disk cannot pass this fixture.
    tube([.005, 0, .01], [.015, 0, .01], .002, internal_first=True)
    mesh = Mesh(np.asarray(vertices), np.asarray(faces, np.int64), 'analytical-fixture:metres',
                {'normals': np.asarray(normals, np.int16), 'uvs': np.asarray(uvs, np.float32)})
    return mesh, groups, mouths, {2: [0., 0., .01]}


class NativeTreeExteriorTests(unittest.TestCase):
    def test_an_unmatched_declared_genus_fails_without_automatic_relaxation(self):
        source, groups, mouths, anchors = branching_fixture()
        _, proof = construct(source, groups, mouths, anchors, wall_seconds=10, expected_genus=1)
        self.assertFalse(proof['passed'])
        self.assertIn('wrong_declared_genus_topology', proof['failures'])
        self.assertEqual(proof['expected_genus'], 1)

    def test_corner_charts_survive_a_union_at_shared_geometry_controls(self):
        source, groups, mouths, anchors = branching_fixture()
        tri = source.triangles
        face_normals = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
        face_normals /= np.linalg.norm(face_normals, axis=1)[:, None]
        words = np.repeat(np.rint(face_normals*32767).astype(np.int16)[:, None], 3, axis=1)
        uv = source.attributes['uvs'][source.faces].copy()
        uv[:, :, 0] += np.arange(len(uv), dtype=np.float32)[:, None]
        arrays, proof = construct(source, groups, mouths, anchors, wall_seconds=10,
                                  source_corner_normal_words=words, source_corner_uvs=uv)
        self.assertTrue(proof['passed'], proof['failures'])
        witnessed = set()
        for fi, face in enumerate(arrays['faces']):
            parent = int(arrays['source_parent_face_ids'][fi])
            if parent < 0:
                continue
            for corner, point in enumerate(arrays['vertices'][face]):
                matches = np.flatnonzero(np.all(source.triangles[parent] == point, axis=1))
                if not len(matches):
                    continue
                original_corner = matches[0]
                np.testing.assert_array_equal(arrays['corner_normal_words'][fi, corner],
                                              words[parent, original_corner])
                np.testing.assert_array_equal(arrays['corner_uvs'][fi, corner],
                                              uv[parent, original_corner])
                witnessed.add(parent)
        self.assertGreater(len(witnessed), 12)

    def test_face_normal_mode_handles_an_oblique_cone_without_changing_uvs(self):
        angle = np.arange(12) * 2 * np.pi / 12
        rim = np.c_[np.cos(angle), np.sin(angle), np.zeros(12)]
        tip = np.array([2., 0., .1])
        faces = np.asarray([[(j + 1) % 12, j, 12] for j in range(12)])
        mean, uv_mean = _cone_fields(rim, tip, faces)
        flat, uv_flat = _cone_fields(rim, tip, faces, generated_normal_mode='face')
        triangles = np.vstack([rim, tip])[faces]
        normal = np.cross(triangles[:, 1] - triangles[:, 0],
                          triangles[:, 2] - triangles[:, 0])
        self.assertTrue(np.any(np.einsum('fki,fi->fk', mean, normal) <= 0))
        self.assertTrue(np.all(np.einsum('fki,fi->fk', flat, normal) > 0))
        np.testing.assert_array_equal(uv_flat, uv_mean)

    def test_exposed_sleeve_keeps_three_exact_ports_and_source_fields(self):
        source, groups, mouths, anchors = branching_fixture()
        before = (source.vertices.copy(), source.faces.copy(),
                  source.attributes['normals'].copy(), source.attributes['uvs'].copy())
        arrays, proof = construct(source, groups, mouths, anchors, wall_seconds=10)
        self.assertTrue(proof['passed'], proof['failures'])
        self.assertEqual(proof['expected_exterior_mouths'], 3)
        self.assertEqual(proof['topology']['open_edges'], 36)
        self.assertEqual(len(proof['topology']['components']), 1)
        self.assertEqual(proof['topology']['components'][0]['euler'], -1)
        self.assertGreater(proof['derived_fragments'], 0)
        self.assertEqual(proof['self_intersections']['defect_pairs'], 0)
        self.assertTrue(proof['self_intersections']['complete_within_pair_budget'])
        self.assertEqual(proof['temporary_classifier_faces_remaining'], 0)
        self.assertGreater(proof['removed_classifier_fragments'], 0)
        self.assertFalse(proof['anatomical_acceptance'])
        self.assertFalse(proof['finite_wall_or_physiological_patency_claim'])
        actual = set()
        for number, input_index in enumerate([0, 1, 3]):
            ids = arrays[f'external_boundary_{number}_vertex_ids']
            np.testing.assert_array_equal(arrays['vertices'][ids], mouths[input_index]['positions_m'])
            np.testing.assert_array_equal(arrays[f'external_boundary_{number}_source_control_ids'],
                                          mouths[input_index]['source_control_ids'])
            actual.update(tuple(sorted((int(ids[j]), int(ids[(j + 1) % len(ids)])))) for j in range(len(ids)))
        edge_counts = {}
        for a, b, c in arrays['faces']:
            for edge in [(a, b), (b, c), (c, a)]:
                key = tuple(sorted(map(int, edge)))
                edge_counts[key] = edge_counts.get(key, 0) + 1
        self.assertEqual(actual, {edge for edge, count in edge_counts.items() if count == 1})
        for face_id, face in enumerate(arrays['faces']):
            parent = int(arrays['source_parent_face_ids'][face_id])
            derived = int(arrays['derived_parent_face_ids'][face_id])
            self.assertEqual(parent >= 0, derived < 0)
            if parent < 0:
                self.assertGreaterEqual(derived, 0)
                chart = arrays['derived_vertices'][arrays['derived_faces'][derived]]
                weights = arrays['corner_derived_barycentric'][face_id]
                np.testing.assert_allclose(weights @ chart, arrays['vertices'][face], atol=2e-12, rtol=0)
                np.testing.assert_array_equal(arrays['corner_uvs'][face_id],
                                              (weights @ arrays['derived_corner_uvs'][derived]).astype(np.float32))
                continue
            np.testing.assert_allclose(arrays['corner_source_barycentric'][face_id] @ source.triangles[parent],
                                       arrays['vertices'][face], atol=2e-12, rtol=0)
            for corner, vertex in enumerate(face):
                point = arrays['vertices'][vertex]
                matches = np.flatnonzero(np.all(source.triangles[parent] == point, axis=1))
                if len(matches):
                    original = source.faces[parent, matches[0]]
                    np.testing.assert_array_equal(arrays['corner_normal_words'][face_id, corner],
                                                  source.attributes['normals'][original])
                    np.testing.assert_array_equal(arrays['corner_uvs'][face_id, corner],
                                                  source.attributes['uvs'][original])
        for actual_array, original in zip((source.vertices, source.faces, source.attributes['normals'],
                                           source.attributes['uvs']), before):
            np.testing.assert_array_equal(actual_array, original)

    def test_missing_anchor_or_false_source_binding_is_rejected(self):
        source, groups, mouths, anchors = branching_fixture()
        with self.assertRaisesRegex(NativeTreeExteriorError, 'Every internal mouth'):
            construct(source, groups, mouths, {}, wall_seconds=10)
        with self.assertRaisesRegex(NativeTreeExteriorError, 'original exterior mouth'):
            construct(source, groups, mouths, {**anchors, 0: [0, 0, .005]}, wall_seconds=10)
        corrupted = [dict(mouth) for mouth in mouths]
        corrupted[0]['positions_m'] = np.asarray(corrupted[0]['positions_m']) + [0, 0, .0001]
        with self.assertRaisesRegex(NativeTreeExteriorError, 'literal source controls'):
            construct(source, groups, corrupted, anchors, wall_seconds=10)


if __name__ == '__main__':
    unittest.main()
