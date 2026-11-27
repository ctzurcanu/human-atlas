import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_normal_duplicates import resolve_source_normal_duplicates


class SourceNormalDuplicateTest(unittest.TestCase):
    def source(self):
        theta = np.arange(6) * np.pi / 3
        vertices = np.concatenate([np.c_[np.full(6, x), np.cos(theta), np.sin(theta)] for x in [0, 1, 2]])
        faces = []
        for ring in range(2):
            for i in range(6):
                j = (i + 1) % 6
                faces += [[ring * 6 + i, ring * 6 + j, (ring + 1) * 6 + j],
                          [ring * 6 + i, (ring + 1) * 6 + j, (ring + 1) * 6 + i]]
        normals = np.tile(np.c_[np.zeros(6), np.cos(theta), np.sin(theta)], (3, 1))
        normals = np.r_[normals, np.tile([1., 0., 0.], (12, 1))]
        return Mesh(np.r_[vertices, vertices[:12]], np.r_[np.asarray(faces), np.asarray(faces[:12]) + 18],
                    'fixture:metres', dict(normals=normals, normalSigns=np.ones(30, dtype='i1'),
                                           colors=np.tile([150, 60, 50], (30, 1)).astype('u1')))

    def test_existing_smooth_strip_is_retained_without_changing_fields(self):
        mesh = self.source()
        before = {k: v.copy() for k, v in mesh.attributes.items()}
        kept, proof = resolve_source_normal_duplicates(mesh)
        np.testing.assert_array_equal(kept, np.arange(24))
        self.assertEqual(proof['removedSourceFaceIds'], list(range(24, 36)))
        self.assertLess(proof['strips'][0]['sourceNormalDiscontinuityScores'][0],
                        proof['strips'][0]['sourceNormalDiscontinuityScores'][1])
        for key in before:
            np.testing.assert_array_equal(before[key], mesh.attributes[key])

    def test_a_color_seam_cannot_be_discarded_as_a_normal_duplicate(self):
        mesh = self.source()
        fields = {k: v.copy() for k, v in mesh.attributes.items()}
        fields['colors'][18:] = [40, 70, 210]
        changed = Mesh(mesh.vertices, mesh.faces, mesh.frame, fields)
        with self.assertRaises(AssertionError):
            resolve_source_normal_duplicates(changed)

    def test_opposite_winding_is_not_a_redundant_source_copy(self):
        mesh = self.source()
        faces = mesh.faces.copy()
        faces[24:] = faces[24:, [0, 2, 1]]
        changed = Mesh(mesh.vertices, faces, mesh.frame, mesh.attributes)
        with self.assertRaises(ValueError):
            resolve_source_normal_duplicates(changed)


if __name__ == '__main__':
    unittest.main()
