import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.declared_internal_mouth import weld_internal_mouth


def split_tube():
    theta = np.arange(6)*np.pi/3
    circle = np.c_[.001*np.cos(theta), .001*np.sin(theta), np.zeros(6)]
    points = np.concatenate([circle+[0, 0, z] for z in [0, .01, .010005, .02]])
    faces = []
    for a, b in [(0, 6), (12, 18)]:
        for i in range(6):
            j = (i+1)%6
            faces += [[a+i, a+j, b+i], [a+j, b+j, b+i]]
    words = np.tile(np.rint(np.c_[np.cos(theta), np.sin(theta), np.zeros(6)]*32767).astype('<i2'), (4, 1))
    return Mesh(points, np.asarray(faces), 'test:metres', dict(normals=words, normalSigns=np.ones(24, dtype='i1')))


class InternalMouthTests(unittest.TestCase):
    def test_complete_seam_joins_components_without_caps_or_removed_faces(self):
        old = split_tube()
        new, proof = weld_internal_mouth(old, np.arange(6, 12), np.arange(12, 18))
        np.testing.assert_array_equal(old.faces, new.faces)
        self.assertEqual(proof['lumenDisksCreated'], 0)
        self.assertEqual(proof['topologyBefore']['components'], 2)
        self.assertEqual(proof['topologyAfter']['components'], 1)
        self.assertEqual(proof['topologyAfter']['mouths'], 2)
        self.assertTrue(proof['topologyAfter']['topological_continuity_pass'])
        self.assertLess(proof['maximumMotion_m'], .0000025001)
        self.assertLess(max(proof['relativeSectionAreaChanges']), 1e-12)
        np.testing.assert_array_equal(old.vertices[np.r_[0:6, 18:24]], new.vertices[np.r_[0:6, 18:24]])

    def test_same_facing_source_rims_are_not_joined(self):
        with self.assertRaisesRegex(ValueError, 'opposed'):
            weld_internal_mouth(split_tube(), np.arange(0, 6), np.arange(12, 18), maximum_gap_m=.02)

    def test_large_unmatched_gap_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'tolerance'):
            weld_internal_mouth(split_tube(), np.arange(6, 12), np.arange(12, 18), maximum_gap_m=.000001)


if __name__ == '__main__':
    unittest.main()
