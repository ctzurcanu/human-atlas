import unittest

import numpy as np

from geometry_tools.algorithms.declared_shared_disk import remove_shared_disk
from geometry_tools.algorithms.tube_continuity import continuity_audit


class DeclaredSharedDiskTests(unittest.TestCase):
    def fixture(self):
        angle = np.arange(6) * np.pi / 3
        points = np.concatenate([
            np.c_[np.cos(angle), np.sin(angle), np.full(6, z)]
            for z in [-1, 0, 1]
        ]) * .001
        faces = []
        for lower, upper in [(0, 6), (6, 12)]:
            for i in range(6):
                j = (i + 1) % 6
                faces.extend([[lower+i, lower+j, upper+i],
                              [lower+j, upper+j, upper+i]])
        faces.extend([[6, 6+i, 7+i] for i in range(1, 5)])
        return points, np.array(faces)

    def test_removal_opens_profile_without_changing_source_walls(self):
        points, faces = self.fixture()
        original = points.copy()
        keep, proof = remove_shared_disk(points, faces, range(6, 12), range(24, 28))
        np.testing.assert_array_equal(points, original)
        np.testing.assert_array_equal(keep, np.arange(24))
        audit = continuity_audit(points, faces[keep], expected_components=1,
                                 expected_mouths=2, expected_euler=0)
        self.assertTrue(audit['topological_continuity_pass'])
        self.assertTrue(proof['allProfileVertexLinksManifold'])
        self.assertTrue(proof['externalBoundaryEdgesLiteral'])
        self.assertFalse(proof['sourceCoordinatesChanged'])

    def test_same_side_source_walls_are_rejected(self):
        points, faces = self.fixture()
        points[:6, 2] = .001
        with self.assertRaisesRegex(ValueError, 'opposite sides'):
            remove_shared_disk(points, faces, range(6, 12), range(24, 28))

    def test_partial_disk_or_wrong_explicit_root_is_rejected(self):
        points, faces = self.fixture()
        with self.assertRaises(ValueError):
            remove_shared_disk(points, faces, range(6, 12), range(24, 27))
        with self.assertRaises(ValueError):
            remove_shared_disk(points, faces, [6, 7, 8, 9, 10, 10], range(24, 28))


if __name__ == '__main__':
    unittest.main()
