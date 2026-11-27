import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.internal_port_disk import remove_internal_port_disk
from geometry_tools.algorithms.tube_continuity import continuity_audit
from geometry_tools.tests import test_source_collar_union as source_fixtures


class InternalPortDiskTests(unittest.TestCase):
    def source(self, same_side=False):
        m = source_fixtures.SourceCollarUnionTests().source()
        # Keep the first source branch, reflect a literal copy across its root.
        vertices = m.vertices[:18].copy()
        vertices = np.r_[vertices, vertices[6:]*([1.4, 1.4, 1] if same_side else [1, 1, -1])]
        first = m.faces[:24]
        mapping = np.r_[np.arange(6), np.arange(18, 30)]
        second = mapping[first[:, [0, 2, 1]]]
        disk = np.array([[0, i, i+1] for i in range(1, 5)])
        return Mesh(vertices, np.r_[first, second, disk], m.frame)

    def test_only_blocking_disk_is_removed(self):
        m = self.source()
        q, proof = remove_internal_port_disk(m)
        np.testing.assert_array_equal(q.vertices, m.vertices)
        np.testing.assert_array_equal(q.faces, m.faces[:48])
        self.assertEqual(proof['removedSourceFaceIds'], list(range(48, 52)))
        self.assertTrue(continuity_audit(q.vertices, q.faces)['topological_continuity_pass'])

    def test_a_terminal_disk_is_not_a_join(self):
        m = self.source()
        terminal = Mesh(m.vertices, np.r_[m.faces[:24], m.faces[-4:]], m.frame)
        with self.assertRaises(ValueError):
            remove_internal_port_disk(terminal)

    def test_same_side_branches_do_not_authorize_disk_removal(self):
        with self.assertRaisesRegex(ValueError, 'opposite sides'):
            remove_internal_port_disk(self.source(same_side=True))

    def test_a_partial_patch_is_not_a_full_port_disk(self):
        m = self.source()
        with self.assertRaises(ValueError):
            remove_internal_port_disk(Mesh(m.vertices, m.faces[:-1], m.frame))


if __name__ == '__main__':
    unittest.main()
