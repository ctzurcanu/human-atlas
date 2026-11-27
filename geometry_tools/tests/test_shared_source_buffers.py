import unittest
import numpy as np
from geometry_tools.algorithms.shared_source_buffers import merge_source_updates


class SharedSourceBuffersTests(unittest.TestCase):
    def test_two_regions_share_the_same_edited_seam_and_retain_other_controls(self):
        source = {'positions': np.arange(18, dtype='<f4').reshape(6, 3),
                  'normals': np.zeros((6, 3), dtype='<i2')}
        a = source['positions'][[0, 1, 2]].copy()
        b = source['positions'][[2, 3, 4]].copy()
        a[2] += 1
        b[0] += 1
        b[1] += 2
        result = merge_source_updates(source, [([0, 1, 2], {'positions': a}),
                                                ([2, 3, 4], {'positions': b})])
        np.testing.assert_array_equal(result['positions'][2], a[2])
        np.testing.assert_array_equal(result['positions'][3], b[1])
        np.testing.assert_array_equal(result['positions'][5], source['positions'][5])
        np.testing.assert_array_equal(result['normals'], source['normals'])
        np.testing.assert_array_equal(source['positions'], np.arange(18, dtype='<f4').reshape(6, 3))

    def test_conflicting_seam_updates_fail_before_packing(self):
        source = {'positions': np.zeros((4, 3), dtype='<f4')}
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            merge_source_updates(source, [([1], {'positions': np.zeros((1, 3), '<f4')}),
                                           ([1], {'positions': np.ones((1, 3), '<f4')})])

    def test_source_attributes_cannot_change_dtype(self):
        with self.assertRaisesRegex(ValueError, 'dtype'):
            merge_source_updates({'normals': np.zeros((2, 3), '<i2')},
                                 [([0], {'normals': np.zeros((1, 3), '<f4')})])


if __name__ == '__main__':
    unittest.main()
