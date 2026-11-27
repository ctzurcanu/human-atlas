import unittest
import numpy as np
from geometry_tools.algorithms.paired_sheet_pinch import restore_paired_pinch


class PairedSheetPinchTests(unittest.TestCase):
    def fixture(self):
        v = np.array([[0,0,0], [1,0,0], [0,1,0], [-1,-1,0],
                      [0,0,0], [1,0,-.01], [0,1,-.01], [-1,-1,-.01]], float)
        f = np.array([[0,1,2],[0,2,3],[0,3,1],[4,6,5],[4,7,6],[4,5,7]])
        return v, f, [[1,5],[2,6],[3,7]]

    def test_complete_opposed_fans_restore_only_collapsed_inner_point(self):
        v,f,p = self.fixture()
        q, proof = restore_paired_pinch(v,f,0,4,p,maximum_motion_m=.02)
        np.testing.assert_allclose(q[4], [0,0,-.01])
        np.testing.assert_array_equal(q[np.arange(len(v)) != 4],v[np.arange(len(v)) != 4])
        self.assertFalse(proof['anatomicalThicknessEstablished'])

    def test_incomplete_fan_rejected(self):
        v,f,p = self.fixture()
        with self.assertRaises(AssertionError):
            restore_paired_pinch(v,f[:-1],0,4,p,maximum_motion_m=.02)

    def test_wrong_correspondence_rejected(self):
        v,f,p = self.fixture()
        p[0][1],p[1][1] = p[1][1],p[0][1]
        with self.assertRaises(AssertionError):
            restore_paired_pinch(v,f,0,4,p,maximum_motion_m=.02)

    def test_outward_or_over_limit_rejected(self):
        v,f,p = self.fixture()
        with self.assertRaises(AssertionError):
            restore_paired_pinch(v,f,0,4,p,maximum_motion_m=.005)
        v[5:,2] *= -1
        with self.assertRaises(AssertionError):
            restore_paired_pinch(v,f,0,4,p,maximum_motion_m=.02)


if __name__ == '__main__':
    unittest.main()
