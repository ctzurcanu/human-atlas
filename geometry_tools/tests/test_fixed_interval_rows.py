import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.fixed_interval_rows import preservation_rows


def sheets():
    triangle = np.array([[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.]])
    return Mesh(np.r_[triangle, triangle + [0., 0., .0003]],
                np.array([[0, 1, 2], [3, 4, 5]]), 'test:metres')


class FixedIntervalRowsTests(unittest.TestCase):
    def test_moving_higher_face_keeps_separation_from_fixed_lower_face(self):
        mesh = sheets()
        station = np.array([-1, -1, -1, 0, 0, 0])
        A, b, proof = preservation_rows(mesh, mesh, station, station, 1,
                                        maximum_motion_m=.001, self_pair=True)
        self.assertEqual(proof['preservedIntervals'], 1)
        self.assertTrue(np.all(A @ np.zeros(3) >= b))
        self.assertTrue(np.any(A @ np.array([0., 0., -.0006]) < b))
        self.assertTrue(np.all(A @ np.array([0., 0., .0006]) >= b))

    def test_two_moving_sections_accept_common_translation_reject_closure(self):
        mesh = sheets()
        station = np.array([0, 0, 0, 1, 1, 1])
        A, b, proof = preservation_rows(mesh, mesh, station, station, 2,
                                        maximum_motion_m=.001, self_pair=True)
        self.assertEqual(proof['preservedIntervals'], 1)
        common = np.array([[.0001, 0., .0002], [.0001, 0., .0002]]).ravel()
        self.assertTrue(np.all(A @ common >= b))
        closing = np.array([[0., 0., .0002], [0., 0., -.0002]]).ravel()
        self.assertTrue(np.any(A @ closing < b))

    def test_all_fixed_source_has_no_motion_rows(self):
        mesh = sheets()
        fixed = np.full(6, -1)
        A, b, proof = preservation_rows(mesh, mesh, fixed, fixed, 1,
                                        maximum_motion_m=.001, self_pair=True)
        self.assertEqual(A.shape, (0, 3))
        self.assertEqual(len(b), 0)
        self.assertEqual(proof['sourceMovingFaceIds'], [])


if __name__ == '__main__':
    unittest.main()
