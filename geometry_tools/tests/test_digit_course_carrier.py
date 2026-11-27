import unittest
import numpy as np
from geometry_tools.algorithms.digit_course_carrier import DigitCourseCarrier
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores


class DigitCourseTests(unittest.TestCase):
    def test_identity_shared_courses(self):
        a = np.array([[0., 0., 0.], [0., .03, .002], [.005, .07, .004]])
        carrier = DigitCourseCarrier([dict(source=a, target=a),
                                     dict(source=a+[.02, 0., 0.], target=a+[.02, 0., 0.])])
        points = np.array([[.001, .02, .005], [.015, .05, -.002]])
        q, R, _ = carrier.poses(points)
        np.testing.assert_allclose(q, points, atol=1e-15)
        np.testing.assert_allclose(R, np.repeat(np.eye(3)[None], len(points), axis=0), atol=1e-14)

    def test_source_root_literal_distal_translation(self):
        a = np.array([[0., 0., 0.], [0., .05, 0.], [0., .1, 0.]])
        carrier = DigitCourseCarrier([dict(source=a, target=a+[.02, 0., 0.])])
        p = np.array([[.003, 0., 0.], [.003, .08, 0.]])
        q, _, _ = carrier.poses(p)
        np.testing.assert_array_equal(q[0], p[0])
        np.testing.assert_allclose(q[1], p[1]+[.02, 0., 0.], atol=1e-15)

    def test_curved_course_rigid_sections_keep_calibre(self):
        a = np.array([[0., 0., 0.], [0., .05, 0.], [0., .1, 0.]])
        b = np.array([[0., 0., 0.], [.02, .05, .005], [.04, .09, .01]])
        centres = np.array([[0., .045, 0.], [0., .07, 0.]])
        q, R, _ = DigitCourseCarrier([dict(source=a, target=b)]).poses(centres)
        theta = np.arange(6)*np.pi/3
        ring = np.c_[.001*np.cos(theta), np.zeros(6), .0007*np.sin(theta)]
        v = np.concatenate([c+ring for c in centres])
        faces = []
        for i in range(6):
            j = (i+1)%6
            faces.extend([[i, j, j+6], [i, j+6, i+6]])
        result = transport_profile_cores(v, np.array(faces), [np.arange(6), np.arange(6, 12)], q, R)
        for ids in [np.arange(6), np.arange(6, 12)]:
            old = np.linalg.norm(v[ids, None]-v[ids][None], axis=2)
            new = np.linalg.norm(result['vertices'][ids, None]-result['vertices'][ids][None], axis=2)
            np.testing.assert_allclose(old, new, atol=2e-17)
        np.testing.assert_allclose(np.linalg.det(R), 1., atol=1e-14)


if __name__ == '__main__':
    unittest.main()
