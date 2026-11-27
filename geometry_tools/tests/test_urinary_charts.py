"""Urinary chart seams must survive the renderer's actual Float32 encoding."""
import unittest

import numpy as np

from geometry_tools.urinary_delivery import shared_encoding


class UrinaryChartTests(unittest.TestCase):
    def test_explicit_attachment_chart_keeps_a_literal_native_bone_contact(self):
        origin = np.array([0., .9140625, 0.])
        points = np.array([[.0115, .902, .0286], [.012, .903, .029], [.01, .9025, .029]], dtype='<f4').astype(float)
        surfaces = {'bone-contact': np.array([[0, 1, 2]])}
        cells = {'named-muscle': [dict(surface='bone-contact', side=1)]}
        encoded, _ = shared_encoding(points, surfaces, cells, {'attachment': origin}, chart_for={'named-muscle': 'attachment'})
        np.testing.assert_array_equal(encoded, points)

    def test_shared_seam_is_exact_in_both_origins_without_moving_the_renal_root(self):
        origins = {'left': np.array([.0625, 1.1875, -.015625]),
                   'bladder': np.array([0., .916015625, 0.])}
        seam = np.array([[.0103456379, .9363124395, -.0233343528],
                         [.0115013579, .9392384913, -.0238035291],
                         [.0127506379, .9356953047, -.0239713971]])
        root = origins['left']+np.array([[-.03, -.05, .02], [-.029, -.05, .02], [-.03, -.049, .02]], dtype='<f4').astype(float)
        naive = [(seam-o).astype('<f4').astype(float)+o for o in origins.values()]
        self.assertFalse(np.array_equal(*naive))
        renal, bladder = 'RENALDETAIL:male:left:ureter-wall:mucosa', 'BLADDER:male:fundus:mucosa'
        surfaces = {'shared':np.array([[0,1,2]]), 'protected-root':np.array([[3,4,5]])}
        cells = {renal:[dict(surface='shared',side=1,neighbor=bladder), dict(surface='protected-root',side=1,neighbor=None)],
                 bladder:[dict(surface='shared',side=-1,neighbor=renal)]}
        points = np.r_[seam, root]
        encoded, evidence = shared_encoding(points, surfaces, cells, origins)
        for origin in origins.values():
            self.assertTrue(np.array_equal((encoded[:3]-origin).astype('<f4').astype(float)+origin, encoded[:3]))
        self.assertTrue(np.array_equal(encoded[3:], root))
        self.assertLess(np.linalg.norm(encoded-points,axis=1).max(), 3e-8)
        self.assertEqual(evidence['shared_controls'], 3)


if __name__ == '__main__':
    unittest.main()
