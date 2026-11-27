import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.directional_clearance import directional_exit_planes


def box(origin=(0., 0., 0.), scale=1.):
    vertices = np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],
                         [0,0,1],[1,0,1],[1,1,1],[0,1,1]], float)
    faces = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],
                      [0,1,5],[0,5,4],[3,7,6],[3,6,2],
                      [0,4,7],[0,7,3],[1,2,6],[1,6,5]])
    return SurfaceIndex(Mesh(vertices*scale+origin, faces, 'test:metres'))


class DirectionalClearanceTest(unittest.TestCase):
    def test_inside_and_behind_use_same_exit_side(self):
        result = directional_exit_planes(box(), [[.25,.4,.1],[.25,.4,.9],[.25,.4,-2]], [0,0,1])
        self.assertTrue(result['valid'].all())
        np.testing.assert_allclose(result['points'][:,2], 1, atol=1e-7)
        np.testing.assert_allclose(result['normals'], [[0,0,1]]*3)
        np.testing.assert_allclose(result['gap_m'], [-.9,-.1,-3], atol=1e-7)

    def test_opposite_exit_and_projection_miss(self):
        result = directional_exit_planes(box(), [[.3,.4,.6],[2,.4,.6]], [0,0,-2])
        np.testing.assert_array_equal(result['valid'], [True,False])
        self.assertAlmostEqual(result['points'][0,2], 0, places=6)
        np.testing.assert_allclose(result['normals'][0], [0,0,-1])
        self.assertTrue(np.isnan(result['gap_m'][1]))

    def test_world_frame_and_small_units(self):
        origin = np.array([.03,1.14,-.06]); scale = .004
        points = origin + scale*np.array([[.3,.25,.6],[.3,1.5,.6]])
        result = directional_exit_planes(box(origin,scale), points, [0,1,0])
        self.assertTrue(result['valid'].all())
        np.testing.assert_allclose(result['points'][:,1], origin[1]+scale, atol=1e-8)
        np.testing.assert_allclose(result['gap_m'], [-.003,.002], atol=1e-8)

    def test_explicit_direction_and_empty_query(self):
        index=box()
        with self.assertRaises(ValueError):
            directional_exit_planes(index, [[0,0,0]], [0,0,0])
        result=directional_exit_planes(index, np.empty((0,3)), [1,0,0])
        self.assertEqual(len(result['valid']),0)
