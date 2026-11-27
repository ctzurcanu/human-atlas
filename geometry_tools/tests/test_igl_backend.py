import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.igl_backend import closest_points,SurfaceIndex
from geometry_tools.spatial import BVH


class NativeDoubleDistanceTests(unittest.TestCase):
    def test_thin_triangle_and_bvh_agree(self):
        # This 10 pm edge disappears after conversion near y=1 to Float32.
        mesh = Mesh(np.array([[0.,1.,0.],[1.,1.,0.],[0.,1.+1e-11,0.],
                              [0.,0.,2.],[1.,0.,2.],[0.,1.,2.]]),
                    np.array([[0,1,2],[3,4,5]]), 'test:metres')
        points=np.array([[.25,1.,.1],[.2,.2,1.9],[-.1,1.,0.],[.3,.4,-.2]])
        native, meta=closest_points(mesh,points)
        reference,_=BVH(mesh).closest_points(points)
        np.testing.assert_allclose(native['distance_m'],reference['distance_m'],atol=1e-13)
        np.testing.assert_allclose(native['points'],reference['points'],atol=1e-13)
        self.assertEqual(meta['precision'],'Float64')
        self.assertTrue(meta['original_source_triangles'])
        self.assertAlmostEqual(native['distance_m'][0],.1)
        index=SurfaceIndex(mesh)
        np.testing.assert_allclose(index.closest_points(points[::-1])[0]['distance_m'],native['distance_m'][::-1],atol=1e-13)

    def test_invalid_queries_rejected(self):
        mesh=Mesh(np.eye(3),np.array([[0,1,2]]),'test:metres')
        with self.assertRaises(ValueError):closest_points(mesh,[[np.nan,0,0]])
