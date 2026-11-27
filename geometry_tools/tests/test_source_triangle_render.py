import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.source_triangle_render import SourceTriangleRender

class NativeRenderPrecision(unittest.TestCase):
 def test_render_preserves_sliver_indices_without_relaxing_queries(self):
  vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[.9,.9,.2],[.9+1e-10,.9,.2],[.9,.9+1e-10,.2]])
  mesh=Mesh(vertices,np.array([[0,1,2],[3,4,5]]),'test:metres')
  with self.assertRaisesRegex(ValueError,'Float32'):SurfaceIndex(mesh)
  renderer=SourceTriangleRender(mesh);self.assertEqual(renderer.render_evidence['sourceTriangleCount'],2);self.assertEqual(renderer.render_evidence['float32NewlyCollapsedFaces'],1);self.assertFalse(renderer.render_evidence['clearancePredicate']);np.testing.assert_array_equal(mesh.vertices,vertices)

if __name__=='__main__':unittest.main()
