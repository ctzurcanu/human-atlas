import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.parallel_surface_rays import cast_axis_rays


class ParallelSurfaceRayTests(unittest.TestCase):
    def test_both_directions_choose_nearest_layer_and_report_misses(self):
        v=np.array([[0,0,0],[1,0,0],[0,0,1],[0,1,0],[1,1,0],[0,1,1]],float)
        m=Mesh(v,np.array([[0,1,2],[3,4,5]]),'test')
        upper,proof=cast_axis_rays(m,[[.2,2,.2],[2,2,2]],sign=-1)
        self.assertEqual(upper['face_ids'].tolist(),[1,-1])
        np.testing.assert_allclose(upper['points'][0],[.2,1,.2],atol=1e-15)
        self.assertTrue(np.isnan(upper['points'][1]).all())
        lower,_=cast_axis_rays(m,[[.2,-1,.2]],sign=1)
        self.assertEqual(lower['face_ids'][0],0)
        self.assertFalse(upper['boundary_ambiguous'][0])
        self.assertTrue(proof['source_positions_preserved'])

    def test_literal_thin_triangle_is_not_float32_collapsed(self):
        v=np.array([[1,1,1],[1+1e-10,1,1],[1,1,1+1e-10]],float)
        m=Mesh(v,np.array([[0,1,2]]),'test')
        rays,_=cast_axis_rays(m,[[1+2e-11,2,1+2e-11]])
        self.assertTrue(rays['hit'][0])
        self.assertEqual(rays['points'][0,1],1)
        np.testing.assert_array_equal(m.vertices,v)

    def test_boundary_and_finite_work_budget_are_explicit(self):
        m=Mesh([[0,0,0],[1,0,0],[0,0,1],[0,1,0],[1,1,0],[0,1,1]],np.array([[0,1,2],[3,4,5]]),'test')
        rays,_=cast_axis_rays(m,[[0,2,0]])
        self.assertTrue(rays['boundary_ambiguous'][0])
        with self.assertRaisesRegex(RuntimeError,'budget'):
            cast_axis_rays(m,[[.2,2,.2]],max_triangle_tests=1)


if __name__=='__main__':unittest.main()
