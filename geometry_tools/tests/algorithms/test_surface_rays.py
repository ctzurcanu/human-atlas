import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.surface_rays import first_surface_hits


class OriginalSurfaceRayTests(unittest.TestCase):
    def planes(self):
        vertices=np.array([[0,0,1],[1,0,1],[1,1,1],[0,1,1],
                           [0,0,2],[1,0,2],[1,1,2],[0,1,2]],float)
        return Mesh(vertices,np.array([[0,1,2],[0,2,3],[4,5,6],[4,6,7]]),'test:metres')

    def test_first_positive_surface_and_actual_barycentric_witness(self):
        mesh=self.planes()
        result,evidence=first_surface_hits(mesh,[[.2,.3,0],[.2,.3,3]],[[0,0,5],[0,0,-7]])
        np.testing.assert_allclose(result['distance_m'],[1,1],atol=1e-15)
        np.testing.assert_allclose(np.einsum('ij,ijk->ik',result['barycentric'],mesh.triangles[result['face_ids']]),result['points'],atol=1e-15)
        self.assertTrue(evidence['original_coordinates_preserved'])
        self.assertTrue(evidence['complete_within_budget'])

    def test_self_hit_is_excluded_and_range_is_respected(self):
        mesh=self.planes()
        result,_=first_surface_hits(mesh,[[.2,.3,1],[.2,.3,0]],[[0,0,1],[0,0,1]],minimum_m=1e-8,maximum_m=.75)
        self.assertFalse(result['hit'].any())

    def test_parallel_and_outside_rays_miss(self):
        result,_=first_surface_hits(self.planes(),[[.2,.3,0],[2,2,0]],[[1,0,0],[0,0,1]])
        self.assertFalse(result['hit'].any())

    def test_original_sub_float32_separation_survives(self):
        z=1.25
        vertices=np.array([[0,0,z],[1,0,z],[0,1,z],
                           [0,0,z+2e-10],[1,0,z+2e-10],[0,1,z+2e-10]],float)
        mesh=Mesh(vertices,np.array([[0,1,2],[3,4,5]]),'test:metres')
        result,_=first_surface_hits(mesh,[[.2,.2,z+1e-10]],[[0,0,1]],maximum_m=1e-8)
        self.assertEqual(result['face_ids'][0],1)
        self.assertLess(abs(result['distance_m'][0]-1e-10),1e-16)

    def test_budget_exhaustion_cannot_return_an_incomplete_success(self):
        with self.assertRaisesRegex(RuntimeError,'triangle budget'):
            first_surface_hits(self.planes(),[[.2,.3,0]],[[0,0,1]],max_triangle_tests=1)

    def test_invalid_direction_is_rejected(self):
        with self.assertRaises(ValueError):
            first_surface_hits(self.planes(),[[.2,.3,0]],[[0,0,0]])


if __name__=='__main__':unittest.main()
