import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.bone_radial_field import last_surface_hit_before


class FiniteSurfaceRayTests(unittest.TestCase):
    def setUp(self):
        v=np.array([[x,y,z]for z in [.003,.007]for x,y in [(-.01,-.01),(.01,-.01),(0,.01)]])
        self.index=SurfaceIndex(Mesh(v,np.array([[0,1,2],[3,4,5]]),'test:m'))

    def test_hits_past_skin_limit_do_not_close_the_corridor(self):
        origins=np.zeros((3,3));directions=np.tile([0.,0.,1.],(3,1))
        result=last_surface_hit_before(self.index,origins,directions,[.005,.009,.002])
        np.testing.assert_allclose(result['distance_m'],[.003,.007,0],atol=1e-12)
        np.testing.assert_array_equal(result['face_ids'],[0,1,-1])
        self.assertFalse(result['solidClearanceCertified'])

    def test_invalid_ray_intervals_and_scaled_directions_are_rejected(self):
        for direction,limit in [([0,0,2],.01),([0,0,1],0),([0,0,1],float('nan'))]:
            with self.assertRaises(ValueError):last_surface_hit_before(self.index,[[0,0,0]],[direction],[limit])
