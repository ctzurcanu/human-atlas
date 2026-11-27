import unittest
import numpy as np
from geometry_tools.algorithms.surface_ray_corridor import literal_surface_rays


class SurfaceRayCorridor(unittest.TestCase):
    def test_known_layer_depth_and_outside_layer(self):
        plane=np.array([[[-.01,-.01,-.002],[.01,-.01,-.002],[0.,.01,-.002]]])
        origins=np.array([[0.,0.,.003],[.02,0.,.003]])
        dirs=np.array([[0.,0.,-1.],[0.,0.,-1.]])
        result,proof=literal_surface_rays(plane,origins,dirs)
        self.assertAlmostEqual(result['distance_m'][0]-.003,.002,places=14)
        self.assertTrue(result['hit'][0]);self.assertFalse(result['hit'][1])
        self.assertTrue(np.isinf(result['distance_m'][1]))
        self.assertTrue(proof['unhitDirectionsRemainUnknown'])
        raised=plane.copy();raised[:,:,2]=.001
        result,_=literal_surface_rays(raised,origins[:1],dirs[:1])
        self.assertAlmostEqual(result['distance_m'][0]-.003,-.001,places=14)

    def test_parallel_ray_is_unknown_and_world_translation_preserves_depth(self):
        plane=np.array([[[-.01,-.01,-.002],[.01,-.01,-.002],[0.,.01,-.002]]])
        o=np.array([[0.,0.,.003]])
        result,_=literal_surface_rays(plane,o,np.array([[1.,0.,0.]]))
        self.assertFalse(result['hit'][0])
        offset=np.array([.1,1.7,-.04])
        result,_=literal_surface_rays(plane+offset,o+offset,np.array([[0.,0.,-1.]]))
        self.assertAlmostEqual(result['distance_m'][0],.005,places=14)


if __name__=='__main__':unittest.main()
