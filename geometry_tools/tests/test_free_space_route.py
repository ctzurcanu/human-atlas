import unittest
import numpy as np
from geometry_tools.algorithms.free_space_route import free_space_route,sampled_clearance
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.tests.test_sparse_group_clearance import cube


class FreeRouteTests(unittest.TestCase):
    def test_conduit_routes_around_fixed_cube_without_cutting_it(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        start=np.array([-.004,.005,.005]);end=np.array([.014,.005,.005])
        path,proof=free_space_route(start,end,bone,skin,np.array([start,end]),
            radius_m=.0002,clearance_m=.0001,grid_m=.001,padding_m=.008,maximum_grid_points=100000)
        np.testing.assert_array_equal(path[[0,-1]],[start,end])
        self.assertGreater(len(path),2)
        samples=np.concatenate([np.linspace(a,b,100)for a,b in zip(path[:-1],path[1:])])
        self.assertGreater(sampled_clearance(samples,bone,skin).min(),.00029)
        self.assertTrue(proof['fullTriangleAuditRequired'])


if __name__=='__main__':unittest.main()
