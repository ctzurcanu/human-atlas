import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.commands.clear_hand_components import signed_proposal
from geometry_tools.commands.render_hand_clearance import raster


class HandSkinProposalTests(unittest.TestCase):
    def test_open_skin_uses_local_outward_normal_without_volume_occupancy(self):
        mesh = Mesh(np.array([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]],float)*.01,
                    np.array([[0,1,2],[0,2,3]]), 'fixture-metres')
        index = SurfaceIndex(mesh)
        points = np.array([[.001,.002,.003],[.001,.002,-.004]])
        np.testing.assert_allclose(signed_proposal(index,points,open_oriented_surface=True),
                                   [.003,-.004],atol=1e-12)

    def test_skin_edge_side_is_only_a_distance_proposal(self):
        mesh = Mesh(np.array([[0,0,0],[.01,0,0],[0,.01,0]]),np.array([[0,1,2]]),'fixture-metres')
        index=SurfaceIndex(mesh)
        # Beyond the finite triangle, the local normal does not establish a volume.
        point=np.array([[.02,0,-.002]])
        value=signed_proposal(index,point,open_oriented_surface=True)
        self.assertAlmostEqual(value[0],-np.sqrt(.01**2+.002**2),places=12)

    def test_diagnostic_depth_keeps_the_front_triangle(self):
        front=np.array([[-.004,-.004,.003],[.004,-.004,.003],[0,.004,.003]])
        back=front-[0,0,.005]
        _,depth=raster(np.array([front,back]),[255,0,0],np.eye(3),np.zeros(3),1000,32)
        self.assertAlmostEqual(depth[16,16],.003,places=12)


if __name__=='__main__': unittest.main()
