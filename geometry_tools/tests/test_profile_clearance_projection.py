import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.profile_clearance_projection import project_profile_centres


class ProfileClearanceTests(unittest.TestCase):
    def test_closed_box_projection_keeps_protected_control(self):
        v = np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],
                      [-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]], float)*.005
        f = np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],
                      [0,1,5],[0,5,4],[1,2,6],[1,6,5],
                      [2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        index = SurfaceIndex(Mesh(v, f, 'fixture-metres'))
        p = np.array([[.001,0,0],[.006,0,0],[0,0,0]])
        q, evidence = project_profile_centres(index, p, np.full(3,.001), [True, True, False])
        np.testing.assert_array_equal(q[2], p[2])
        self.assertGreaterEqual(float(q[:2,0].min()), .00624999)
        self.assertTrue(evidence['protectedCentresLiteral'])
        self.assertFalse(evidence['wholeSurfaceClearanceCertified'])
        self.assertGreater(evidence['steps'][0]['occupiedProfiles'], 0)


if __name__ == '__main__':
    unittest.main()
