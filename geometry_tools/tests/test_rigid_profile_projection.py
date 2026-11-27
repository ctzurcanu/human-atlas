import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.rigid_profile_projection import project_rigid_profiles


class RigidProfileProjectionTests(unittest.TestCase):
    def test_aliases_and_overlapping_profile_groups_preserve_every_profile(self):
        p = np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,0]])
        m = Mesh(p, np.array([[0,1,2],[0,2,3]]), 'fixture')
        desired = p+np.arange(15).reshape(5,3)*.001
        desired[4] = desired[1]
        result, proof = project_rigid_profiles([m], [desired], [[np.array([0,1,2]),np.array([0,2,3])]], [np.zeros(5,bool)])
        np.testing.assert_allclose(result[0]-p, np.broadcast_to(result[0][0]-p[0],p.shape))
        self.assertEqual(proof['coupledRigidGroups'], 1)

    def test_one_fixed_control_fixes_the_whole_profile(self):
        p = np.array([[0.,0,0],[1,0,0],[0,1,0]])
        m = Mesh(p, np.array([[0,1,2]]), 'fixture')
        result, _ = project_rigid_profiles([m], [p+1], [[np.arange(3)]], [np.array([True,False,False])])
        np.testing.assert_array_equal(result[0],p)

    def test_two_meshes_share_one_exact_source_interface(self):
        p = np.array([[0.,0,0],[1,0,0],[0,1,0]])
        m = Mesh(p, np.array([[0,1,2]]), 'fixture')
        result, _ = project_rigid_profiles([m,m], [p+.002,p+.002], [[np.arange(3)],[np.arange(3)]], [np.zeros(3,bool),np.zeros(3,bool)])
        np.testing.assert_array_equal(result[0],result[1])

    def test_proper_rotation_follows_the_guide_without_rescaling_a_profile(self):
        p = np.array([[0.,0,0],[1,0,0],[0,1,0]])
        m = Mesh(p, np.array([[0,1,2]]), 'fixture')
        rotation = np.array([[0.,1,0],[-1,0,0],[0,0,1]])
        target = p@rotation+[.2,.3,.4]
        result, proof = project_rigid_profiles([m], [target], [[np.arange(3)]], [np.zeros(3,bool)], allow_rotation=True)
        np.testing.assert_allclose(result[0],target,atol=1e-14)
        self.assertTrue(proof['properRotationsAllowed'])


if __name__ == '__main__':
    unittest.main()
