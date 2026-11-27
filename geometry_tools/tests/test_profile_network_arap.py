import unittest
import numpy as np
from geometry_tools.algorithms.profile_network_arap import fit_profile_network


class ProfileNetworkArapTest(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0,0],[1.,0,0],[1.,1,0],[0.,1,0]])
        self.e=[[0,1],[1,2],[2,3],[3,0],[0,2]]

    def test_unchanged_graph_is_retained(self):
        r=fit_profile_network(self.p,self.p,self.e,fixed=[0])
        np.testing.assert_allclose(r['centres'],self.p,atol=1e-14)
        self.assertTrue(r['fixed_nodes_literal'])

    def test_whole_proper_pose_is_exact(self):
        R=np.array([[0.,-1,0],[1.,0,0],[0,0,1]])
        target=self.p@R.T+[.3,.5,1]
        r=fit_profile_network(self.p,target,self.e,reference_rotations=np.repeat(R[None],4,axis=0))
        np.testing.assert_allclose(r['centres'],target,atol=1e-14)
        np.testing.assert_allclose(r['rotations'],np.repeat(R[None],4,axis=0),atol=1e-14)

    def test_fixed_recipient_nodes_do_not_move_with_target(self):
        target=self.p+[.2,.3,.4]
        r=fit_profile_network(self.p,target,self.e,fixed=[0,1])
        np.testing.assert_array_equal(r['centres'][:2],self.p[:2])
        self.assertLess(np.max(abs(np.linalg.det(r['rotations'])-1)),1e-12)
        self.assertFalse(r['adjacency_inferred'])

    def test_reflected_reference_is_rejected(self):
        with self.assertRaises(ValueError):
            fit_profile_network(self.p,self.p,self.e,reference_rotations=np.repeat(np.diag([-1.,1,1])[None],4,axis=0))

    def test_explicit_mixed_source_and_moving_attachment_anchors_are_literal(self):
        R=np.array([[0.,-1,0],[1.,0,0],[0,0,1]])
        target=self.p.copy();target[2:]=self.p[2:]@R.T+[.2,.3,.4]
        rotations=np.repeat(np.eye(3)[None],4,axis=0);rotations[2:]=R
        r=fit_profile_network(self.p,target,self.e,fixed=[0,2,3],
            reference_rotations=rotations,fixed_at_targets=True)
        np.testing.assert_array_equal(r['centres'][[0,2,3]],target[[0,2,3]])
        np.testing.assert_array_equal(r['rotations'][[0,2,3]],rotations[[0,2,3]])
        self.assertTrue(r['fixed_nodes_literal']);self.assertTrue(r['fixed_rotations_literal'])
        self.assertEqual(r['fixed_node_reference'],'explicit-target')
        self.assertTrue(np.isfinite(r['energies']).all())

    def test_wholly_anchored_target_pose_is_retained_without_a_free_solve(self):
        R=np.array([[0.,-1,0],[1.,0,0],[0,0,1]])
        target=self.p@R.T+[.3,.5,1];rotations=np.repeat(R[None],4,axis=0)
        r=fit_profile_network(self.p,target,self.e,fixed=range(4),
            reference_rotations=rotations,fixed_at_targets=True)
        np.testing.assert_array_equal(r['centres'],target)
        np.testing.assert_array_equal(r['rotations'],rotations)


if __name__=='__main__':unittest.main()
