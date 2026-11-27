import unittest
import numpy as np
from geometry_tools.algorithms.section_axis_phase import phase_rigid_groups


class SectionAxisPhaseTests(unittest.TestCase):
    def test_source_groups_are_congruent_and_attachment_is_literal(self):
        section=np.array([[.003,0.,0.],[.003,.0004,0.],[.003,0.,.0004]])
        v=np.r_[section+[0,0,.005],section+[0,0,.025]]
        axis=np.c_[np.zeros(9),np.zeros(9),np.linspace(0,.04,9)]
        tangent=np.tile([0.,0.,1.],(9,1));labels=np.array([0,0,0,1,1,1])
        q,R,proof=phase_rigid_groups(v,labels,axis,tangent,20,[0,1,2])
        np.testing.assert_array_equal(q[:3],v[:3])
        self.assertGreater(np.linalg.norm(q[3:]-v[3:]),.0001)
        for ids in [np.arange(3),np.arange(3,6)]:
            np.testing.assert_allclose(np.linalg.norm(q[ids,None]-q[None,ids],axis=2),np.linalg.norm(v[ids,None]-v[None,ids],axis=2),atol=1e-17)
        np.testing.assert_allclose(np.linalg.det(R),1.,atol=1e-12)
        self.assertFalse(proof['axisMoved'])

    def test_out_of_bound_angle_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'30 degrees'):
            phase_rigid_groups(np.zeros((3,3)),[0,0,0],np.array([[0.,0.,0.],[0.,0.,.01]]),np.tile([0.,0.,1.],(2,1)),40,[0])


if __name__=='__main__':unittest.main()
