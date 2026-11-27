import unittest
import numpy as np
from geometry_tools.algorithms.tube_axis_map import axis_transport,endpoint_course

class AxisTransportTests(unittest.TestCase):
    def test_constant_proper_rotation_preserves_all_pairwise_distances(self):
        old=np.array([[0.,0.,0.],[0.,0.,1.],[0.,0.,2.]])
        r=np.array([[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]])
        new=old@r.T+[2.,3.,4.];points=np.random.default_rng(37).normal(size=(25,3))
        actual,_=axis_transport(points,old,new,np.tile(r,(3,1,1)))
        np.testing.assert_allclose(actual,points@r.T+[2.,3.,4.],rtol=0,atol=2e-15)
    def test_fixed_tail_replays_exactly_and_nearby_branch_gets_same_offset(self):
        old=np.column_stack([np.zeros(5),np.zeros(5),np.arange(5.)]);new=old.copy();new[:2,0]=[.2,.1]
        r=np.tile(np.eye(3),(5,1,1));points=np.array([[.01,.04,.5],[.017,-.027,3.14159]])
        actual,_=axis_transport(points,old,new,r)
        np.testing.assert_allclose(actual[0],points[0]+[.15,0,0],rtol=0,atol=1e-16)
        np.testing.assert_array_equal(actual[1],points[1])
    def test_reflection_is_rejected(self):
        old=np.array([[0.,0.,0.],[0.,0.,1.]])
        with self.assertRaises(ValueError):axis_transport(old,old,old,np.tile(np.diag([-1.,1.,1.]),(2,1,1)))
    def test_unit_determinant_shear_cannot_change_wall_sections(self):
        old=np.array([[0.,0.,0.],[0.,0.,1.]])
        shear=np.array([[1.,.3,0.],[0.,1.,0.],[0.,0.,1.]])
        with self.assertRaises(ValueError):axis_transport(old,old,old,np.tile(shear,(2,1,1)))

if __name__=='__main__':unittest.main()
