import unittest
import numpy as np
from geometry_tools.algorithms.faithful_float32 import oriented_rounding
from geometry_tools.algorithms.tube_sections import TubeGeometryError


class OrientedFloat32Tests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0.,0.],[.001,.001,0.],[.002,.002+1e-10,0.],[.003,.004,0.]])
        self.f=np.array([[0,1,2],[1,2,3]])
        self.w=np.broadcast_to(np.array([0,0,32767],np.int16),(2,3,3)).copy()
        self.q=self.p.astype(np.float32).astype(float)

    def test_repairs_collapsed_facet_without_moving_protected_controls(self):
        result,evidence=oriented_rounding(self.p,self.q,self.f,self.w,mutable_controls=[2])
        np.testing.assert_array_equal(result[[0,1,3]],self.q[[0,1,3]])
        t=result[self.f];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        self.assertTrue(np.all(n[:,2]>0))
        self.assertEqual(evidence['diagnosed_faces'],1)
        self.assertEqual(evidence['changed_controls'],[2])
        np.testing.assert_array_equal(result.astype(np.float32).astype(float),result)

    def test_cannot_modify_an_immutable_rounding_failure(self):
        with self.assertRaises(TubeGeometryError):
            oriented_rounding(self.p,self.q,self.f,self.w,mutable_controls=[])

    def test_refuses_an_inward_authored_source_field(self):
        with self.assertRaises(TubeGeometryError):
            oriented_rounding(self.p,self.q,self.f,-self.w,mutable_controls=[2])


if __name__=='__main__':unittest.main()
