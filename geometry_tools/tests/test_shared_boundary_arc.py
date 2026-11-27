import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.shared_boundary_arc import reconcile_shared_arc


class SharedArcTests(unittest.TestCase):
    def fixture(self,offset=2e-7):
        theta=np.arange(6)*np.pi/3;a=np.column_stack([.001*np.cos(theta),.001*np.sin(theta),np.zeros(6)])
        b=a[[3,2,1,0,5,4]].copy();b[[1,2],2]+=offset;b[[4,5],2]-=.0007
        def annulus(rim,height):
            v=np.r_[rim,rim+[0,0,height]];f=np.array([t for i in range(6)for t in [[i,(i+1)%6,6+i],[(i+1)%6,6+(i+1)%6,6+i]]]);return v,f
        av,af=annulus(a,.002);bv,bf=annulus(b,-.002);m=Mesh(np.r_[av,bv],np.r_[af,bf+12],'fixture')
        return m,list(range(12)),list(range(12,24)),list(range(6)),list(range(12,18))

    def test_two_split_seam_points_join_the_annuli_without_caps(self):
        args=self.fixture();old=args[0];new,proof=reconcile_shared_arc(*args,maximum_motion_m=5e-7)
        np.testing.assert_array_equal(new.faces,old.faces);self.assertEqual(proof['newFaces'],0);self.assertEqual(proof['newCaps'],0)
        self.assertEqual(proof['topology']['boundaryLoopLengths'],[6,6,6]);self.assertEqual(proof['topology']['eulerCharacteristic'],-1)
        self.assertEqual(proof['topology']['edgeConnectedComponents'],1)
        moved=np.any(new.vertices!=old.vertices,axis=1);self.assertEqual(int(moved.sum()),2)
        np.testing.assert_array_equal(new.vertices[~moved],old.vertices[~moved])

    def test_distinct_source_curves_are_not_welded(self):
        with self.assertRaises(ValueError):reconcile_shared_arc(*self.fixture(offset=2e-6),maximum_motion_m=5e-7)

    def test_partial_source_partition_is_rejected(self):
        args=list(self.fixture());args[2]=args[2][1:]
        with self.assertRaises(ValueError):reconcile_shared_arc(*args,maximum_motion_m=5e-7)

    def test_wrong_winding_cannot_form_an_internal_shared_edge(self):
        args=list(self.fixture());m=args[0];f=m.faces.copy();f[12:]=f[12:, [0,2,1]];args[0]=Mesh(m.vertices,f,m.frame)
        with self.assertRaises(ValueError):reconcile_shared_arc(*args,maximum_motion_m=5e-7)


if __name__=='__main__':unittest.main()
