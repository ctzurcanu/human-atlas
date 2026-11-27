import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.reference_interval_rows import reference_preservation_rows


class ReferenceIntervalTest(unittest.TestCase):
    def test_gap_is_based_on_current_pose_but_variables_are_absolute(self):
        v=np.array([[0.,0.,0.],[0.,.001,0.],[0.,0.,.001]])
        ref=Mesh(v,np.array([[0,1,2]]),'test')
        current=Mesh(v+[.0002,0,0],ref.faces,'test')
        obstacle=Mesh(v+[.00021,0,0],ref.faces,'test')
        A,b,p=reference_preservation_rows(current,obstacle,ref,obstacle,
            np.zeros(3,int),np.full(3,-1,int),1,maximum_motion_m=.0003)
        self.assertGreater(len(b),0)
        self.assertTrue(np.all(A@np.array([.0002,0,0])>=b-1e-12))
        self.assertFalse(np.all(A@np.array([.00022,0,0])>=b-1e-12))
        self.assertTrue(p['immutableReferenceCoordinates'])

    def test_only_both_repair_faces_are_released(self):
        v=np.array([[0.,0.,0.],[0.,.001,0.],[0.,0.,.001]])
        a=Mesh(v,np.array([[0,1,2]]),'test');b=Mesh(v+[.00001,0,0],a.faces,'test')
        args=(a,b,a,b,np.zeros(3,int),np.full(3,-1,int),1)
        A,_,_=reference_preservation_rows(*args,maximum_motion_m=.0003,
            released_first_faces=[0],released_second_faces=[])
        C,_,p=reference_preservation_rows(*args,maximum_motion_m=.0003,
            released_first_faces=[0],released_second_faces=[0])
        self.assertGreater(A.shape[0],0);self.assertEqual(C.shape[0],0)
        self.assertGreater(p['releasedRepairPatchPairs'],0)


if __name__=='__main__':unittest.main()
