import unittest
import numpy as np
from geometry_tools.algorithms.shared_planar_mouth import stitch_planar_mouth
from geometry_tools.algorithms.shared_triangle_edges import boundary

def box(z0,z1):
    v=np.array([[x,y,z]for z in [z0,z1]for x,y in [(0,0),(1,0),(1,1),(0,1)]],float)
    f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
    return v,f

class SharedPlanarMouthTests(unittest.TestCase):
    def test_exact_mouth_removes_only_internal_disks(self):
        a,af=box(0,1);b,bf=box(1,2);v,f,p=stitch_planar_mouth(a,af,b,bf,[0,0,1],1)
        self.assertFalse(boundary(f));self.assertEqual(len(f),20)
        t=v[f];volume=np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()/6
        self.assertAlmostEqual(volume,2);self.assertEqual(p['removedConstructionDiskTriangles'],[2,2]);self.assertTrue(np.array_equal(v[:len(a)],a))

    def test_geometrically_different_mouth_is_rejected(self):
        a,af=box(0,1);b,bf=box(1,2);b[:,0]+=.01
        with self.assertRaises(ValueError):stitch_planar_mouth(a,af,b,bf,[0,0,1],1)

if __name__=='__main__':unittest.main()
