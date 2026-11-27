import unittest
import numpy as np
from geometry_tools.algorithms.float32_diagonal_repair import repair_diagonals
from geometry_tools.algorithms.shared_triangle_edges import boundary

class Float32DiagonalTests(unittest.TestCase):
    def test_explicit_multichart_encoding_controls_the_repair(self):
        v=np.array([[.001,1.001,0.],[.02,1.02,0.],[.0105,1.0105000001,0.],[.02,1.001,0.]])
        f=np.array([[0,1,2],[1,0,3]])
        q=(v-[0,1.,0]).astype('<f4').astype(float)+[0,1.,0]
        out,r=repair_diagonals(v,f,[0,0,0],stored_positions=q)
        self.assertGreater(r['diagonalsChanged'],0)
        self.assertEqual(boundary(out),boundary(f))
        t=q[out];self.assertTrue(np.all(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2]>0))

    def test_planar_sliver_repaired_without_moving_boundary(self):
        v=np.array([[.001,1.001,0.],[.02,1.02,0.],[.0105,1.0105000001,0.],[.02,1.001,0.]])
        f=np.array([[0,1,2],[1,0,3]]);q=(v-[0,1.,0]).astype('<f4').astype(float)+[0,1.,0]
        out,r=repair_diagonals(v,f,[0,1.,0]);self.assertGreater(r['diagonalsChanged'],0);self.assertEqual(boundary(out),boundary(f));self.assertTrue(r['originalControlsLiteral'])
        t=q[out];self.assertTrue(np.all(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])[:,2]>0))
    def test_real_nonplanar_fold_not_flattened(self):
        v=np.array([[.001,1.001,0.],[.02,1.02,0.],[.0105,1.0105000001,0.],[.02,1.001,.01]])
        with self.assertRaises(ValueError):repair_diagonals(v,np.array([[0,1,2],[1,0,3]]),[0,1.,0])

if __name__=='__main__':unittest.main()
