import unittest
import numpy as np
from geometry_tools.algorithms.tetra_edge_refinement import refine_tetrahedral_edges

class TetraEdgeRefinementTests(unittest.TestCase):
 def test_neighboring_tetrahedra_preserve_volume_and_shared_faces(self):
  p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]]);t=np.array([[0,1,2,3],[0,2,1,4]]);r=refine_tetrahedral_edges(p,t);v=r['positions'];q=v[r['tetrahedra']];det=np.linalg.det(np.transpose(q[:,1:]-q[:,:1],(0,2,1)));self.assertTrue(np.all(det>0));self.assertAlmostEqual(det.sum()/6,1/3,places=13);np.testing.assert_array_equal(v[:len(p)],p);faces=np.concatenate([r['tetrahedra'][:,i]for i in [[1,2,3],[0,3,2],[0,1,3],[0,2,1]]]);keys,count=np.unique(np.sort(faces,axis=1),axis=0,return_counts=True);self.assertTrue(np.all(count<=2));self.assertEqual(np.sum(np.all(v[keys][:,:,2]==0,axis=1)&(count==2)),4)
