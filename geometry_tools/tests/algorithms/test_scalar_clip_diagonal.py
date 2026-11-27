import unittest
import numpy as np
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface

class ScalarClipDiagonalTests(unittest.TestCase):
 def test_bank_size_does_not_change_cut_diagonal_or_geometry(self):
  vertices=np.array([[0.,0,0],[2.,0,0],[0.,1,0]])
  faces=np.array([[0,1,2]]);values=np.array([-1.,1.,1.]);a=clip_scalar_surface(vertices,faces,values,0.)
  extra=np.repeat([[3.,4.,5.]],20000,axis=0);b=clip_scalar_surface(np.r_[vertices,extra],faces,np.r_[values,np.ones(len(extra))],0.)
  for key in ['positive_faces','negative_faces']:
   np.testing.assert_array_equal(a['positions'][a[key]],b['positions'][b[key]])
  total=0
  for key in ['positive_faces','negative_faces']:
   q=a['positions'][a[key]];total+=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1).sum()/2
  self.assertAlmostEqual(total,1.);self.assertEqual(a['caps_added'],0);self.assertTrue(a['original_controls_literal'])

if __name__=='__main__':unittest.main()
