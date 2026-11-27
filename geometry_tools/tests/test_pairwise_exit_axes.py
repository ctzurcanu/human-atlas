import unittest
import numpy as np
from geometry_tools.algorithms.pairwise_exit_axes import pairwise_exit_axes
from geometry_tools.distance import triangle_pairs
class PairwiseExitTests(unittest.TestCase):
 def test_different_route_directions_separate_every_actual_triangle_pair(self):
  a=np.array([[[-.01,-.01,0],[.01,-.01,0],[0,.01,0]]]*2);b=np.array([[[-.01,0,-.01],[.01,0,-.01],[0,0,.01]]]*2);rays=np.array([[0.,1.,0.],[0.,-1.,0.]])
  axes,amount=pairwise_exit_axes(a,b,rays,clearance_m=2e-5);q=a+amount[:,None,None]*rays[:,None,:]
  gaps=np.einsum('pvc,pc->pv',q,axes).min(1)-np.einsum('pvc,pc->pv',b,axes).max(1)
  measured=triangle_pairs(q,b)
  self.assertTrue(np.all(gaps>=2e-5-1e-14));self.assertTrue(np.all(measured['distance_m']>=2e-5-1e-14));self.assertFalse(np.any(np.isin(measured['relation'],['crossing','coplanar-overlap'])))
 def test_zero_route_is_rejected(self):
  a=np.array([[[0,0,0],[1,0,0],[0,1,0]]],float)
  with self.assertRaises(ValueError):pairwise_exit_axes(a,a,np.zeros((1,3)))
if __name__=='__main__':unittest.main()
