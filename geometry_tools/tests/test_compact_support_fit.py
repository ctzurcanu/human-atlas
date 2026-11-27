import unittest
import numpy as np
from geometry_tools.algorithms.compact_support_fit import fit_support

class CompleteSupportFit(unittest.TestCase):
 def test_encloses_complete_cohort_and_respects_fixed_lower_extent(self):
  rng=np.random.default_rng(347);q=rng.normal(size=(300,3));q/=np.maximum(np.linalg.norm(q,axis=1)[:,None],1);p=q*np.array([.02,.025,.03])+[0,1.7,.005];c,a,r=fit_support(p,[0,1.7,.005],[.05,.06,.07],minimum_y_m=1.65)
  self.assertLessEqual(np.linalg.norm((p-c)/a,axis=1).max(),.8+1e-8);self.assertGreaterEqual(c[1]-a[1],1.65-1e-9);self.assertLess(r['volumeRatio'],.5);self.assertEqual(c[0],0)

if __name__=='__main__':unittest.main()
