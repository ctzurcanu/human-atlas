import unittest
import numpy as np
from geometry_tools.algorithms.local_translation_cage import LocalTranslationCage

class TranslationCageTests(unittest.TestCase):
    def test_fixed_exterior_and_positive_analytic_derivative(self):
        cage=LocalTranslationCage(np.array([[0,.015,.065,.095],[1.3,1.32,1.42,1.44],[-.04,-.02,.04,.06]]),np.array([.008,0,0]))
        rng=np.random.default_rng(508);q=rng.uniform([-.02,1.28,-.06],[.12,1.46,.08],(500,3));j=cage.jacobians(q);step=1e-7
        measured=np.stack([(cage.map(q+np.eye(3)[i]*step)-cage.map(q-np.eye(3)[i]*step))/(2*step)for i in range(3)],axis=2)
        self.assertLess(float(abs(j-measured).max()),2e-8);self.assertGreaterEqual(float(np.linalg.det(j).min()),cage.determinant_lower_bound-1e-12)
        outside=np.array([[-.01,1.35,0],[.12,1.35,0],[.04,1.48,0],[.04,1.35,.09]])
        np.testing.assert_array_equal(cage.map(outside),outside)
    def test_excessive_translation_is_rejected(self):
        with self.assertRaises(ValueError):LocalTranslationCage(np.tile([0.,1.,2.,3.],(3,1)),np.array([1.,1.,1.]))

if __name__=='__main__':unittest.main()
