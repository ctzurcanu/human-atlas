import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.spectral_tangent_cones import mixed_tangent_cones,mixed_tangent_residuals


class SpectralTangentCones(unittest.TestCase):
    def test_current_inherited_two_axis_strain_is_feasible_without_extra_budget(self):
        o=np.array([[.12,0.,0.,0.,-.04,0.]])
        A,b,cones,mask,proof=mixed_tangent_cones(sparse.csr_matrix((6,3)),o,bound=.12,reserve=0)
        self.assertTrue(mask[0]);self.assertGreater(np.linalg.norm(o),.12)
        p=0;M=np.zeros((5,5))
        for j in range(5):
            for i in range(j+1):M[i,j]=M[j,i]=b[p]/(np.sqrt(2) if i!=j else 1);p+=1
        np.testing.assert_allclose(M[:3,3:],o.reshape(2,3).T,atol=1e-16)
        self.assertGreaterEqual(float(np.linalg.eigvalsh(M).min()),-1e-15)
        self.assertLessEqual(float(mixed_tangent_residuals(o,mask)[0]),.12+1e-15)
        self.assertTrue(proof['originalReferencesAndBudgetsUnchanged'])

    def test_svec_linear_sign_and_exact_spectral_violation(self):
        C=sparse.csr_matrix(np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]]))
        o=np.array([[.12,0,0,0,-.04,0.]])
        A,b,_,mask,_=mixed_tangent_cones(C,o,bound=.12,reserve=0,extra_columns=1)
        x=np.array([.01,0,0,0]);slack=b-A@(x/.001);M=np.zeros((5,5));p=0
        for j in range(5):
            for i in range(j+1):M[i,j]=M[j,i]=slack[p]/(np.sqrt(2) if i!=j else 1);p+=1
        np.testing.assert_allclose(M[:3,3:],(o+(C@x[:3])[None]).reshape(2,3).T,atol=1e-16)
        self.assertLess(float(np.linalg.eigvalsh(M).min()),0)
        self.assertGreater(float(mixed_tangent_residuals(o+(C@x[:3])[None],mask)[0]),.12)


if __name__=='__main__':unittest.main()
