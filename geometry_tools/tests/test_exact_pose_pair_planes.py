import unittest
import numpy as np
from geometry_tools.algorithms.exact_pose_pair_planes import ExactPosePairPlanes

class ExactPosePairPlaneTests(unittest.TestCase):
    def model(self):
        return ExactPosePairPlanes([[.001,.002,.003],[.003,.001,.002]],[[.002,.004,.001],[.001,0.,.003]],np.array([0,1]),np.array([1,-1]),[[1.,0.,0.],[0.,1.,0.]],[0.,0.],[[0.,0.,0.],[.002,.003,0.]],[.02,.03])
    def test_analytic_local_row_matches_independent_finite_difference(self):
        model=self.model();x=np.array([[.0002,0.,.0003],[0.,-.0001,.0002],[.0001,.0002,0.],[.0001,0.,.0003]])
        G,b,_=model.local_rows(x,np.arange(2));norm=np.sqrt(np.asarray(G.power(2).sum(1)).ravel())
        # Scaling is >=1; recover the derivative by evaluating the local row
        # difference at two points and compare with the actual replay slope.
        eps=1e-8;direction=np.array([[.1,-.2,.3],[-.1,.2,.1],[.2,.1,-.3],[.1,-.1,.2]])
        actual=(model.values(x+eps*direction)-model.values(x-eps*direction))/(2*eps)
        baseG=[]
        for j in range(x.size):
            dx=np.zeros_like(x);dx.ravel()[j]=eps
            baseG.append((model.values(x+dx)-model.values(x-dx))/(2*eps))
        derivative=np.array(baseG).T
        scales=np.maximum(1.,1/np.linalg.norm(derivative,axis=1))
        np.testing.assert_allclose(G@direction.ravel(),actual*scales,atol=2e-8,rtol=1e-7)
    def test_local_halfspace_does_not_certify_a_violating_current_pose(self):
        model=self.model();x=np.zeros((4,3));model.minimum=model.values(x)+.00001
        bad,proof=model.violations(x);self.assertEqual(len(bad),2)
        A,b,_=model.local_rows(x,bad)
        self.assertTrue(np.all(b-A@x.ravel()>0))
