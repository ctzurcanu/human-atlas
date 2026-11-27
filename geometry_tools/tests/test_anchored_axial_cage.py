import unittest
import numpy as np
from geometry_tools.algorithms.anchored_axial_cage import AnchoredAxialCage


class AnchoredAxialCageTests(unittest.TestCase):
    def cage(self,axis=1):
        return AnchoredAxialCage([0,.3,.7,1],[[1,1],[.7,.8],[.8,.6],[1,1]],
            [[1,1],[.8,.7],[.6,.9],[1,1]],[0,.4,1],[[.01,-.03],[.02,-.02],[-.01,0]],
            transition_width_m=.002,longitudinal_axis=axis,minimum_stretch=.5)

    def test_analytic_jacobian_matches_independent_difference_and_bound(self):
        rng=np.random.default_rng(140);c=self.cage();p=rng.uniform(-.1,.1,(200,3));p[:,1]=rng.uniform(.01,.99,len(p))
        actual=c.jacobians(p);numerical=np.empty_like(actual)
        for k in range(3):
            d=np.zeros(3);d[k]=1e-7;numerical[:,:,k]=(c.map(p+d)-c.map(p-d))/(2e-7)
        np.testing.assert_allclose(actual,numerical,atol=2e-9,rtol=2e-8)
        self.assertGreaterEqual(np.linalg.det(actual).min(),c.evidence()['continuous_determinant_lower_bound'])

    def test_anchor_curve_and_outside_identity_bands_stay_literal(self):
        c=self.cage();y=c.ay;a=c.av
        points=np.c_[a[:,0],y,a[:,1]];np.testing.assert_array_equal(c.map(points),points)
        outside=np.array([[.2,-1,-.3],[-.5,2,.8]])
        np.testing.assert_array_equal(c.map(outside),outside)
        np.testing.assert_array_equal(c.jacobians(outside),np.tile(np.eye(3),(2,1,1)))

    def test_other_longitudinal_axis_and_inverse_transpose_normal(self):
        c=self.cage(2);p=np.array([[.04,-.05,.44],[.03,.02,.73]])
        j=c.jacobians(p);normal=np.array([.3,.5,.8]);normal/=np.linalg.norm(normal)
        tangent=np.cross(normal,[0,0,1.]);mapped_tangent=np.einsum('nij,j->ni',j,tangent)
        mapped_normal=np.linalg.solve(j.transpose(0,2,1),np.tile(normal,(len(p),1))[...,None])[...,0]
        np.testing.assert_allclose(np.einsum('ni,ni->n',mapped_tangent,mapped_normal),0,atol=1e-15)
        np.testing.assert_array_equal(c.map(p)[:,2],p[:,2])

    def test_invalid_domains_and_nonpositive_stretch_are_rejected(self):
        for knots,negative in [([0,0],[[1,1],[1,1]]),([0,1],[[1,1],[0,1]]),([0,1],[[1,1],[float('nan'),1]])]:
            with self.assertRaises(ValueError):AnchoredAxialCage(knots,negative,[[1,1],[1,1]],[0,1],[[0,0],[0,0]])

if __name__=='__main__':unittest.main()
