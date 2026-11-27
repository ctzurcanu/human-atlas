import unittest
import numpy as np
from geometry_tools.algorithms.patch_metric_constraints import PatchMetricConstraints


class PatchMetricConstraintTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0,0],[.0002,0,0],[.0001,.00003,0],[.0001,.00002,.0001]])+[.03,1.72,.07]
        self.f=np.array([[0,1,2],[0,2,3]])

    def test_analytic_jacobian_matches_centered_directional_difference(self):
        c=PatchMetricConstraints(self.p,self.f,[1,2,3],scale_m=1e-5,strain=.075,rotation_degrees=45)
        z=np.array([.12,-.08,.05,-.07,.1,.04,.03,.02,-.04]);d=np.array([1.,2,3,4,5,6,7,8,9]);d/=np.linalg.norm(d)
        value,jac=c.evaluate(z);h=1e-3
        fd=(c.evaluate(z+h*d)[0]-c.evaluate(z-h*d)[0])/(2*h)
        np.testing.assert_allclose(jac@d,fd,rtol=2e-4,atol=2e-7)

    def test_source_feasible_and_fixed_controls_literal(self):
        c=PatchMetricConstraints(self.p,self.f,[1,2],scale_m=1e-5,strain=.075,rotation_degrees=45)
        values,_=c.evaluate(np.zeros(6));self.assertTrue(np.all(values>0))
        q=c.construct(np.ones(6)*.1);np.testing.assert_array_equal(q[[0,3]],self.p[[0,3]])

    def test_large_material_stretch_is_rejected(self):
        c=PatchMetricConstraints(self.p,self.f,[1,2,3],scale_m=1e-5,strain=.075,rotation_degrees=45)
        z=((self.p[1:]-self.p[0])*.3/1e-5).ravel()
        self.assertLess(c.evaluate(z)[0].min(),0)
