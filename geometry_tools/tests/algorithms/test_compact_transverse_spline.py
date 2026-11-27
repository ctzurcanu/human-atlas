import unittest
import numpy as np
from geometry_tools.algorithms.compact_transverse_spline import CompactTransverseSpline


class TransverseSplineTests(unittest.TestCase):
    def fixture(self):
        knots=np.r_[[-.1]*4,np.arange(-.075,.1,.025),[.1]*4]
        field=CompactTransverseSpline([dict(y_m=-.1,x_m=-.002,z_m=.001),dict(y_m=.1,x_m=.002,z_m=-.001)],knots,core_xz_m=[.006,.006],support_xz_m=[.02,.02])
        c=np.zeros((field.count,2));c[3:-3]=np.random.default_rng(84).uniform(-.0006,.0006,(field.count-6,2))
        return field,c

    def test_jacobian_calibre_and_literal_fixed_support(self):
        field,c=self.fixture();p=np.random.default_rng(1304).uniform([-.03,-.12,-.03],[.03,.12,.03],(300,3));q,j,proof=field.apply(p,c);fd=np.zeros_like(j)
        for axis in range(3):
            step=np.eye(3)[axis]*1e-7;fd[:,:,axis]=(field.apply(p+step,c)[0]-field.apply(p-step,c)[0])/(2e-7)
        np.testing.assert_allclose(j,fd,atol=6e-8);self.assertTrue(np.all(np.linalg.det(j)>=proof['analytic_determinant_lower_bound']-1e-12));np.testing.assert_array_equal(q[:,1],p[:,1])
        outside=np.array([[.1,0,0],[0,.2,0],[0,0,.1]]);np.testing.assert_array_equal(field.apply(outside,c)[0],outside)
        ring=np.array([[-.003,.01,-.003],[.003,.01,-.003],[.003,.01,.003],[-.003,.01,.003]])
        changed=field.apply(ring,c)[0];np.testing.assert_allclose(changed-changed[0],ring-ring[0],atol=1e-17)

    def test_source_sample_operator_matches_actual_triangle_barycentrics(self):
        field,c=self.fixture();v=np.array([[-.002,-.01,0],[.002,.02,.001],[0,.04,-.001]]);w=np.array([.2,.3,.5]);changed=field.apply(v,c)[0]
        sample=v.T@w;actual=changed.T@w;displacement=(w@field.design(v))@c
        np.testing.assert_allclose(actual,sample+np.array([displacement[0],0,displacement[1]]),atol=1e-17)

    def test_endpoint_motion_and_nonpositive_certificate_are_rejected(self):
        field,c=self.fixture();c[0,0]=.0001
        with self.assertRaises(ValueError):field.apply(np.zeros((1,3)),c)
        c[:]=0;c[4,0]=.1
        with self.assertRaises(ValueError):field.apply(np.zeros((1,3)),c)
