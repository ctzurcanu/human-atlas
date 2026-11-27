import unittest
import numpy as np
from geometry_tools.algorithms.convex_core_transport import ConvexCoreTransport,swept_face_separation


class ConvexCoreTransportTests(unittest.TestCase):
    def field(self):
        return ConvexCoreTransport(np.r_[np.eye(3),-np.eye(3)],np.full(6,-1.),np.array([.1,.02,0.]),.2,2e-5,64)

    def test_flat_core_is_bitwise_rigid_and_remote_points_are_exact(self):
        f=self.field();p=np.array([[0.,0.,0.],[1.,.2,.1],[4.,0.,0.]])
        normal=np.tile([0.,1.,0.],(3,1));q,n,s=f.apply(p,normal)
        np.testing.assert_array_equal(q[:2],p[:2]+f.translation)
        np.testing.assert_array_equal(q[2],p[2]);np.testing.assert_array_equal(n,normal)
        self.assertGreater(s['stepDeterminantLowerBound'],.9)

    def test_gradient_matches_derivative_and_covectors_stay_unit(self):
        f=self.field();p=np.array([[1.1,.1,.1]]);_,g,_=f.field(p);eps=1e-7
        numeric=np.array([(f.field(p+eps*a)[0]-f.field(p-eps*a)[0])/(2*eps)for a in np.eye(3)]).ravel()
        np.testing.assert_allclose(g[0],numeric,rtol=2e-7,atol=1e-9)
        _,n,s=f.apply(p,np.array([[1.,1.,0.]])/np.sqrt(2))
        self.assertGreater(s['stepDeterminantMinimumObserved'],s['stepDeterminantLowerBound'])
        self.assertAlmostEqual(float(np.linalg.norm(n[0])),1.)

    def test_filled_triangle_swept_bound_accounts_for_translation(self):
        f=self.field();t=np.array([[[1.3,0,0],[1.3,.2,0],[1.3,0,.2]]])
        gap,_=swept_face_separation(t,f.normals,f.offsets,f.translation)
        self.assertAlmostEqual(float(gap[0]),.2)

    def test_uncertified_step_count_is_rejected(self):
        with self.assertRaises(ValueError):
            ConvexCoreTransport(np.r_[np.eye(3),-np.eye(3)],np.full(6,-1.),np.array([1.,0.,0.]),.01,steps=1)


if __name__=='__main__':unittest.main()
