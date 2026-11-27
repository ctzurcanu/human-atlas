import unittest
import numpy as np

from geometry_tools.algorithms.joint_kernel_clearance import GaussianTranslationBasis,propose_joint_step


class JointKernelClearanceTests(unittest.TestCase):
    def basis(self):
        return GaussianTranslationBasis([[-.004,0,0],[.004,0,0],[0,0,.004],[0,0,-.004]],.008,frame='test:metres',free_below=.02,fixed_above=.03)

    def test_fixed_root_and_analytic_jacobian(self):
        b=self.basis();p=np.array([[.001,.01,.002],[0,.04,0]])
        c=np.array([[.0002,0,0],[-.0001,.0001,0],[0,0,.0003],[0,.0002,0]])
        q,J=b.map_and_jacobians(p,c);np.testing.assert_array_equal(q[1],p[1]);np.testing.assert_array_equal(J[1],np.eye(3))
        for axis in range(3):
            h=np.zeros(3);h[axis]=1e-7
            a=b.map_and_jacobians(p+h,c)[0];d=b.map_and_jacobians(p-h,c)[0]
            np.testing.assert_allclose((a-d)/(2e-7),J[:,:,axis],atol=1e-8,rtol=0)

    def test_joint_corridor_improves_violation_without_crossing_other_wall(self):
        b=self.basis();p=np.array([[-.0002,0,0],[.0004,0,.001]])
        # Permitted interval is x in [0, .0006]. Both increasing-clearance
        # normals are required; optimizing just one could cross the other wall.
        gaps=np.r_[p[:,0],.0006-p[:,0]];normals=np.array([[1.,0,0],[1.,0,0],[-1.,0,0],[-1.,0,0]])
        c,proof=propose_joint_step(p,b,gaps,normals,np.tile(np.arange(2),2),maximum_step_m=.0003,wall_seconds=3)
        q,_=b.map_and_jacobians(p,c)
        self.assertTrue(proof['complete']);self.assertGreater(q[0,0],p[0,0]);self.assertTrue((q[:,0]<=.0006).all());self.assertLessEqual(np.linalg.norm(q-p,axis=1).max(),.0003)

    def test_stationary_boundary_keeps_value_and_derivative(self):
        rng=np.random.default_rng(47)
        b=GaussianTranslationBasis(rng.uniform(-.008,.008,(16,3)),.008,frame='test:metres',free_below=.02,fixed_above=.03)
        p=np.array([[.002,0,0],[-.0002,.002,.001]])
        c,proof=propose_joint_step(p,b,[-.0002],[[1.,0,0]],[1],maximum_step_m=.0003,wall_seconds=3,stationary_point_ids=[0])
        q,J=b.map_and_jacobians(p,c)
        np.testing.assert_allclose(q[0],p[0],atol=1e-13,rtol=0);np.testing.assert_allclose(J[0],np.eye(3),atol=1e-10,rtol=0)
        self.assertGreater(q[1,0],p[1,0]);self.assertEqual(proof['stationaryControls'],1)
