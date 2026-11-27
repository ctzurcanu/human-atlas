import unittest
import numpy as np
from geometry_tools.algorithms.hinge_isometry import project_hinge_isometry,intrinsic_facet_strains


class HingeIsometryTests(unittest.TestCase):
    def test_rotation_preserves_edges_and_intrinsic_metric(self):
        p=np.array([[0.,0,0],[2,0,0],[.8,.1,0]])
        q=p.copy();q[2]=[.9,.1,.2]
        q[2]=project_hinge_isometry(p,q,2,[0,1])
        np.testing.assert_allclose(np.linalg.norm(q[2]-q[:2],axis=1),np.linalg.norm(p[2]-p[:2],axis=1),atol=1e-13)
        np.testing.assert_allclose(intrinsic_facet_strains(p,q,np.array([[0,1,2]])),0,atol=1e-13)
        np.testing.assert_array_equal(q[:2],p[:2])

    def test_translation_covariance(self):
        p=np.array([[0.,0,0],[1,0,0],[.4,.2,0]]);q=p.copy();q[2,2]=.3
        shift=np.array([.5,1.8,-.2])
        np.testing.assert_allclose(project_hinge_isometry(p+shift,q+shift,2,[0,1])-shift,project_hinge_isometry(p,q,2,[0,1]),atol=1e-14)

    def test_moving_or_degenerate_hinge_rejected(self):
        p=np.array([[0.,0,0],[1,0,0],[.4,.2,0]]);q=p.copy();q[0,0]=.1
        with self.assertRaises(ValueError):project_hinge_isometry(p,q,2,[0,1])
        p[1]=p[0]
        with self.assertRaises(ValueError):project_hinge_isometry(p,p.copy(),2,[0,1])
