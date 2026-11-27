import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.identity_target_sections import solve_identity_target_sections
from geometry_tools.algorithms.hybrid_section_clearance import solve_hybrid_sections


class IdentityTargetSections(unittest.TestCase):
    def test_slack_elimination_matches_original_constrained_optimization(self):
        S=sparse.csr_matrix([[1.,-1.]])
        A=sparse.csr_matrix([[-1.,0.,0.,0.,0.,0.],[0.,0.,0.,0.,-1.,0.]])
        b=np.array([-.0001,-.0002]);I=sparse.eye(6,format='csr')
        B=sparse.vstack([I,-I]);target=np.array([.0008,.0002,-.0005,-.0006,.0007,.0003]);c=np.r_[target,-target]
        args=(S,A,b,B,c)
        old,p0=solve_hybrid_sections(*args,maximum_motion_m=.001,time_limit_seconds=2.)
        new,p1=solve_identity_target_sections(*args,maximum_motion_m=.001,time_limit_seconds=2.)
        self.assertIsNotNone(old);self.assertIsNotNone(new)
        # Independent analytic optimum: x0 and y1 are active; the other
        # coordinates solve the two-cell ridge normal equations. The older
        # slack solver deliberately stops at a certified approximate point.
        w=10000.;expected=np.empty((2,3))
        expected[0,0]=.0001;expected[1,1]=.0002
        expected[1,0]=(w*target[3]+expected[0,0])/(w+2)
        expected[0,1]=(w*target[1]+expected[1,1])/(w+2)
        expected[:,2]=np.linalg.solve(np.array([[w+2,-1.],[-1.,w+2]]),w*target[[2,5]])
        np.testing.assert_allclose(new,expected,atol=1e-10,rtol=1e-7)
        def objective(d):
            x=d/.001
            return .5*(np.linalg.norm(x)**2+np.linalg.norm(S@x)**2+w*np.linalg.norm(x.ravel()-target/.001)**2)
        gap=objective(old)-objective(new)
        self.assertGreaterEqual(gap,-1e-6)
        self.assertLessEqual(gap,p0['terminationHistory'][-1]['absoluteGap']+1e-6)
        self.assertLessEqual(float(np.max(b-A@new.ravel())),1e-10)
        self.assertEqual(p1['eliminatedSlackVariables'],12)

    def test_unpaired_target_is_not_silently_reinterpreted(self):
        with self.assertRaises(ValueError):
            solve_identity_target_sections(sparse.eye(1),sparse.csr_matrix((0,3)),np.array([]),
                sparse.eye(3),np.zeros(3),maximum_motion_m=.001)


if __name__=='__main__':unittest.main()
