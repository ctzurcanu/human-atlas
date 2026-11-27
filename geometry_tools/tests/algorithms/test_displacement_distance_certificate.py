import importlib.util
from pathlib import Path
import unittest
import numpy as np

P=Path(__file__).resolve().parents[3]/'server/audit-biodigital-displacement-skin-guard.py'
S=importlib.util.spec_from_file_location('bounded_displacement_skin_guard',P);M=importlib.util.module_from_spec(S);S.loader.exec_module(M)

class DisplacementDistanceCertificateTests(unittest.TestCase):
    def test_triangle_vertex_error_bounds_every_barycentric_correspondence(self):
        vertices=np.array([[0.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        native=vertices+np.array([[-.125,0,0],[.0625,0,0],[0,.0625,0]])
        weights=np.array([[1.,0,0],[0,1.,0],[0,0,1.],[.25,.25,.5],[.125,.25,.625]])
        epsilon=np.linalg.norm(native-vertices,axis=1).max()
        self.assertTrue(np.all(np.linalg.norm(weights@native-weights@vertices,axis=1)<=epsilon))

    def test_conditional_lipschitz_bound_in_analytical_worst_direction(self):
        # Original plane x=0, native plane x=-epsilon; query rounds +delta.
        true_distance=np.array([.5,.75]);epsilon=.125;delta=.0625;eta=.03125
        native_computed=true_distance+epsilon+delta+eta
        lower,certified=M.conditional_distance_lower_bound(native_computed,epsilon,delta,eta)
        self.assertTrue(certified);np.testing.assert_array_equal(lower,true_distance)

    def test_unproved_native_witness_can_never_prune(self):
        # A literal witness on x=10 is 9.5 away, while the actual mesh also
        # contains x=0 only .5 away. A witness is an upper bound, not a lower.
        witness=np.array([9.5]);lower,certified=M.conditional_distance_lower_bound(witness,0.,0.,None)
        self.assertFalse(certified);np.testing.assert_array_equal(lower,[0.])
        self.assertGreater(float(witness[0]),.5)

    def test_negative_error_rejected_and_lower_clamped(self):
        with self.assertRaises(ValueError):M.conditional_distance_lower_bound([1.],-.1,[0.],0.)
        with self.assertRaises(ValueError):M.conditional_distance_lower_bound([1.],0.,[0.],-.1)
        lower,ok=M.conditional_distance_lower_bound([.1],.2,[.1],0.)
        self.assertTrue(ok);np.testing.assert_array_equal(lower,[0.])

if __name__=='__main__':unittest.main()
