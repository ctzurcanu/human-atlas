import unittest
import numpy as np
from geometry_tools.algorithms.constrained_radial_fit import (
    fit_radial_shell,evaluate_radial_fit,RadialFitInfeasible)


class ConstrainedRadialFitTests(unittest.TestCase):
    def test_uniform_shell_matches_known_linear_map(self):
        r=np.linspace(.01,.015,31);target=.95*r
        fit=fit_radial_shell(r,target,target)
        np.testing.assert_allclose(evaluate_radial_fit(fit,r),target,atol=2e-11,rtol=0)
        self.assertTrue(fit['analyticProof']['jacobianDeterminantLower']>0)
        self.assertLessEqual(fit['elapsedSeconds'],3.)

    def test_nonuniform_fit_when_no_uniform_scale_exists(self):
        r=np.array([.01,.012,.015]);lower=np.array([.0099,.0116,.014]);upper=lower.copy()
        self.assertGreater((lower/r).max(),(upper/r).min())
        fit=fit_radial_shell(r,lower,upper)
        np.testing.assert_allclose(evaluate_radial_fit(fit,r),lower,atol=2e-11,rtol=0)
        self.assertGreater(fit['derivativeCoefficientL1Variation'],0)
        self.assertEqual(len(fit['nativeStages']),2)

    def test_contradictory_host_constraints_are_rejected(self):
        with self.assertRaises(RadialFitInfeasible):
            fit_radial_shell([.01,.01],[.0095,.0099],[.0096,.0100])

    def test_required_extension_beyond_derivative_bound_is_rejected(self):
        with self.assertRaises(RadialFitInfeasible):
            fit_radial_shell([.01,.0101],[.01,.0103],[.01,.0103])

    def test_derivative_and_ratio_proof_prevents_between_sample_overshoot(self):
        r=np.linspace(.01,.015,11);fit=fit_radial_shell(r,.85*r,.9*r)
        x=np.linspace(*fit['domain_m'],3001);y=evaluate_radial_fit(fit,x);d=evaluate_radial_fit(fit,x,derivative=1)
        self.assertGreaterEqual(float((y/x).min()),.8-2e-10)
        self.assertLessEqual(float((y/x).max()),1.05+2e-10)
        self.assertGreaterEqual(float(d.min()),.75-2e-10)
        self.assertLessEqual(float(d.max()),1.5+2e-10)
        self.assertGreater(fit['analyticProof']['jacobianDeterminantLower'],0)

    def test_no_extrapolation_or_nonfinite_input(self):
        fit=fit_radial_shell([.01,.015],[.009,.014],[.010,.015])
        with self.assertRaises(ValueError):evaluate_radial_fit(fit,[.009])
        with self.assertRaises(ValueError):evaluate_radial_fit(fit,[np.nan])
        with self.assertRaises(ValueError):fit_radial_shell([.01,np.nan],[.009,.014],[.010,.015])

    def test_resources_and_input_arrays_remain_bounded(self):
        r=np.array([.01,.015]);low=np.array([.009,.014]);high=np.array([.010,.015]);copies=[x.copy()for x in [r,low,high]]
        fit_radial_shell(r,low,high)
        for a,b in zip([r,low,high],copies):np.testing.assert_array_equal(a,b)
        with self.assertRaises(ValueError):fit_radial_shell(r,low,high,max_knots=65)
        with self.assertRaises(ValueError):fit_radial_shell(r,low,high,wall_seconds=0)


if __name__=='__main__':unittest.main()
