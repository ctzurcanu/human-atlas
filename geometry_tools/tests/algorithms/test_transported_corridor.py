import unittest
import numpy as np
from geometry_tools.algorithms.compact_positive_fields import transported_corridor_translation


class TransportedCorridorTests(unittest.TestCase):
    def setUp(self):
        self.kw = dict(centerline_stations=[dict(y_m=.55, x_m=.02, z_m=0.), dict(y_m=.72, x_m=.024, z_m=.01)],
            core_axes_m=[.004, .01, .004], support_axes_m=[.009, .1, .009],
            translation_m=[.011, 0., 0.], y_core_interval_m=[.605, .66], y_support_interval_m=[.555, .715], steps=48)

    def test_full_core_profiles_are_rigid_and_y_fixed(self):
        theta = np.arange(24)*2*np.pi/24
        for y in [.58, .625, .69]:
            cx = .02+(y-.55)/.17*.004
            cz = (y-.55)/.17*.01
            points = np.c_[cx+.002*np.cos(theta), np.full(len(theta), y), cz+.002*np.sin(theta)]
            q, j, proof = transported_corridor_translation(points, **self.kw)
            np.testing.assert_allclose(q-q[0], points-points[0], atol=2e-15)
            np.testing.assert_array_equal(q[:, 1], points[:, 1])
            self.assertGreater(proof['analytic_determinant_lower_bound'], 0)
            np.testing.assert_allclose(np.linalg.det(j), 1., atol=2e-13)

    def test_chain_rule_matches_finite_differences_and_determinants_positive(self):
        rng = np.random.default_rng(4)
        p = np.c_[rng.uniform(.010, .05, 100), rng.uniform(.56, .71, 100), rng.uniform(-.006, .015, 100)]
        q, j, proof = transported_corridor_translation(p, **self.kw)
        h = 1e-8
        for axis in range(3):
            delta = np.eye(3)[axis]*h
            hi = transported_corridor_translation(p+delta, **self.kw)[0]
            lo = transported_corridor_translation(p-delta, **self.kw)[0]
            np.testing.assert_allclose((hi-lo)/(2*h), j[:, :, axis], atol=2e-6, rtol=2e-6)
        self.assertGreater(float(np.linalg.det(j).min()), proof['analytic_determinant_lower_bound']*.999)

    def test_outside_support_and_neighbouring_deep_vessel_fixed_exactly(self):
        p = np.array([[.06, .64, 0.], [.02, .5, 0.], [-.02, .64, 0.]])
        q, j, _ = transported_corridor_translation(p, **self.kw)
        np.testing.assert_array_equal(p, q)
        np.testing.assert_array_equal(j, np.broadcast_to(np.eye(3), j.shape))

    def test_rejects_uncertified_single_step(self):
        with self.assertRaises(ValueError):
            transported_corridor_translation(np.zeros((1, 3)), **dict(self.kw, steps=1))


if __name__ == '__main__':
    unittest.main()
