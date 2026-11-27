import unittest
import numpy as np
from geometry_tools.algorithms.common_origin_encoding import encode_positions


class CommonOriginEncodingTests(unittest.TestCase):
    def test_shared_atlas_origins_decode_identically(self):
        rng = np.random.default_rng(92)
        p = rng.uniform([-.05, .86, -.04], [.05, .94, .04], (400, 3))
        origins = np.array([[0, .9140625, -.0244140625], [0, .9140625, 0], [0, 0, 0]])
        masks = rng.integers(1, 8, len(p), dtype=np.uint64)
        q, evidence = encode_positions(p, masks, origins)
        for i, origin in enumerate(origins):
            used = (masks & (1 << i)) != 0
            np.testing.assert_array_equal((q[used] - origin).astype('<f4').astype(float) + origin, q[used])
        self.assertTrue(evidence['all_owning_decoders_exact'])
        self.assertLessEqual(evidence['maximum_world_rounding_error_m'], evidence['maximum_derived_rounding_bound_m'])

    def test_source_interface_already_in_both_lattices_is_literal(self):
        p = np.array([[0, .890625, -.015625], [.03125, .875, 0.]])
        q, _ = encode_positions(p, np.array([3, 3]), [[0, .9140625, -.0244140625], [0, 0, 0]])
        np.testing.assert_array_equal(q, p)

    def test_single_bank_keeps_its_finer_local_precision(self):
        p = np.array([[.001, .914063001, -.024413876]])
        origin = np.array([[0, .9140625, -.0244140625]])
        q, _ = encode_positions(p, [1], origin)
        np.testing.assert_array_equal(q, (p-origin).astype('<f4').astype(float)+origin)

    def test_incompatible_origins_are_rejected(self):
        with self.assertRaisesRegex(ValueError, 'incompatible'):
            encode_positions([[.5, .5, .5]], [3], [[0, 0, 0], [1e-9, 0, 0]])

    def test_invalid_masks_and_nonfinite_input_are_rejected(self):
        with self.assertRaises(ValueError):
            encode_positions([[0, 0, 0]], [2], [[0, 0, 0]])
        with self.assertRaises(ValueError):
            encode_positions([[np.nan, 0, 0]], [1], [[0, 0, 0]])

    def test_unused_controls_stay_literal(self):
        p = np.array([[.123456789, .987654321, .011223344]])
        q, _ = encode_positions(p, [0], [[0, 0, 0]])
        np.testing.assert_array_equal(q, p)


if __name__ == '__main__':
    unittest.main()
