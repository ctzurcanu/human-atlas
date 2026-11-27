import unittest
from unittest.mock import patch
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.anchor_path_feasibility import (
    anchor_path_feasibility, mesh_anchor_path_feasibility,
)


class AnchorPathFeasibilityTests(unittest.TestCase):
    def chain(self, target, **kwargs):
        p = np.array([[0., 0., 0.], [.5, 0., 0.], [1., 0., 0.]])
        return anchor_path_feasibility(p, np.array([[0, 1], [1, 2]]),
                                       [0, 2], target,
                                       maximum_edge_extension=1.15, **kwargs)

    def test_straight_chain_cannot_double_its_anchor_separation(self):
        r = self.chain([[0., 0., 0.], [2., 0., 0.]])
        self.assertFalse(r['passedNecessaryBound'])
        self.assertEqual(r['violatingPairs'], 1)
        w = r['largestRequiredExtensionWitness']
        self.assertEqual((w['firstRawId'], w['secondRawId']), (0, 2))
        self.assertEqual(w['originalShortestSourceEdgePath_m'], 1.)
        self.assertEqual(w['desiredAnchorSeparation_m'], 2.)
        self.assertEqual(w['requiredEdgeExtensionLowerBound'], 2.)

    def test_bent_target_pass_is_explicitly_not_a_surface_feasibility_claim(self):
        r = self.chain([[0., 0., 0.], [.5, .5, 0.]])
        self.assertTrue(r['passedNecessaryBound'])
        self.assertFalse(r['sufficientFeasibilityEstablished'])
        self.assertFalse(r['anatomicalAcceptance'])

    def test_shortest_source_edge_path_is_stronger_than_a_longer_route(self):
        p = np.array([[0., 0., 0.], [1., 0., 0.], [1., 1., 0.]])
        target = [[0., 0., 0.], [1.8, 0., 0.]]
        long = anchor_path_feasibility(p, np.array([[0, 1], [1, 2]]), [0, 2], target, maximum_edge_extension=1.15)
        direct = mesh_anchor_path_feasibility(Mesh(p, np.array([[0, 1, 2]]), 'test:metres'), [0, 2], target, maximum_edge_extension=1.15)
        self.assertTrue(long['passedNecessaryBound'])
        self.assertFalse(direct['passedNecessaryBound'])
        self.assertAlmostEqual(direct['largestRequiredExtensionWitness']['originalShortestSourceEdgePath_m'], np.sqrt(2))

    def test_disconnected_pairs_do_not_imply_impossibility_or_feasibility(self):
        r = anchor_path_feasibility(np.array([[0., 0., 0.], [1., 0., 0.], [2., 0., 0.]]), np.array([[0, 1]]), [0, 2], [[0., 0., 0.], [100., 0., 0.]], maximum_edge_extension=1.15)
        self.assertTrue(r['passedNecessaryBound'])
        self.assertEqual(r['checkedConnectedPairs'], 0)
        self.assertEqual(r['disconnectedPairsIgnored'], 1)
        self.assertIsNone(r['largestRequiredExtensionWitness'])
        self.assertFalse(r['sufficientFeasibilityEstablished'])

    def test_exact_aliases_couple_targets_and_reduce_zero_length_edges(self):
        p = np.array([[0., 0., 0.], [1., 0., 0.], [0., 0., 0.]])
        r = anchor_path_feasibility(p, np.array([[0, 2], [2, 1]]), [2, 0, 1], [[0., 0., 0.], [0., 0., 0.], [1., 0., 0.]], maximum_edge_extension=1.)
        self.assertTrue(r['passedNecessaryBound'])
        self.assertEqual(r['physicalAnchors'], 2)
        self.assertEqual(r['declaredRawAnchors'], 3)
        self.assertEqual(r['collapsedExactAliasEdges'], 1)
        self.assertEqual(r['largestRequiredExtensionWitness']['firstRawId'], 2)
        with self.assertRaisesRegex(ValueError, 'inconsistent target'):
            anchor_path_feasibility(p, np.array([[0, 2], [2, 1]]), [0, 2], [[0., 0., 0.], [1e-15, 0., 0.]], maximum_edge_extension=1.)

    def test_nearby_points_are_not_welded(self):
        p = np.array([[0., 0., 0.], [1e-9, 0., 0.]])
        r = anchor_path_feasibility(p, np.array([[0, 1]]), [0, 1], [[0., 0., 0.], [1e-8, 0., 0.]], maximum_edge_extension=1.15)
        self.assertFalse(r['passedNecessaryBound'])
        self.assertEqual(r['sourcePhysicalControls'], 2)
        self.assertEqual(r['collapsedExactAliasEdges'], 0)

    def test_small_positive_edges_do_not_underflow_to_an_alias(self):
        p = np.array([[0., 0., 0.], [1e-250, 0., 0.]])
        r = anchor_path_feasibility(p, np.array([[0, 1]]), [0, 1], [[0., 0., 0.], [1e-249, 0., 0.]], maximum_edge_extension=1.15)
        self.assertFalse(r['passedNecessaryBound'])
        self.assertGreater(r['largestRequiredExtensionWitness']['originalShortestSourceEdgePath_m'], 0.)

    def test_deadline_cannot_return_a_partial_passing_result(self):
        with patch('geometry_tools.algorithms.anchor_path_feasibility.time.monotonic', side_effect=[0., 10.]):
            with self.assertRaises(TimeoutError):self.chain([[0., 0., 0.], [1., 0., 0.]])

    def test_limits_and_invalid_anchor_inputs_are_rejected(self):
        p = np.array([[0., 0., 0.], [1., 0., 0.]])
        edges = np.array([[0, 1]])
        for factor in [.9, 0., float('inf'), float('nan')]:
            with self.subTest(factor=factor), self.assertRaises(ValueError):
                anchor_path_feasibility(p, edges, [0, 1], p, maximum_edge_extension=factor)
        for wall in [0., 5.01, float('nan')]:
            with self.subTest(wall=wall), self.assertRaises(ValueError):
                anchor_path_feasibility(p, edges, [0, 1], p, maximum_edge_extension=1., wall_seconds=wall)
        with self.assertRaises(ValueError):anchor_path_feasibility(p, edges.astype(float), [0, 1], p, maximum_edge_extension=1.)
        with self.assertRaises(ValueError):anchor_path_feasibility(p, edges, [0, 2], p, maximum_edge_extension=1.)
        with self.assertRaises(ValueError):anchor_path_feasibility(p, edges, list(range(129)), np.zeros((129, 3)), maximum_edge_extension=1.)
        with self.assertRaises(ValueError):anchor_path_feasibility(p, edges, [0, 1], [[0., 0., 0.], [float('nan'), 0., 0.]], maximum_edge_extension=1.)

    def test_source_and_target_arrays_stay_literal(self):
        p = np.array([[0., 0., 0.], [.5, .2, 0.], [1., 0., 0.]])
        edges = np.array([[0, 1], [1, 2]])
        raw = np.array([0, 2]);target=p[raw].copy()
        before=[a.copy() for a in [p, edges, raw, target]]
        anchor_path_feasibility(p, edges, raw, target, maximum_edge_extension=1.15)
        for a, original in zip([p, edges, raw, target], before):np.testing.assert_array_equal(a, original)


if __name__ == '__main__':unittest.main()
