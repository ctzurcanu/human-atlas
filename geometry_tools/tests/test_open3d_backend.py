"""Native acceleration must preserve literal witnesses and fitting safeguards."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex, diagnostics, icp_transform, run_query
from geometry_tools.spatial import BVH
from geometry_tools.audit import winding_contains
from geometry_tools.fit import transform_mesh
from test_core import tetra


def organ_fixture():
    m = tetra()
    return Mesh(m.vertices * [.08, .12, .06], m.faces, m.frame,
                {'uvs': np.arange(8, dtype=np.float32).reshape(4, 2),
                 'normalSigns': np.array([1, -1, 1, 1], dtype=np.int8)})


def icp_controls():
    return dict(target_frame='fixture', max_correspondence_m=.015,
                voxel_sizes_m=[.004, .002, .001], max_correction_m=.02,
                max_rotation_degrees=10, sample_count=4000, minimum_fitness=.95)


class Open3DTests(unittest.TestCase):
    def test_batch_witnesses_match_literal_triangles_and_float64_search(self):
        m = organ_fixture(); origin = np.array([10000., -20000, 30000])
        m = Mesh(m.vertices + origin, m.faces, m.frame, m.attributes)
        q = np.random.default_rng(77).uniform([-.02, -.02, -.02], [.1, .14, .08], (100, 3)) + origin
        index = SurfaceIndex(m); fast, evidence = index.closest_points(q)
        certified, proof = index.closest_points(q, certify=True)
        reference, _ = BVH(m).closest_points(q, backend='bvh')
        np.testing.assert_allclose(certified['distance_m'], reference['distance_m'], atol=1e-12, rtol=1e-12)
        actual = np.einsum('ni,nij->nj', fast['barycentric'], m.triangles[fast['face_ids']])
        np.testing.assert_allclose(actual, fast['points'], atol=1e-10, rtol=0)
        np.testing.assert_allclose(fast['barycentric'].sum(1), 1, atol=1e-14)
        self.assertTrue(np.all(fast['distance_m'] >= reference['distance_m'] - 1e-12))
        self.assertFalse(evidence['nearest_triangle_verified']); self.assertTrue(proof['nearest_triangle_verified'])
        self.assertEqual(index.closest_points(np.empty((0, 3)))[0]['points'].shape, (0, 3))

    def test_bad_native_seed_is_corrected_by_full_bvh_search(self):
        m = organ_fixture(); q = np.array([[.02, .02, -.01]])
        expected = BVH(m).closest_point(q[0])
        actual = BVH(m).closest_point(q[0], seed_face=3)
        self.assertAlmostEqual(actual['distance_m'], expected['distance_m'], places=14)

    def test_thin_triangle_rounding_uses_exact_fallback(self):
        v = np.array([[0., 0, 0], [1, 0, 0], [1, 1e-12, 0], [0, 1, 0]])
        # Normalization translates tiny separation onto a large Float32 word.
        m = Mesh(v, np.array([[0, 1, 2], [0, 2, 3]]), 'fixture')
        with self.assertRaisesRegex(ValueError, 'collapses'):
            SurfaceIndex(m)
        q = np.tile([.5, .3, .01], (40, 1))
        a, info = BVH(m).closest_points(q, backend='auto')
        self.assertEqual(info['backend'], 'bvh'); self.assertTrue(info['nearest_triangle_verified'])
        np.testing.assert_allclose(a['distance_m'], .01)

    def test_rays_retain_face_ids_and_explicit_misses(self):
        m = Mesh(np.array([[0., 0, 0], [1, 0, 0], [0, 1, 0]]), np.array([[0, 1, 2]]), 'fixture')
        result, info = SurfaceIndex(m).cast_rays([[.2, .3, 1], [2, 2, 1]], [[0, 0, -5], [0, 0, -1]])
        np.testing.assert_array_equal(result['hit'], [True, False]); np.testing.assert_array_equal(result['face_ids'], [0, -1])
        self.assertAlmostEqual(result['distance_m'][0], 1, places=7)
        self.assertLess(info['maximum_ray_residual_m'], 1e-7)
        self.assertFalse(info['nearest_hit_verified'])
        with self.assertRaises(ValueError): SurfaceIndex(m).cast_rays([[0, 0, 0]], [[0, 0, 0]])

    def test_native_diagnostics_keep_seams_and_do_not_mutate_sources(self):
        m = organ_fixture()
        # Every face has its own attribute seam; weld equal coordinates only.
        split = Mesh(m.triangles.reshape(-1, 3), np.arange(12).reshape(4, 3), m.frame)
        before = (split.geometry_sha, split.attributes_sha)
        info = diagnostics(split)
        self.assertTrue(info['vertex_manifold']); self.assertEqual(info['self_intersection_pairs'], 0)
        self.assertFalse(info['acceptance_gate']); self.assertEqual(before, (split.geometry_sha, split.attributes_sha))
        crossing = Mesh(np.array([[-1., 0, 0], [1, 0, 0], [0, 1, 0], [0, .2, -1], [0, .2, 1], [0, -.5, 0]]), np.arange(6).reshape(2, 3), m.frame)
        self.assertEqual(diagnostics(crossing)['self_intersection_pairs'], 1)

    def test_certified_batch_containment_matches_existing_winding(self):
        m = tetra(); q = np.random.default_rng(19).uniform(-.2, 1.2, (40, 3)); q[0] = 0
        a = winding_contains(m, q, distance_backend='bvh')
        b = winding_contains(m, q, distance_backend='open3d')
        self.assertEqual(a['inside'], b['inside']); self.assertEqual(a['boundary_ambiguous'], b['boundary_ambiguous'])
        self.assertTrue(b['distance_query']['nearest_triangle_verified'])

    def test_full_defect_export_does_not_truncate_repair_input(self):
        pair=np.array([[-1.,0,0],[1,0,0],[0,1,0],[0,.2,-1],[0,.2,1],[0,-.5,0]])
        vertices=np.concatenate([pair+[5*i,0,0] for i in range(23)])
        mesh=Mesh(vertices,np.arange(len(vertices)).reshape(-1,3),'fixture')
        before=mesh.geometry_sha
        info=diagnostics(mesh,include_pairs=True,maximum_export_pairs=23)
        self.assertEqual(info['self_intersection_pairs'],23)
        self.assertEqual(len(info['source_face_pair_examples']),20)
        self.assertTrue(info['source_face_pairs_complete'])
        self.assertEqual({tuple(sorted(p)) for p in info['source_face_pairs']},{(2*i,2*i+1) for i in range(23)})
        self.assertEqual(mesh.geometry_sha,before)
        with self.assertRaisesRegex(ValueError,'export budget'):
            diagnostics(mesh,include_pairs=True,maximum_export_pairs=22)
        with self.assertRaises(ValueError):
            diagnostics(mesh,self_intersections=False,include_pairs=True)

    def test_native_near_coplanar_false_positive_never_becomes_a_defect(self):
        # Real capsule pair: near-parallel triangles separated by 34 micrometres.
        v = np.array([[.05968920700252056, 1.1788093550130725, -.009751982055604458],
                      [.06045939959585667, 1.177353659644723, -.010102082509547472],
                      [.060113002080470324, 1.1773926774039865, -.010092698503285646],
                      [.06043651863001287, 1.1772470995783806, -.01012771064415574],
                      [.06008202629163861, 1.1772470995783806, -.01012771064415574],
                      [.06010606023482978, 1.177360052242875, -.010100544895976782]])
        info = diagnostics(Mesh(v, np.arange(6).reshape(2, 3), 'fixture'))
        self.assertEqual(info['self_intersection_pairs'], 0)
        self.assertEqual(info['native_candidate_pairs'], info['native_candidates_rejected_by_source_predicates'])
        if info['native_pair_replay_examples']:
            self.assertGreater(info['native_pair_replay_examples'][0]['distance_m'], 3e-5)

    def test_icp_refinement_preserves_scale_faces_attributes_and_limits(self):
        m = organ_fixture(); angle = np.radians(3)
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        initial = np.diag([.8, .8, .8, 1.]); expected = initial.copy()
        expected[:3, :3] = .8 * rotation; expected[:3, 3] = [.002, -.001, .003]
        target = transform_mesh(m, expected)
        matrix, evidence = icp_transform(m, target, initial, **icp_controls())
        actual = transform_mesh(m, matrix)
        np.testing.assert_allclose(actual.vertices, target.vertices, atol=8e-5, rtol=0)
        np.testing.assert_array_equal(actual.faces, m.faces); np.testing.assert_array_equal(actual.attributes['uvs'], m.attributes['uvs'])
        np.testing.assert_array_equal(actual.attributes['normalSigns'], m.attributes['normalSigns'])
        np.testing.assert_allclose(matrix[:3, :3].T @ matrix[:3, :3], np.eye(3) * .64, atol=1e-12)
        self.assertLess(evidence['inlier_rmse_m'], evidence['initial_inlier_rmse_m'])
        with self.assertRaisesRegex(ValueError, 'bounds'):
            icp_transform(m, target, initial, **{**icp_controls(), 'max_correction_m': 1e-5})
        with self.assertRaises(ValueError):
            icp_transform(m, target, np.diag([-1., 1, 1, 1]), **icp_controls())
        with self.assertRaises(ValueError):
            icp_transform(m, target, initial, **{**icp_controls(), 'target_frame': 'different'})

    def test_cli_frame_and_overwrite_guards(self):
        from argparse import Namespace
        with tempfile.TemporaryDirectory() as td:
            td = Path(td); m = organ_fixture(); m.save(td / 'mesh.npz')
            np.savez(td / 'queries.npz', points=np.array([[.01, .01, -.01]]), frame=np.array('fixture'))
            args = Namespace(mesh=str(td / 'mesh.npz'), queries=str(td / 'queries.npz'),
                             output=str(td / 'result.npz'), model=None, frame=None,
                             rays=False, certify=True, threads=2)
            info = run_query(args, None)
            self.assertTrue(info['nearest_triangle_verified']); self.assertTrue((td / 'result.npz').exists())
            with self.assertRaises(FileExistsError): run_query(args, None)
            args.output = str(td / 'other.npz')
            np.savez(td / 'bad-frame.npz', points=np.array([[0., 0, 0]]), frame=np.array('other'))
            args.queries = str(td / 'bad-frame.npz')
            with self.assertRaisesRegex(ValueError, 'frame'): run_query(args, None)
            with patch('geometry_tools.open3d_backend.available', return_value=False):
                _, fallback = BVH(m).closest_points(np.zeros((40, 3)), backend='auto')
                self.assertEqual(fallback['backend'], 'bvh')


if __name__ == '__main__':
    unittest.main()
