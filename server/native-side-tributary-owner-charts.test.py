"""Focused regressions for mixed literal/transported charts and stale audits.

These are nonanatomical attribute/receipt fixtures. No wall construction,
packing, source catalogue, atlas or geometry acceptance is performed.
"""
from pathlib import Path
import importlib.util
import json
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file  # noqa: E402
from geometry_tools.model import Mesh  # noqa: E402

module_path = ROOT/'server/pack-native-side-tributary-owners.py'
loader = importlib.util.spec_from_file_location('two_owner_chart_test_subject', module_path)
subject = importlib.util.module_from_spec(loader)
loader.loader.exec_module(subject)


class OwnerChartTests(unittest.TestCase):
    def test_fixed_word_survives_mixed_float_chart_and_flipped_mask_rejects(self):
        with tempfile.TemporaryDirectory(prefix='atlas-owner-charts-') as name:
            base = Path(name)
            vertices = np.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
            faces = np.array([[0, 1, 2]], dtype=np.int64)
            # Renormalizing this legitimate outward word changes its leading
            # component. The fixed corner must retain all three original words.
            raw = np.array([[32760, 100, 200]]*3, dtype=np.int16)
            normals = raw.astype(float)/np.linalg.norm(raw.astype(float), axis=1)[:, None]
            source = base/'source.npz'
            atomic_npz(source, vertices=vertices, faces=faces,
                profiles=vertices[None], profile_vertex_indices=np.array([[0, 1, 2]]),
                attribute_normals=normals, attribute_uvs=np.array([[0., 0.], [1., 0.], [0., 1.]], dtype=np.float32),
                literal_normal_words=raw, literal_normal_mask=np.array([True, False, False]))
            proof_path = base/'transport.json'
            proof = dict(proper_normal_transport=True, normal_sign_flips=0,
                target_source_sha256=digest_file(source),
                literal_normal_words_array='literal_normal_words', literal_normal_mask_array='literal_normal_mask')
            atomic_json(proof_path, proof)
            binding = lambda p: dict(path=str(p), sha256=digest_file(p))
            config = dict(encoding='transported-float64', transport_proof=binding(proof_path),
                literal_normal_words_array='literal_normal_words', literal_normal_mask_array='literal_normal_mask')
            spec = dict(normal_charts=dict(parent=config, tributary=config))
            plan = dict(parent=binding(source), tributary=binding(source))
            sources = subject.source_charts(plan, base/'unused-plan.json', spec, {})
            surfaces = {'collector:outer': faces, 'tributary:outer': faces}
            masks = dict(collector=np.array([True, False]), tributary=np.array([False, True]))
            data = dict(outer_component=np.array([0, 1]), outer_source_face=np.array([0, 0]))
            words, uvs, _, counts = subject.charts(vertices, vertices, surfaces, masks, data, sources)
            np.testing.assert_array_equal(words['collector:outer'][0, 0], raw[0])
            self.assertFalse(np.array_equal(words['collector:outer'][0, 1], raw[1]),
                             'Fixture must distinguish literal words from requantized vectors')
            np.testing.assert_array_equal(uvs['collector:outer'][0],
                                          np.array([[0., 0.], [1., 0.], [0., 1.]], dtype=np.float32))
            self.assertEqual(counts['parent']['transported_chart_literal_raw_words'], 1)
            self.assertEqual(counts['parent']['float64_literal_chart_quantized_once'], 2)
            # Keep the target binding valid while deliberately violating the
            # promised fixed-normal direction: the actual vectors must catch it.
            flipped = base/'flipped-source.npz'
            with np.load(source, allow_pickle=False) as z:
                arrays = {k: z[k].copy() for k in z.files}
            arrays['attribute_normals'][0] *= -1
            atomic_npz(flipped, **arrays)
            bad_proof_path = base/'flipped-transport.json'
            atomic_json(bad_proof_path, dict(proof, target_source_sha256=digest_file(flipped)))
            bad_config = dict(config, transport_proof=binding(bad_proof_path))
            with self.assertRaisesRegex(ValueError, 'does not preserve the literal native word direction'):
                subject.source_charts(dict(parent=binding(flipped), tributary=binding(source)),
                    base/'unused-plan.json', dict(normal_charts=dict(parent=bad_config, tributary=config)), {})

    def test_passed_audit_for_other_actual_geometry_is_rejected(self):
        vertices = np.array([[0., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        faces = np.array([[0, 1, 2]], dtype=np.int64)
        mesh = Mesh(vertices, faces, 'nonanatomical:audit-identity:metres')
        row = dict(passed=True, geometry_sha256=mesh.geometry_sha,
                   attributes_sha256=mesh.attributes_sha, frame=mesh.frame)
        subject.verify_matching_audit(row, mesh, 'Fixture')
        moved = Mesh(vertices+np.array([.001, 0., 0.]), faces, mesh.frame)
        with self.assertRaisesRegex(ValueError, 'different actual geometry'):
            subject.verify_matching_audit(row, moved, 'Fixture')
        with self.assertRaisesRegex(ValueError, 'frame or attributes'):
            subject.verify_matching_audit(dict(row, frame='another-frame'), mesh, 'Fixture')


if __name__ == '__main__':
    unittest.main()
