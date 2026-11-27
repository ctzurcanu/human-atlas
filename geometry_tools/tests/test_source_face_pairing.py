import unittest
import tempfile
from pathlib import Path
import numpy as np
from geometry_tools.algorithms.source_face_pairing import remap_source_pairs
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.course_io import sha


class SourceFacePairingTests(unittest.TestCase):
    def setUp(self):
        self.faces = np.array([[0,1,2],[0,2,3],[3,2,4],[3,4,5],[5,4,6],[5,6,7]])

    def test_removed_faces_and_control_compaction_replay(self):
        kept = np.array([4,5,0,1])
        current = (self.faces[kept]+10)
        current = np.vstack([current, [20,21,22]])
        pairs, proof = remap_source_pairs(self.faces, current, kept, np.arange(6).reshape(-1,2))
        np.testing.assert_array_equal(pairs, [[2,3],[0,1]])
        self.assertEqual(proof['constructedFaces'], 1)
        self.assertFalse(proof['constructedQuadPairingInferred'])

    def test_stale_or_reversed_face_mapping_is_rejected(self):
        bad = self.faces.copy()
        bad[1] = bad[1, ::-1]
        with self.assertRaisesRegex(ValueError, 'control map'):
            remap_source_pairs(self.faces, bad, np.arange(6), np.arange(6).reshape(-1,2))

    def test_aliasing_original_controls_is_rejected(self):
        bad = self.faces.copy()
        bad[bad == 7] = 6
        with self.assertRaisesRegex(ValueError, 'control map'):
            remap_source_pairs(self.faces, bad, np.arange(6), np.arange(6).reshape(-1,2))

    def test_missing_or_stale_bank_cannot_silently_infer_pairs(self):
        with self.assertRaisesRegex(ValueError, 'Missing explicit'):
            load_source_pairing(None, self.faces)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'bank.npz'
            np.savez(path, source_faces=self.faces,
                     retained_source_triangle_pairs=np.arange(6).reshape(-1,2))
            binding=dict(path=str(path),sha256=sha(path))
            self.assertEqual(load_source_pairing(binding,self.faces).shape,(3,2))
            with self.assertRaisesRegex(ValueError, 'Stale source'):
                load_source_pairing(binding,self.faces[::-1])
