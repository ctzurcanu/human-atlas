import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from geometry_tools import review
from geometry_tools.common import digest_file
from geometry_tools.model import Mesh


class ReviewMetadataTests(unittest.TestCase):
    def test_request_retains_part_identity_and_only_declared_framing_markers(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            frame = 'model:test:assembled-metres'
            vertices = np.array([[0., 0., 0.], [.001, 0., 0.], [0., .001, 0.]])
            faces = np.array([[0, 1, 2]])
            entries = []
            cases = [
                ('TEST:anatomical-part', False, False),
                ('TEST:framing-from-entry', True, False),
                ('TEST:framing-from-source', False, True),
            ]
            for index, (part_id, entry_marker, source_marker) in enumerate(cases):
                path = root / f'source-{index}.npz'
                source = {'nonanatomical_framing_only': True} if source_marker else {}
                Mesh(vertices + index * .01, faces, frame, source=source).save(path)
                entry = {'mesh': str(path), 'name': part_id, 'part_id': part_id}
                if entry_marker:
                    entry['nonanatomical_framing_only'] = True
                entries.append(entry)

            # The temporary queue cannot dispatch a request to the live bridge.
            with patch.object(review, 'QUEUE', root / 'queue'), patch.object(review, 'ROOT', root):
                result = review.prepare(None, {'meshes': entries, 'wall_seconds': 1})

            request_path = Path(result['folder']) / 'request.json'
            request = json.loads(request_path.read_text())
            rows = request['meshes']
            self.assertEqual([row['part_id'] for row in rows], [case[0] for case in cases])
            self.assertNotIn('nonanatomical_framing_only', rows[0])
            self.assertIs(rows[1]['nonanatomical_framing_only'], True)
            self.assertIs(rows[2]['nonanatomical_framing_only'], True)
            self.assertEqual(request['frame'], frame)
            self.assertEqual(digest_file(request_path), result['request_sha256'])

            pending = root / 'queue' / 'pending' / f"{result['request_id']}.json"
            receipt = json.loads(pending.read_text())
            self.assertEqual(Path(receipt['request']), request_path)
            self.assertEqual(receipt['sha256'], result['request_sha256'])
            for row in rows:
                exported = Mesh.load(row['path'])
                self.assertEqual(exported.frame, frame)
                self.assertEqual(exported.geometry_sha, row['geometry_sha256'])
                self.assertEqual(digest_file(row['path']), row['sha256'])


if __name__ == '__main__':
    unittest.main()
