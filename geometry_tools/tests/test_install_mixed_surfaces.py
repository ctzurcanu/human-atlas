import unittest
import tempfile,os
from pathlib import Path
from unittest.mock import patch
import numpy as np
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.common import digest_file
from geometry_tools.audit import audit
from geometry_tools.install import emitted_requires_closed


class MixedSurfaceInstallation(unittest.TestCase):
    def test_open_component_view_does_not_relax_finite_wall_gate(self):
        vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float)
        closed=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        open_view=Mesh(vertices,closed[:-1],'fixture:metres')
        spec=dict(require_closed=False,closed_parts=['wall'])
        self.assertTrue(audit(open_view,require_closed=emitted_requires_closed(spec,'view'))['passed'])
        self.assertFalse(audit(open_view,require_closed=emitted_requires_closed(spec,'wall'))['passed'])
        wall=Mesh(vertices,closed,'fixture:metres')
        self.assertTrue(audit(wall,require_closed=emitted_requires_closed(spec,'wall'))['passed'])

    def test_unspecified_policy_still_requires_closed_geometry(self):
        self.assertTrue(emitted_requires_closed({},'wall'))

    def test_shared_bank_hash_reuse_rechecks_a_changed_file(self):
        with tempfile.TemporaryDirectory()as directory:
            path=Path(directory)/'bank.bin';path.write_bytes(b'first');reader=AtlasReader()
            with patch('geometry_tools.model.digest_file',wraps=digest_file)as hashed:
                expected=reader.file_sha(path)
                for _ in range(100):self.assertEqual(reader.file_sha(path),expected)
                self.assertEqual(hashed.call_count,1)
                stamp=path.stat().st_mtime_ns;path.write_bytes(b'other')
                os.utime(path,ns=(stamp+1_000_000,stamp+1_000_000))
                self.assertNotEqual(reader.file_sha(path),expected)
                self.assertEqual(hashed.call_count,2)


if __name__=='__main__':unittest.main()
