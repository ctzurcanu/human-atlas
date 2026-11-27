import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from geometry_tools.algorithms.revision_mouth_ledger import inherit_or_retain_mouth_ledger


class RevisionMouthLedgerTests(unittest.TestCase):
    def fixture(self, root):
        old,new,review=[root/name for name in ['old.json','new.json','review.json']]
        first=dict(first='a',second='b',firstRing=[0,1,2],secondRing=[3,4,5])
        second=dict(first='b',second='c',firstRing=[6,7,8],secondRing=[0,1,2])
        old.write_text(json.dumps(dict(interfaces=[first],frame='metres')))
        new.write_text(json.dumps(dict(interfaces=[first,second],frame='metres')))
        sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
        review.write_text(json.dumps(dict(complete=True,passed=True,catalogueSha256='catalogue',
            sourceReportSha256='report',previousLedgerSha256=sha(old),sourcePortLedgerSha256=sha(new),
            addedInterfaces=[second])))
        return old,new,review

    def call(self, paths):
        return inherit_or_retain_mouth_ledger(*paths,catalogue_sha='catalogue',source_report_sha='report')

    def test_reviewed_connection_survives_refresh(self):
        with tempfile.TemporaryDirectory()as d:
            paths=self.fixture(Path(d));before=paths[1].read_bytes()
            result=self.call(paths)
            self.assertEqual(result['added'],1)
            self.assertEqual(paths[1].read_bytes(),before)

    def test_unknown_or_stale_extension_is_never_silently_replaced(self):
        with tempfile.TemporaryDirectory()as d:
            paths=self.fixture(Path(d));paths[2].unlink();before=paths[1].read_bytes()
            with self.assertRaises(ValueError):self.call(paths)
            self.assertEqual(paths[1].read_bytes(),before)
            paths=self.fixture(Path(d));paths[1].write_text(paths[1].read_text()+' ')
            with self.assertRaises(ValueError):self.call(paths)

    def test_missing_ledger_inherits_exact_previous_bytes(self):
        with tempfile.TemporaryDirectory()as d:
            paths=self.fixture(Path(d));paths[1].unlink()
            self.assertEqual(self.call(paths)['added'],0)
            self.assertEqual(paths[1].read_bytes(),paths[0].read_bytes())


if __name__=='__main__':unittest.main()
