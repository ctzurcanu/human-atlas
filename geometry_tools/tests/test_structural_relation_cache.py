import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location('relation_cache', Path(__file__).resolve().parents[2] / 'scripts/structural_relation_cache.py')
cache = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cache)


class RelationshipCacheTests(unittest.TestCase):
    def test_only_hierarchy_checksum_may_change(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / cache.PARENT_INDEX
            path.parent.mkdir(parents=True)
            original = {'source': 'authority', 'parents': {'child': 'parent'}, 'hierarchySha256': 'before'}
            raw = json.dumps(original).encode()
            path.write_bytes(raw)
            consumer = root / 'consumer.ts'
            consumer.write_text('read .parents')
            prepared = {'resolverInputs': {cache.PARENT_INDEX: cache.sha(raw), 'consumer.ts': cache.sha(consumer.read_bytes())}}
            cache.bind_parent_snapshot(root, prepared)
            self.assertTrue(cache.resolver_inputs_match(root, prepared))
            path.write_text(json.dumps({**original, 'hierarchySha256': 'after'}))
            # The synthetic consumer stands in for an explicitly audited graph.
            with patch.object(cache,'AUDITED_PARENT_ONLY_INPUT_SIGNATURE',cache.parent_only_input_signature(prepared['resolverInputs'])):
                self.assertTrue(cache.resolver_inputs_match(root, prepared))
                extra=root/'new-consumer.ts';extra.write_text('uses hierarchy checksum')
                prepared['resolverInputs']['new-consumer.ts']=cache.sha(extra.read_bytes())
                self.assertFalse(cache.resolver_inputs_match(root,prepared))
                del prepared['resolverInputs']['new-consumer.ts']
                path.write_text(json.dumps({**original, 'parents': {'child': 'different'}}))
                self.assertFalse(cache.resolver_inputs_match(root, prepared))
                path.write_text(json.dumps({**original, 'hierarchySha256': 'after'}))
                consumer.write_text('changed consumer')
                self.assertFalse(cache.resolver_inputs_match(root, prepared))

    def test_unbound_snapshot_cannot_authorize_reuse(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            path = root / cache.PARENT_INDEX
            path.parent.mkdir(parents=True)
            raw = json.dumps({'source': 'authority', 'parents': {}, 'hierarchySha256': 'new'}).encode()
            path.write_bytes(raw)
            prepared = {'resolverInputs': {cache.PARENT_INDEX: 'not-the-original-hash'}}
            cache.bind_parent_snapshot(root, prepared)
            self.assertNotIn('resolverInputSnapshots', prepared)
            self.assertFalse(cache.resolver_inputs_match(root, prepared))


if __name__ == '__main__':
    unittest.main()
