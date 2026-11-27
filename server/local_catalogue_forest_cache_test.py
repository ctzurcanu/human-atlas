"""Unchanged-model reuse must not hide changed inputs or broken sidecars."""
import importlib.util,json,pathlib,tempfile,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('catalogue_revision',ROOT/'scripts/local_catalogue_revision.py');revision=importlib.util.module_from_spec(spec);spec.loader.exec_module(revision)
class ForestCacheTest(unittest.TestCase):
 def fixture(self,root):
  folder=root/'reports/structural-dependencies/current';folder.mkdir(parents=True)
  for name in ['scripts/build-spatial-dependency-forest.py','scripts/export-structural-relation-edges.mjs','scripts/atlas_mesh_geometry.py','geometry_tools/model.py','geometry_tools/ta98_authority.py','geometry_tools/source_port_dependencies.py','app/data/ta98-metadata.json']:
   path=root/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text('source')
  parent=root/'app/data/ta98-parents.json';parent.write_text('{}')
  (folder/'female-spatial-index.npz').write_bytes(b'index');(folder/'female-anatomical-edges.json').write_text(json.dumps(dict(resolverInputs={'app/data/ta98-parents.json':revision.SHA(parent)})))
  (folder/'female-forest.json').write_text(json.dumps(dict(sourceSha256='current',spatialIndex=dict(file='reports/structural-dependencies/current/female-spatial-index.npz',sha256=revision.SHA(folder/'female-spatial-index.npz')),anatomicalRelationEvidence=dict(sha256=revision.SHA(folder/'female-anatomical-edges.json')))))
  return folder
 def test_unchanged_complete_fresh_model_is_reused(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);self.fixture(root);self.assertTrue(revision.reusable_unchanged_forest('female','current','current',root));self.assertFalse(revision.reusable_unchanged_forest('female','old','current',root))
 def test_corrupt_sidecar_is_not_reused(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);folder=self.fixture(root);(folder/'female-spatial-index.npz').write_bytes(b'changed');self.assertFalse(revision.reusable_unchanged_forest('female','current','current',root))
 def test_changed_relation_code_is_not_reused(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);folder=self.fixture(root);path=root/'app/relations.ts';path.write_text('changed');import os;stamp=(folder/'female-forest.json').stat().st_mtime_ns+1000000000;os.utime(path,ns=(stamp,stamp));self.assertFalse(revision.reusable_unchanged_forest('female','current','current',root))
 def test_changed_parent_json_is_not_hidden_by_unchanged_model_or_mtime(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);self.fixture(root);path=root/'app/data/ta98-parents.json';stamp=path.stat().st_mtime_ns;path.write_text('{"new-native-selector":"A12.2.01.102"}');import os;os.utime(path,ns=(stamp,stamp));self.assertFalse(revision.reusable_unchanged_forest('female','current','current',root))
 def test_old_unbound_relation_snapshot_is_not_reused(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);folder=self.fixture(root);sidecar=folder/'female-anatomical-edges.json';sidecar.write_text('{}');forest=folder/'female-forest.json';value=json.loads(forest.read_text());value['anatomicalRelationEvidence']['sha256']=revision.SHA(sidecar);forest.write_text(json.dumps(value));self.assertFalse(revision.reusable_unchanged_forest('female','current','current',root))
 def test_new_source_port_ledger_invalidates_a_preexisting_male_forest(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);folder=self.fixture(root)
   for suffix in ['forest.json','anatomical-edges.json']:(folder/('male-'+suffix)).write_bytes((folder/('female-'+suffix)).read_bytes())
   ledger=root/'.local-models/geometry-workbench/biodigital-male-circulation-v1/source-port-dependencies.json';ledger.parent.mkdir(parents=True);ledger.write_text('{}')
   self.assertFalse(revision.reusable_unchanged_forest('male','current','current',root))
 def test_source_port_hash_change_is_not_hidden_by_an_unchanged_mtime(self):
  with tempfile.TemporaryDirectory()as d:
   root=pathlib.Path(d);folder=self.fixture(root)
   for suffix in ['forest.json','anatomical-edges.json']:(folder/('male-'+suffix)).write_bytes((folder/('female-'+suffix)).read_bytes())
   ledger=root/'.local-models/geometry-workbench/biodigital-male-circulation-v1/source-port-dependencies.json';ledger.parent.mkdir(parents=True);ledger.write_text('{}')
   path=folder/'male-forest.json';value=json.loads(path.read_text());value['sourcePortDependencyEvidence']={'available':True,'sha256':revision.SHA(ledger)};path.write_text(json.dumps(value));self.assertTrue(revision.reusable_unchanged_forest('male','current','current',root))
   stamp=ledger.stat().st_mtime_ns;ledger.write_text('{"new-port":true}');import os;os.utime(ledger,ns=(stamp,stamp));self.assertFalse(revision.reusable_unchanged_forest('male','current','current',root))
if __name__=='__main__':unittest.main()
