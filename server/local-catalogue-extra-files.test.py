"""Sidecar payloads must share the catalogue transaction and its rollback."""
import gzip,importlib.util,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('revision_extra',Path(__file__).resolve().parents[1]/'scripts/local_catalogue_revision.py')
revision=importlib.util.module_from_spec(spec);spec.loader.exec_module(revision)

class ExtraFilesTest(unittest.TestCase):
    def test_path_hash_and_output_guards(self):
        with tempfile.TemporaryDirectory()as d,patch.object(revision,'ROOT',Path(d).resolve()):
            root=Path(d).resolve();p=root/'labels.json';p.write_bytes(b'old');sha=revision.SHA(p)
            changes,expected=revision.guarded_extra_changes({p:b'new'},{p:sha},[])
            self.assertEqual(changes,{p:b'new'});self.assertEqual(expected,{p:sha})
            for edits,hashes,reserved in [({p:b'new'},{},[]),({p:b'new'},{p:'0'*64},[]),({p:b'new'},{p:sha},[p]),({root.parent/'outside':b'x'},{root.parent/'outside':'0'*64},[])]:
                with self.assertRaises(ValueError):revision.guarded_extra_changes(edits,hashes,reserved)
            alias=root/'alias';alias.symlink_to(p)
            with self.assertRaises(ValueError):revision.guarded_extra_changes({p:b'a',alias:b'b'},{p:sha,alias:sha},[])

    def test_failed_derived_build_restores_sidecar_and_catalogues(self):
        with tempfile.TemporaryDirectory()as d,patch.object(revision,'ROOT',Path(d).resolve()):
            root=Path(d).resolve();routes=[];original={}
            for sex in ['male','female']:
                cp=root/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json';cp.parent.mkdir(parents=True)
                cp.write_bytes(revision.ENC(dict(parts=[dict(id='tube',chunk=0)],chunks=[dict(url='/source.bin')])));original[cp]=cp.read_bytes()
                z=Path(str(cp)+'.gz');z.write_bytes(gzip.compress(cp.read_bytes(),mtime=0));original[z]=z.read_bytes()
                vp=root/f'.local-models/ta98-review/visits/{sex}.json';vp.parent.mkdir(exist_ok=True)
                vp.write_bytes(revision.ENC(dict(bindings=[dict(id='tube')],sourceCatalogueSha256=revision.SHA(cp))));original[vp]=vp.read_bytes()
                routes.append(dict(sex=sex,sourceCatalogueSha256=revision.SHA(cp),routeSha256=revision.SHA(vp)))
            rp=root/'reports/alimentary-visit-paths.json';rp.parent.mkdir();rp.write_bytes(revision.ENC(routes));original[rp]=rp.read_bytes()
            labels=root/'.local-models/ta98-review/bowel-wall-local-atlas/male/cardiac-label-sites.json';labels.write_bytes(b'old labels');original[labels]=labels.read_bytes()
            def failed_build(*args,**kwargs):
                self.assertEqual(labels.read_bytes(),b'new labels')
                raise RuntimeError('simulated derived build failure')
            with patch.object(revision.subprocess,'run',side_effect=failed_build):
                with self.assertRaisesRegex(RuntimeError,'simulated'):
                    revision.install_revision('extra-rollback-test',lambda sex,a,before:a.update(version='new'),extra_changes={labels:b'new labels'},expected_extra_sha256={labels:revision.SHA(labels)})
            for path,data in original.items():self.assertEqual(path.read_bytes(),data,str(path))
            manifest=json.loads((root/'.local-models/ta98-review/extra-rollback-test/manifest.json').read_bytes())
            self.assertIn(str(labels.relative_to(root)),[row['path']for row in manifest])
            self.assertFalse((root/'reports/extra-rollback-test.json').exists())

if __name__=='__main__':unittest.main()
