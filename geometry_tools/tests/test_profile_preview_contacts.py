"""Acceptance-policy guards for a partial source-preserving assembly preview."""
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

from geometry_tools.course_io import sha

spec=importlib.util.spec_from_file_location('profile_preview','scripts/install-profile-clearance-preview.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class ProfilePreviewContactsTest(unittest.TestCase):
 def fixture(self,path):
  source=path/'source.json';source.write_text('{}');atlas=path/'atlas.json';atlas.write_text('{}')
  expected={'artery':{'sha256':'artery-sha'},'vein':{'sha256':'vein-sha'}};paths=[]
  for pid,peer,system in [('artery','vein','venous'),('vein','artery','arterial')]:
   r=dict(complete=True,candidateContextConsistent=True,contextStart=0,contextLimit=3,activeContextCount=3,candidateReportSha256=sha(source),catalogueSha256=sha(atlas),partId=pid,candidateSha256=expected[pid]['sha256'],rows=[dict(partId=peer,system=system,before=3,after=2),dict(partId='bone',system='skeletal',before=20,after=19),dict(partId='skin',system='integumentary',before=8,after=7)])
   p=path/(pid+'.json');p.write_text(json.dumps(r));paths.append(p)
  return paths,source,atlas,expected

 def test_counts_mutual_blood_contacts_once(self):
  with tempfile.TemporaryDirectory()as tmp:
   args=self.fixture(Path(tmp));r=module.contact_evidence(*args)
   self.assertEqual(r['criticalContacts']['blood'],[3,2]);self.assertEqual(r['newCriticalObjectPairs'],0)

 def test_rejects_new_critical_pair_despite_lower_total(self):
  with tempfile.TemporaryDirectory()as tmp:
   paths,source,atlas,expected=self.fixture(Path(tmp));r=json.loads(paths[0].read_text());r['rows'].append(dict(partId='other-bone',system='skeletal',before=0,after=1));paths[0].write_text(json.dumps(r))
   with self.assertRaisesRegex(AssertionError,'New critical'):module.contact_evidence(paths,source,atlas,expected)

 def test_rejects_omitted_owner_stale_report_and_worse_system(self):
  with tempfile.TemporaryDirectory()as tmp:
   paths,source,atlas,expected=self.fixture(Path(tmp))
   with self.assertRaises(AssertionError):module.contact_evidence(paths[:1],source,atlas,expected)
   r=json.loads(paths[0].read_text());r['rows'][1]['after']=30;paths[0].write_text(json.dumps(r))
   with self.assertRaisesRegex(AssertionError,'Worse whole-assembly skeletal'):module.contact_evidence(paths,source,atlas,expected)
   source.write_text('{"changed":true}')
   with self.assertRaises(AssertionError):module.contact_evidence(paths,source,atlas,expected)


if __name__=='__main__':unittest.main()
