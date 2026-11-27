import importlib.util
import unittest

spec=importlib.util.spec_from_file_location('blood_ownership_refresh','scripts/refresh-blood-review-ownership.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class BloodMetadataPreservationTest(unittest.TestCase):
 def test_geometry_history_can_only_append_its_actual_revision(self):
  old={'parts':[],'chunks':[],'concepts':[{'id':'parent'}],'bloodProfileClearanceHistory':[{'revision':'old'}]}
  new={**old,'bloodProfileClearanceHistory':old['bloodProfileClearanceHistory']+[{'revision':'new','sourceCatalogueSha256':'actual-source'}]}
  module.verify_preserved_metadata(old,new,'actual-source','new')
  with self.assertRaises(AssertionError):module.verify_preserved_metadata(old,new,'stale-source','new')
  with self.assertRaises(AssertionError):module.verify_preserved_metadata(old,new,'actual-source','other')

 def test_cannot_rewrite_history_or_change_anatomical_identity(self):
  old={'concepts':[{'id':'parent'}],'bloodProfileClearanceHistory':[{'revision':'old'}]}
  new={**old,'bloodProfileClearanceHistory':[{'revision':'new','sourceCatalogueSha256':'source'}]}
  with self.assertRaises(AssertionError):module.verify_preserved_metadata(old,new,'source','new')
  new={**old,'concepts':[{'id':'changed'}]}
  with self.assertRaises(AssertionError):module.verify_preserved_metadata(old,new,'source','new')


if __name__=='__main__':unittest.main()
