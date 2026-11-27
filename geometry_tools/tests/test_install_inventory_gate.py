import copy
import json
import sqlite3
import unittest
from pathlib import Path
from unittest.mock import patch

from geometry_tools import install, inventory_validation
from geometry_tools.tests import test_inventory_validation as inventory_fixtures

digest=inventory_fixtures.digest


class InstallInventoryGate(unittest.TestCase):
    def setUp(self):
        self.fixture=inventory_fixtures.AuthoritativeInventoryValidation('test_full_inventory_allows_truthful_missing_local_preview')
        self.fixture.setUp();self.addCleanup(self.fixture.doCleanups)
        self.fixture.catalogue.write_text(json.dumps(dict(parts=[],concepts=[],chunks=[])))
        self.fixture.inventory['source_catalogue_sha256']=digest(self.fixture.catalogue)

    def inventory_file(self,inventory=None):
        p=self.fixture.root/'inventory.json';p.write_text(json.dumps(inventory or self.fixture.inventory));return p

    def test_prepare_rejects_omitted_authority_rows_before_any_geometry_or_db_work(self):
        f=self.fixture
        current=f.root/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
        current.parent.mkdir(parents=True);current.write_bytes(f.catalogue.read_bytes())
        proposed=f.root/'proposed.json';proposed.write_bytes(current.read_bytes())
        inv=copy.deepcopy(f.inventory);inv['terms'].pop();inv['term_count']=3
        spec=dict(model='local-ta98',expected_catalogue_sha256=digest(current),proposed_catalogue=str(proposed),
                  changed_parts=[],inventory=str(self.inventory_file(inv)))
        actual=inventory_validation.validate_authoritative_inventory
        def validate(value,**kwargs):
            return actual(value,**kwargs,authority_path=f.authority,corrections_path=f.corrections)
        with patch.object(install,'ROOT',f.root),patch.object(inventory_validation,'validate_authoritative_inventory',side_effect=validate),patch.object(install,'AtlasReader')as reader:
            with self.assertRaisesRegex(ValueError,'full declared authority closure'):
                install.prepare(None,spec,f.root/'output')
            reader.assert_not_called()
        self.assertFalse((f.root/'output').exists())

    def proposal(self,with_validation=True):
        f=self.fixture;inv=self.inventory_file();proposed=f.root/'proposed.json';proposed.write_bytes(f.catalogue.read_bytes())
        p=dict(model='local-ta98',current=str(f.catalogue),expected_catalogue_sha256=digest(f.catalogue),
               proposed_catalogue=str(proposed),proposed_sha256=digest(proposed),inventory=str(inv),inventory_sha256=digest(inv),
               job_id='fixture-reviewed-job',bank_bindings={})
        if with_validation:
            p['inventory_validation']=f.validate();p['inventory_bindings']=p['inventory_validation']['input_bindings']
        path=f.root/'install.json';path.write_text(json.dumps(p));return path

    def test_apply_requires_new_integrity_evidence_even_for_reviewed_job(self):
        proposal=self.proposal(False)
        with patch.object(install,'job_row',return_value={'state':'reviewed','id':'fixture-reviewed-job'}),patch.object(install,'artifacts',return_value=[]):
            with self.assertRaisesRegex(ValueError,'inventory integrity gate; prepare'):
                install.apply(None,proposal,'never-installed')

    def test_apply_rejects_authority_changed_after_preparation(self):
        proposal=self.proposal()
        with sqlite3.connect(self.fixture.authority)as conn:
            conn.execute("UPDATE ta98 SET name_en='later authority edit' WHERE id='A12.2.01.101'")
        with patch.object(install,'job_row',return_value={'state':'reviewed','id':'fixture-reviewed-job'}),patch.object(install,'artifacts',return_value=[]):
            with self.assertRaisesRegex(ValueError,'Bound input/artifact changed'):
                install.apply(None,proposal,'never-installed')


if __name__=='__main__':unittest.main()
