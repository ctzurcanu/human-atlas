import copy
import hashlib
import sqlite3
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from geometry_tools.inventory_validation import validate_authoritative_inventory


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


class AuthoritativeInventoryValidation(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name);self.authority=self.root/'authority.sqlite'
        self.corrections=self.root/'corrections.py';self.corrections.write_text('# fixture corrections\n')
        self.catalogue=self.root/'current.json';self.catalogue.write_text('{"revision":"current"}\n')
        with sqlite3.connect(self.authority)as conn:
            conn.execute('CREATE TABLE ta98(id TEXT PRIMARY KEY,parent_id TEXT,name_en TEXT,name_la TEXT,type_of_entity TEXT,female_gender INTEGER,male_gender INTEGER,immaterial INTEGER)')
            conn.executemany('INSERT INTO ta98 VALUES(?,?,?,?,?,?,?,?)',[
                ('A12.2.01.001',None,'root','root','material entity',0,0,0),
                ('A12.2.01.101','A12.2.01.001','branch','branch','material entity',0,0,0),
                ('A12.2.01.102','A12.2.01.101','opening','opening','immaterial entity',0,0,1),
                ('A12.2.01.103','A12.2.01.001','female branch','branch','material entity',1,0,0),
                ('A12.2.02.001',None,'other root','other','material entity',0,0,0),
                ('A12.2.02.002','A12.2.02.001','other child','other child','material entity',0,0,0)])
            conn.row_factory=sqlite3.Row;self.rows={r['id']:dict(r)for r in conn.execute('SELECT * FROM ta98')}
        self.inventory=dict(schema='ta98-geometry-inventory/v1',model='local-ta98',model_sex='male',
            root_term='A12.2.01.001',authority_path=str(self.authority),authority_sha256=digest(self.authority),
            authority_corrections_path=str(self.corrections),authority_corrections_sha256=digest(self.corrections),
            source_catalogue_sha256=digest(self.catalogue),term_count=4,applicable_term_count=3,
            terms=[self.term(code)for code in list(self.rows)[:4]])

    def term(self,code):
        return dict(authority=copy.deepcopy(self.rows[code]),representation_class='opening'if code.endswith('.102')else'tissue-region',
                    applicable_to_model=code!='A12.2.01.103',authority_sex='female'if code.endswith('.103')else None,
                    source_meshes=[],missing_distinct_representation=True,distinct_representation_verified=False)

    def validate(self,inventory=None):
        return validate_authoritative_inventory(inventory or self.inventory,model='local-ta98',current_catalogue=self.catalogue,
                                                authority_path=self.authority,corrections_path=self.corrections)

    def test_full_inventory_allows_truthful_missing_local_preview(self):
        result=self.validate()
        self.assertTrue(result['authority_inventory_complete']);self.assertEqual(result['term_count'],4)
        self.assertEqual(result['applicable_term_count'],3);self.assertEqual(result['declared_missing_representation_count'],3)
        self.assertEqual(result['unverified_representation_count'],3)
        self.assertFalse(result['anatomical_coverage_complete']);self.assertFalse(result['anatomical_coverage_evaluated'])
        self.assertFalse(result['ordinary_view_representation_coverage_evaluated'])

    def test_blood_can_be_classified_as_luminal_content_without_a_wall_claim(self):
        self.inventory['terms'][1]['representation_class']='luminal-content'
        result=self.validate()
        self.assertEqual(self.inventory['terms'][1]['representation_class'],'luminal-content')
        self.assertEqual(result['unverified_representation_count'],3)
        self.assertFalse(result['anatomical_coverage_complete'])

    def test_omitted_duplicate_and_falsely_counted_terms_fail(self):
        changed=copy.deepcopy(self.inventory);changed['terms'].pop();changed['term_count']=3
        with self.assertRaisesRegex(ValueError,'full declared authority closure'):self.validate(changed)
        changed=copy.deepcopy(self.inventory);changed['terms'][1]=changed['terms'][0]
        with self.assertRaisesRegex(ValueError,'duplicate'):self.validate(changed)
        changed=copy.deepcopy(self.inventory);changed['term_count']=1
        with self.assertRaisesRegex(ValueError,'term_count'):self.validate(changed)

    def test_authority_codes_parents_and_names_are_not_caller_metadata(self):
        for field,value in [('id','A12.2.09.999'),('parent_id','A12.2.02.001'),('name_en','invented anatomy')]:
            changed=copy.deepcopy(self.inventory);changed['terms'][1]['authority'][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):self.validate(changed)

    def test_missing_unclassified_and_invented_classes_fail(self):
        for value in [None,'unclassified','anything']:
            changed=copy.deepcopy(self.inventory)
            if value is None:changed['terms'][0].pop('representation_class')
            else:changed['terms'][0]['representation_class']=value
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,'legal classified'):self.validate(changed)

    def test_current_catalogue_binding_supports_existing_shapes_but_rejects_stale(self):
        for key,value in [('current_catalogue_sha256',digest(self.catalogue)),('source_catalogue',dict(path=str(self.catalogue),sha256=digest(self.catalogue)))]:
            changed=copy.deepcopy(self.inventory);changed.pop('source_catalogue_sha256');changed[key]=value
            self.assertTrue(self.validate(changed)['passed'])
        for value in [None,'0'*64]:
            changed=copy.deepcopy(self.inventory)
            if value is None:changed.pop('source_catalogue_sha256')
            else:changed['source_catalogue_sha256']=value
            with self.assertRaisesRegex(ValueError,'exact current catalogue'):self.validate(changed)

    def test_sex_applicability_cannot_reduce_the_denominator(self):
        changed=copy.deepcopy(self.inventory);changed['terms'][1]['applicable_to_model']=False
        with self.assertRaisesRegex(ValueError,'applicability'):self.validate(changed)
        changed=copy.deepcopy(self.inventory);changed['applicable_term_count']=1
        with self.assertRaisesRegex(ValueError,'applicable_term_count'):self.validate(changed)

    def test_multiple_roots_and_explicit_context_keep_canonical_parents(self):
        changed=copy.deepcopy(self.inventory);changed['additional_context_root']='A12.2.02.001'
        changed['terms'] += [self.term(c)for c in ['A12.2.02.001','A12.2.02.002']]
        changed['term_count']=6;changed['applicable_term_count']=5
        self.assertEqual(self.validate(changed)['term_count'],6)
        changed=copy.deepcopy(self.inventory)
        with patch.dict('geometry_tools.inventory_validation.CONTEXT_REFERENCES',{'A12.2.01.001':{'A12.2.02.002':'primary fixture'}}):
            with self.assertRaisesRegex(ValueError,'full declared authority closure'):self.validate(changed)
            changed['terms'].append(self.term('A12.2.02.002'));changed['term_count']=5;changed['applicable_term_count']=4
            self.assertTrue(self.validate(changed)['passed'])

    def test_hidden_sources_and_candidate_ids_do_not_prove_coverage(self):
        changed=copy.deepcopy(self.inventory)
        for row in changed['terms']:
            row.update(source_meshes=[dict(part_id='hidden-donor',hidden=True,suppressed=True)],candidate_part_ids=['missing-candidate'],
                       missing_distinct_representation=False,distinct_representation_verified=False)
        changed['terms'][0]['planned_representation']={'boundaryProvisional':True}
        result=self.validate(changed)
        self.assertEqual(result['unverified_representation_count'],3)
        self.assertEqual(result['explicitly_provisional_count'],1)
        self.assertEqual(result['declared_verified_representation_count'],0)
        self.assertFalse(result['anatomical_coverage_complete'])

    def test_changed_authority_or_corrections_rejected_at_replay(self):
        first=self.validate();self.assertIn(str(self.authority.resolve()),first['input_bindings'])
        with sqlite3.connect(self.authority)as conn:conn.execute("UPDATE ta98 SET name_en='changed' WHERE id='A12.2.01.101'")
        with self.assertRaisesRegex(ValueError,'actual current authoritative cache'):self.validate()
        self.inventory['authority_sha256']=digest(self.authority)
        self.corrections.write_text('# changed corrections\n')
        with self.assertRaisesRegex(ValueError,'stale authority corrections'):self.validate()


if __name__=='__main__':unittest.main()
