"""Exercise the metadata builder's actual resolver without executing its writes."""
import ast,collections,json,pathlib,re,subprocess,sys,tempfile,unittest,unicodedata
ROOT=pathlib.Path(__file__).resolve().parents[1]
source=ast.parse((ROOT/'scripts/build-ta98-metadata.py').read_text())
functions={'declared_identity','resolve','matched','record','normalized','cleaned','variants','orthographic'}
module=ast.Module(body=[n for n in source.body if isinstance(n,ast.FunctionDef) and n.name in functions],type_ignores=[])

class CatalogueIdentityTest(unittest.TestCase):
 def setUp(self):
  ids=['A01.1.00.001','A01.1.00.012','A02.5.07.001','A02.5.07.002','A02.5.07.005','A02.5.05.005','A04.3.01.002','A05.6.02.002','A09.1.03.008','A09.1.02.001F']
  terms={key:dict(id=key,source_id=key,name_en=key,name_la=key,entity_id_number=None,fma_id=None,parent_id=None) for key in ids}
  self.env=dict(re=re,unicodedata=unicodedata,terms=terms,reviewed=[])
  for key in ['by_primary','by_synonym','by_base_primary','by_base_synonym','by_qualified','by_fma','by_orthographic']:self.env[key]=collections.defaultdict(set)
  self.env['by_primary'].update({'head':{'A01.1.00.001'},'neck':{'A01.1.00.012'},'superior part':{'A04.3.01.002'},'anterior surface':{'A09.1.03.008'}})
  exec(compile(module,'actual metadata resolver','exec'),self.env)
 def resolve(self,name,cid,system,identity=None):return self.env['resolve'](name,cid,system,identity)
 def test_provisional_regions_keep_the_authored_term_and_scope(self):
  for name,cid,system in [('head','TA98FIBULA:A02.5.07.002','skeletal'),('neck','TA98FIBULA:A02.5.07.005','skeletal'),('anterior surface','TA98PATELLA:A02.5.05.005','skeletal'),('superior part','TA98DUODENUM:A05.6.02.002','digestive')]:
   with self.subTest(cid=cid):
    term=cid.split(':')[-1];r=self.resolve(name,cid,system,{'ta98Term':term,'ta98Kind':'parent'});self.assertEqual((r['_term'],r['_kind'],r['_method']),(term,'parent','catalogue-tag'));self.assertIsNone(r['fma'])
 def test_declared_code_alone_does_not_claim_exact_geometry(self):
  r=self.resolve('head','TA98FIBULA:A02.5.07.002','skeletal');self.assertEqual((r['_term'],r['_kind']),('A02.5.07.002','parent'))
 def test_generic_review_cannot_overwrite_authored_identity(self):
  self.env['reviewed']=[dict(term='A01.1.00.001',kind='exact',systems=['skeletal'],names=['head'])]
  r=self.resolve('head','TA98FIBULA:A02.5.07.002','skeletal',dict(ta98Term='A02.5.07.002',ta98Kind='parent'));self.assertEqual(r['_term'],'A02.5.07.002')
 def test_source_scoped_review_can_correct_a_stale_tag(self):
  cid='TA98FIBULA:A02.5.07.002';self.env['reviewed']=[dict(term='A02.5.07.001',kind='parent',systems=['skeletal'],names=['head'],conceptIds=[cid])]
  r=self.resolve('head',cid,'skeletal',dict(ta98Term='A02.5.07.002',ta98Kind='exact'));self.assertEqual((r['_term'],r['_kind'],r['_method']),('A02.5.07.001','parent','reviewed'))
 def test_authored_broader_scope_precedes_a_planned_id(self):
  r=self.resolve('head','TA98FIBULA:A02.5.07.002','skeletal',dict(ta98Term='A02.5.07.001',ta98Kind='parent'));self.assertEqual(r['_term'],'A02.5.07.001')
 def test_author_entity_sex_suffix_is_retained(self):
  r=self.resolve('unqualified','TA98REGION:A09.1.02.001F','reproductive',dict(ta98Term='A09.1.02.001F',ta98Kind='exact'));self.assertEqual((r['_term'],r['_kind']),('A09.1.02.001F','exact'))
 def test_arbitrary_local_id_does_not_create_a_ta98_declaration(self):
  r=self.env['declared_identity']('LOCAL:notes:A02.5.07.002',None,self.env['terms']);self.assertIsNone(r)
 def test_cell_and_attachment_annotations_remain_excluded(self):
  for system in ['attachments','cell-membranes']:
   r=self.resolve('head','TA98FIBULA:A02.5.07.002',system,dict(ta98Term='A02.5.07.002',ta98Kind='exact'));self.assertIsNone(r['ta98'])
 def test_invalid_authored_code_is_reported_instead_of_name_matched(self):
  with self.assertRaisesRegex(ValueError,'Invalid authored TA98 identity'):
   self.resolve('head','TA98FIBULA:A02.5.07.002','skeletal',dict(ta98Term='A99.0.00.001',ta98Kind='exact'))
 def test_untagged_body_head_keeps_existing_resolution(self):
  r=self.resolve('head','body:head','regions');self.assertEqual(r['_term'],'A01.1.00.001')

class InstalledCatalogueTest(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  cls.catalogues={sex:ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json' for sex in ['male','female']}
  if not all(path.exists() for path in cls.catalogues.values()):
   raise unittest.SkipTest('Installed local review catalogues are not present')
  with tempfile.TemporaryDirectory(prefix='ta98-current-catalogues-') as folder:
   output=pathlib.Path(folder)/'metadata.json'
   subprocess.run([sys.executable,str(ROOT/'scripts/build-ta98-metadata.py'),'--output',str(output)],cwd=ROOT,check=True,capture_output=True,text=True)
   cls.metadata=json.loads(output.read_text())
 def test_installed_catalogue_tags_survive_the_complete_builder(self):
  for sex,path in self.catalogues.items():
   catalogue=json.loads(path.read_text())
   active={p['conceptId'] for p in catalogue['parts'] if not p.get('suppressed') and p['system']!='attachments' and not p['system'].startswith('cell-')}
   for concept in catalogue['concepts']:
    cid=concept['id']
    if cid not in active or not concept.get('ta98Term'):continue
    expected=self.metadata['reviewedCatalogueIdentities'].get(cid) or dict(term=concept['ta98Term'],kind=concept.get('ta98Kind') or 'parent')
    with self.subTest(sex=sex,cid=cid):
     actual=self.metadata['byConceptMatch'].get(cid)
     self.assertIsNotNone(actual)
     self.assertEqual({key:actual[key] for key in expected},expected)
 def test_existing_intrarenal_collectives_keep_broader_review_scope(self):
  for namespace in ['ZA','REF']:
   for side in ['left','right']:
    cid=f'{namespace}:Intrarenal veins of {side} kidney'
    self.assertEqual(self.metadata['byConceptMatch'][cid],dict(term='A08.1.04.001',kind='parent',method='reviewed'))
if __name__=='__main__':unittest.main()
