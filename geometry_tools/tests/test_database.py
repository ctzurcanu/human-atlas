"""Opt-in PostgreSQL integration; every database change is rolled back."""
from contextlib import contextmanager
import os,json,tempfile,unittest,uuid
from pathlib import Path
from unittest.mock import patch
import numpy as np
from psycopg.types.json import Jsonb
from geometry_tools.db import Database
from geometry_tools.common import identity,atomic_json
from geometry_tools.model import Mesh
from geometry_tools.registry import import_catalogue,load_registered,search
from geometry_tools.constraints import import_forest,closure,add_constraint
from geometry_tools.jobs import execute,row,artifacts
from geometry_tools.spatial import cached_bvh
from geometry_tools.contact import contacts
from geometry_tools.fit import fit
from test_core import tetra


@unittest.skipUnless(os.environ.get('ATLAS_GEOMETRY_TEST_DB')=='1','Enable local PostgreSQL integration with ATLAS_GEOMETRY_TEST_DB=1')
class DatabaseTests(unittest.TestCase):
    def setUp(self):
        real=Database();self.context=real.connect();self.conn=self.context.__enter__();self.conn.execute('SAVEPOINT geometry_test')
        class RollbackDatabase(Database):
            @contextmanager
            def connect(inner):yield self.conn
        self.db=RollbackDatabase();self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name).resolve();self.patchers=[]
        for name in ['geometry_tools.spatial.CACHE','geometry_tools.contact.CACHE','geometry_tools.jobs.CACHE']:
            p=patch(name,self.root/'cache');p.start();self.patchers.append(p)

    def tearDown(self):
        for p in self.patchers:p.stop()
        self.conn.execute('ROLLBACK TO SAVEPOINT geometry_test');self.context.__exit__(None,None,None);self.tmp.cleanup()

    def fixture_atlas(self,sex='female',side='left'):
        project=self.root/sex;(project/'public').mkdir(parents=True,exist_ok=True);m=tetra();raw=bytearray(m.vertices.astype('<f4').tobytes());offset=len(raw);raw.extend(m.faces.astype('<u4').tobytes());(project/'public/bank.bin').write_bytes(raw)
        atlas={'sex':sex,'chunks':[{'url':'/bank.bin'}],'concepts':[{'id':'kidney','ta98Term':'A08.1.01.001','ta98Kind':'exact'}],'parts':[{'id':'kidney','conceptId':'kidney','name':'Kidney ('+side+')','system':'urinary','chunk':0,'vertexCount':4,'indexCount':12,'positions':0,'indices':offset,'license':'fixture'}]}
        path=project/'atlas.json';atomic_json(path,atlas);model='fixture-'+sex+'-'+uuid.uuid4().hex
        result=import_catalogue(self.db,path,model_id=model,root=project,frame_id='fixture');return model,path,atlas,result

    def test_registry_stale_bank_filtering_and_shape_search(self):
        model,p,a,r=self.fixture_atlas();model2,p2,a2,r2=self.fixture_atlas('male','right')
        source=load_registered(self.db,'kidney',model)[0]
        result=search(self.db,similar=source['id'],model=model2,sex='male',side='right',term='A08.1.01.001',license_name='fixture')['results'];self.assertEqual(len(result),1);self.assertAlmostEqual(result[0]['shape_distance'],0,places=6)
        self.assertFalse(search(self.db,similar=source['id'],model=model2,side='left')['results'])
        bank=p.parent/'public/bank.bin';old=bank.read_bytes();bank.write_bytes(old+b'extra');new=import_catalogue(self.db,p,model_id=model,root=p.parent,frame_id='fixture');self.assertNotEqual(new['revision_sha256'],r['revision_sha256'])
        with self.assertRaises(ValueError):load_registered(self.db,source['id'])
        a['parts'][0]['name']='Changed current source';atomic_json(p,a,overwrite=True)
        with self.assertRaises(ValueError):search(self.db,model=model)

    def test_identical_catalogue_in_two_models_keeps_both_searchable(self):
        model,path,atlas,first=self.fixture_atlas()
        other='fixture-shared-'+uuid.uuid4().hex
        second=import_catalogue(self.db,path,model_id=other,root=path.parent,frame_id='fixture')
        self.assertNotEqual(first['revision_sha256'],second['revision_sha256'])
        for target,revision in ((model,first),(other,second)):
            matches=search(self.db,query='Kidney',model=target)['results']
            self.assertEqual(len(matches),1)
            self.assertEqual(load_registered(self.db,'kidney',target)[0]['catalogue_sha256'],revision['revision_sha256'])

    def test_complete_forest_expansion_and_fixed_blocker(self):
        model,p,a,r=self.fixture_atlas();a['parts'].append(dict(a['parts'][0],id='other',name='Other'));atomic_json(p,a,overwrite=True);r=import_catalogue(self.db,p,model_id=model,root=p.parent,frame_id='fixture')
        from geometry_tools.common import digest_file
        forest={'sourceSha256':digest_file(p),'forest':[{'id':'frame','children':[{'id':'assembly','name':'Both','partIds':['kidney','other'],'activePartIds':['kidney','other']}]}],'interfaces':[{'from':'kidney','to':'other','kind':'anatomical:innervates','coMove':False}]};fp=self.root/'forest.json';atomic_json(fp,forest)
        result=import_forest(self.db,fp,model,fixed_systems=());self.assertEqual(result['owned_items'],2);self.assertEqual(result['total_forest_relations_preserved'],1);self.assertEqual(result['imported_interface_declarations'],0)
        add_constraint(self.db,{'kind':'fixed','first':'other','evidence':{'source':'test policy'},'parameters':{}},model)
        c=closure(self.db,['kidney'],model);self.assertEqual(len(c['members']),2);self.assertFalse(c['movable'])
        with self.assertRaises(ValueError):fit(self.db,{'model':model,'moving':['kidney'],'mode':'translation','translation_m':[0,0,.01]},self.root/'fit')

    def test_staged_registry_publishes_only_the_complete_immutable_revision(self):
        model,path,atlas,baseline=self.fixture_atlas()
        atlas['parts'].extend(dict(atlas['parts'][0],id=name,name=name)for name in ['second','third'])
        atomic_json(path,atlas,overwrite=True)
        first=import_catalogue(self.db,path,model_id=model,root=path.parent,frame_id='fixture',max_parts=1)
        self.assertFalse(first['complete']);self.assertEqual(first['parts'],1)
        current=self.conn.execute('SELECT catalogue_sha256 FROM atlas_geom.current_catalogues WHERE model_id=%s',(model,)).fetchone()
        self.assertEqual(current['catalogue_sha256'],baseline['revision_sha256'])
        second=import_catalogue(self.db,path,model_id=model,root=path.parent,frame_id='fixture',max_parts=1)
        self.assertEqual(second['parts'],2);self.assertEqual(second['parts_imported_this_run'],1);self.assertFalse(second['complete'])
        final=import_catalogue(self.db,path,model_id=model,root=path.parent,frame_id='fixture',max_parts=1)
        self.assertTrue(final['complete']);self.assertEqual(final['parts'],3);self.assertEqual(final['parts_imported_this_run'],1)
        self.assertEqual(load_registered(self.db,'third',model)[0]['catalogue_sha256'],final['revision_sha256'])
        cached=import_catalogue(self.db,path,model_id=model,root=path.parent,frame_id='fixture',max_parts=1)
        self.assertTrue(cached['unchanged']);self.assertEqual(cached['parts_imported_this_run'],0)

    def test_cache_and_jobs_reject_changed_artifacts(self):
        a=tetra();b=tetra((2,0,0));ba,first=cached_bvh(self.db,a);bb,second=cached_bvh(self.db,a);self.assertFalse(first['cache_hit']);self.assertTrue(second['cache_hit'])
        aa,r=contacts(self.db,a,b,1.1);ab,s=contacts(self.db,a,b,1.1);self.assertTrue(s['cache_hit']);np.testing.assert_array_equal(aa['pairs'],ab['pairs'])
        p=self.root/'mesh.npz';a.save(p);spec={'kind':'audit','parameters':{'mesh':str(p),'controls':{'require_closed':True,'expected_components':1,'expected_genus':[0]}}}
        first=execute(self.db,spec);again=execute(self.db,spec);self.assertEqual(first['state'],'passed');self.assertTrue(again['cache_hit'])
        out=row(self.db,first['job_id']);files=artifacts(self.db,out['id']);Path(files[0]['path']).write_bytes(b'changed result')
        with self.assertRaises(ValueError):execute(self.db,spec)

    def test_non_fixed_group_fit_runs_all_checks_and_never_installs(self):
        p=self.root/'source.npz';m=tetra();m.save(p);mid=self.db.register_artifact_mesh(m,p,part_id='fixture',name='Fixture')
        result=fit(self.db,{'moving':[mid],'mode':'translation','translation_m':[0,.01,0],'protected':[],'require_closed':True},self.root/'fit')
        self.assertTrue(result['passed']);self.assertFalse(result['atlas_modified']);self.assertEqual(len(result['candidates']),1)
        candidate=Mesh.load(result['candidates'][0]['path']);np.testing.assert_allclose(candidate.vertices,m.vertices+[0,.01,0]);np.testing.assert_array_equal(Mesh.load(p).vertices,m.vertices)

    def test_icp_job_binds_target_moves_complete_closure_and_respects_fixed_items(self):
        model,p,a,_=self.fixture_atlas()
        a['parts'].append(dict(a['parts'][0],id='attached',name='Attached tissue'))
        atomic_json(p,a,overwrite=True);import_catalogue(self.db,p,model_id=model,root=p.parent,frame_id='fixture')
        from geometry_tools.common import digest_file
        forest={'sourceSha256':digest_file(p),'forest':[{'id':'organ','name':'Complete organ','partIds':['kidney','attached'],'activePartIds':['kidney','attached']}],'interfaces':[]}
        fp=self.root/'forest.json';atomic_json(fp,forest);import_forest(self.db,fp,model,fixed_systems=())
        source=tetra();target=Mesh(source.vertices+[.01,-.005,.008],source.faces,'recipient')
        tp=self.root/'target.npz';target.save(tp)
        target_id=self.db.register_artifact_mesh(target,tp,part_id='target',name='Target')
        spec={'kind':'fit','parameters':{'model':model,'moving':['kidney'],'mode':'icp',
              'registration_source':'kidney','registration_target':target_id,'target_frame':'recipient',
              'initial_matrix':np.eye(4).tolist(),'require_closed':True,
              'icp':{'max_correspondence_m':.05,'voxel_sizes_m':[.03,.015],
                     'max_correction_m':.04,'max_rotation_degrees':5,'sample_count':2000,'minimum_fitness':.95}}}
        result=execute(self.db,spec);self.assertTrue(result['passed'])
        self.assertFalse(result['result']['atlas_modified']);self.assertEqual(len(result['result']['candidates']),2)
        job=row(self.db,result['job_id']);self.assertIn(target_id,job['parameters']['mesh_bindings'])
        for item in result['result']['candidates']:
            candidate=Mesh.load(item['path']);self.assertEqual(candidate.frame,'recipient')
            np.testing.assert_allclose(candidate.vertices,target.vertices,atol=5e-4)
        add_constraint(self.db,{'kind':'fixed','first':'attached','evidence':{'source':'fixed skeleton fixture'},'parameters':{}},model)
        with self.assertRaisesRegex(ValueError,'fixed'):
            fit(self.db,spec['parameters'],self.root/'blocked-icp')
        # Cached registration must never be reused after its target changes.
        tp.write_bytes(tp.read_bytes()+b'changed target')
        with self.assertRaises(ValueError):execute(self.db,spec)


if __name__=='__main__':unittest.main()
