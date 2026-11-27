"""Real temporary banks/manifests test partial and changed stages fail closed."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'server'))
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh
import tongue_packing_stages as stages


class StageTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory(dir=ROOT/'.local-models/geometry-workbench')
        self.base=Path(self.temp.name);self.output=self.base/'output';self.output.mkdir()
        self.source=self.base/'source.bin';self.source.write_bytes(b'actual small input')
        self.code=self.base/'code.py';self.code.write_text('code = 1\n')
        self.plan=self.base/'inputs.json'
        atomic_json(self.plan,dict(bindings=dict(source=dict(path=str(self.source),sha256=digest_file(self.source)))))
        self.runtime=self.base/'runtime.json'
        atomic_json(self.runtime,dict(schema='bounded-command/v1',pid=123,state='completed',returncode=0,worker_returncode=0))
        self.runtime_binding=dict(path=str(self.runtime),worker_pid=123)
        self.codes=patch.object(stages,'code_files',return_value={self.code});self.codes.start()
        self.identity=patch.object(stages,'runtime_identity',return_value=self.runtime_binding);self.identity.start()
    def tearDown(self):
        self.identity.stop();self.codes.stop();self.temp.cleanup()
    def fixture(self):
        snapshot=stages.begin_geometry(self.plan)
        vertices=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype='<f4')
        faces=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]],dtype='<u4')
        bank=self.output/'geometry.bin';bank.write_bytes(vertices.tobytes()+faces.tobytes())
        parts=[];records=[]
        for i in range(12):
            pid='test:'+str(i);origin=[i*2.,0.,0.]
            mesh=Mesh(vertices.astype(float)+origin,faces,'model:local-ta98:assembled-metres')
            parts.append(dict(id=pid,chunk=0,vertexCount=4,indexCount=12,positions=0,indices=vertices.nbytes,geometryOrigin=origin))
            records.append(dict(part_id=pid,decoded_geometry_sha256=mesh.geometry_sha,topology=topology(mesh)))
        marker=dict(activePhysicalMeshes=12,partIds=[p['id']for p in parts],termGroups=[{}for _ in range(17)])
        atlas=dict(parts=parts,chunks=[dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')))],finiteTongueIntegration=marker)
        atomic_json(self.output/'pending-atlas.json',atlas)
        atomic_json(self.output/'pending-inventory.json',dict(term_count=36))
        atomic_json(self.output/'decoded-packing-proof.json',dict(schema='tongue-coupled-pending-geometry/v1',
            all_input_gates_passed=True,installed=False,records=records))
        stages.seal_geometry(self.output,snapshot,lambda:None)
        return stages.load_checkpoint(self.output)
    def verify_all(self):
        checkpoint,_,_=stages.load_checkpoint(self.output)
        for batch in range(checkpoint['verification_batches']):
            stages.verify_stage(self.output,batch)
    def test_discovery_deduplicates_actual_files_without_hashing(self):
        with patch.object(stages,'digest_file',side_effect=AssertionError('discovery must not hash')):
            found=stages.discover_files([dict(witness=str(self.source))for _ in range(1000)])
        self.assertEqual(found,{str(self.source.resolve())})
    def test_real_bank_verification_and_publication(self):
        self.fixture();self.verify_all()
        result=stages.finalize_stage(self.output)
        self.assertTrue(result['candidate_published'])
        self.assertFalse(result['installed'])
        stages.require_finalized_stage(self.output/'atlas.json')
        evidence=stages.read(self.output/'packing-evidence.json')
        self.assertIn(str(self.source),evidence['input_bindings'])
        self.assertIn(str(self.code),evidence['input_bindings'])
        proof=stages.read(self.output/'verified-batch-0.json')['independent_geometry_read']
        self.assertEqual(len(proof['actual_meshes']),12)
        self.assertEqual((self.output/'atlas.json').read_bytes(),(self.output/'pending-atlas.json').read_bytes())
    def test_missing_batch_cannot_publish(self):
        self.fixture()
        with self.assertRaises(FileNotFoundError):stages.finalize_stage(self.output)
        self.assertFalse((self.output/'atlas.json').exists())
    def test_changed_source_after_verification_cannot_publish(self):
        self.fixture();self.verify_all();self.source.write_bytes(b'changed actual input')
        with self.assertRaisesRegex(ValueError,'snapshotted'):stages.finalize_stage(self.output)
        self.assertFalse((self.output/'atlas.json').exists())
    def test_changed_loaded_code_cannot_publish(self):
        self.fixture();self.verify_all();self.code.write_text('code = 2\n')
        with self.assertRaisesRegex(ValueError,'snapshotted'):stages.finalize_stage(self.output)
    def test_geometry_timeout_even_after_checkpoint_is_rejected(self):
        self.fixture();report=stages.read(self.runtime);report.update(state='timed_out',returncode=124)
        atomic_json(self.runtime,report,overwrite=True)
        with self.assertRaisesRegex(ValueError,'did not finish successfully'):stages.verify_stage(self.output,0)
    def test_timed_out_finalizer_never_qualifies(self):
        self.fixture();self.verify_all();stages.finalize_stage(self.output)
        report=stages.read(self.runtime);report.update(state='timed_out',returncode=124);atomic_json(self.runtime,report,overwrite=True)
        with self.assertRaisesRegex(ValueError,'did not finish successfully'):stages.require_finalized_stage(self.output/'atlas.json')
    def test_bad_geometry_is_independently_rejected(self):
        snapshot=stages.begin_geometry(self.plan)
        (self.output/'geometry.bin').write_bytes(b'not a geometry bank')
        atomic_json(self.output/'pending-atlas.json',dict(finiteTongueIntegration={}))
        atomic_json(self.output/'pending-inventory.json',dict(term_count=36))
        atomic_json(self.output/'decoded-packing-proof.json',dict(schema='tongue-coupled-pending-geometry/v1',all_input_gates_passed=True,installed=False))
        stages.seal_geometry(self.output,snapshot,lambda:None)
        with self.assertRaisesRegex(ValueError,'identities changed'):stages.verify_stage(self.output,0)
    def test_manifest_does_not_reuse_old_partial_output(self):
        (self.output/'geometry.bin').write_bytes(b'unqualified partial output')
        with self.assertRaises(FileNotFoundError):stages.load_checkpoint(self.output)


class MaterialTests(unittest.TestCase):
    def test_four_uv_free_parts_only_and_native_material_unchanged(self):
        spec=importlib.util.spec_from_file_location('tongue_packer_test',ROOT/'server/build-tongue-coupled-candidate.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        mesh=Mesh(vertices,np.array([[0,1,2]]),'test')
        materials={'native':dict(color=[120,80,70],map='surface-grain.png')}
        for family in ('Genioglossus','Hyoglossus'):
            for side in ('l','r'):
                part=dict(id='ZA:'+family+' muscle.'+side,material='native')
                module.selected_extrinsic_material(part,mesh,materials)
                self.assertNotEqual(part['material'],'native')
                self.assertFalse(part['materialCorrection']['sampledTextures'])
                self.assertNotIn('map',materials[part['material']])
        self.assertEqual(materials['native'],dict(color=[120,80,70],map='surface-grain.png'))
        with self.assertRaisesRegex(ValueError,'restricted to the four'):
            module.selected_extrinsic_material(dict(id='unrelated',material='native'),mesh,materials)
        mesh.attributes['uvs']=np.zeros((3,2),np.float32)
        with self.assertRaisesRegex(ValueError,'UV-free'):
            module.selected_extrinsic_material(dict(id='ZA:Genioglossus muscle.l',material='native'),mesh,materials)


if __name__=='__main__':unittest.main()
