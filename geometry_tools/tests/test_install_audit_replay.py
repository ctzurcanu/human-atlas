import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.tissue import validate_cells
from geometry_tools.install import replay_tissue_cell_audit


def tetra():
    return Mesh(np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],float),np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]]),'test:metres')


class InstallAuditReplay(unittest.TestCase):
    def test_positive_instance_replays_literal_source_and_rejects_unrecorded_edits(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);models=root/'.local-models';models.mkdir();t=tetra();bank=models/'tetra.bin';bank.write_bytes(t.vertices.astype('<f4').tobytes()+t.faces.astype('<u4').tobytes());reader=AtlasReader(root=root);atlas={'chunks':[{'url':'/local-models/tetra.bin'}]};record=dict(id='cell',chunk=0,positions=0,indices=48,vertexCount=4,indexCount=12)
            source=reader.read_part(atlas,record,frame=t.frame);checkpoints=root/'checks';self.assertTrue(validate_cells({'cell':source},max_pairs=1000,checkpoint_directory=checkpoints)['passed'])
            placed=reader.read_part(atlas,dict(record,geometryOrigin=[2.,3.,4.],geometryInstance={'rotation':[0.,0.,np.sin(.2),np.cos(.2)],'scale':.8}),frame='recipient:metres')
            self.assertIsNone(replay_tissue_cell_audit(placed,checkpoints,max_pairs=1000))
            replay=replay_tissue_cell_audit(placed,checkpoints,max_pairs=1000,canonical_frame=t.frame);self.assertTrue(replay['passed']);self.assertTrue(replay['geometry_instance']['all_world_positions_and_normal_fields_replayed']);self.assertAlmostEqual(replay['topology']['oriented_triangle_integral_m3'],1/6*.8**3)
            with self.assertRaisesRegex(ValueError,'replay'):
                replay_tissue_cell_audit(Mesh(placed.vertices+.001,placed.faces,placed.frame,source=placed.source),checkpoints,max_pairs=1000,canonical_frame=t.frame)
            with self.assertRaisesRegex(ValueError,'attribute set'):
                replay_tissue_cell_audit(Mesh(placed.vertices,placed.faces,placed.frame,{'colors':np.ones((4,3),dtype='u1')},placed.source),checkpoints,max_pairs=1000,canonical_frame=t.frame)

    def test_literal_geometry_frame_and_attributes_are_required(self):
        with tempfile.TemporaryDirectory() as directory:
            source=tetra();validate_cells({'cell':source},max_pairs=1000,checkpoint_directory=directory)
            self.assertTrue(replay_tissue_cell_audit(source,directory,max_pairs=1000)['passed'])
            moved=Mesh(source.vertices+.1,source.faces,source.frame)
            self.assertIsNone(replay_tissue_cell_audit(moved,directory,max_pairs=1000))
            changed_attributes=Mesh(source.vertices,source.faces,source.frame,{'normals':np.ones((4,3),dtype=np.int16)})
            self.assertIsNone(replay_tissue_cell_audit(changed_attributes,directory,max_pairs=1000))
            self.assertIsNone(replay_tissue_cell_audit(source,directory,max_pairs=999))

    def test_edited_cached_pass_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source=tetra();validate_cells({'cell':source},max_pairs=1000,checkpoint_directory=directory)
            path=next(Path(directory).glob('*.json'));record=json.loads(path.read_text());record['result']['passed']=False;path.write_text(json.dumps(record))
            with self.assertRaisesRegex(ValueError,'edited'):
                replay_tissue_cell_audit(source,directory,max_pairs=1000)

    def test_replay_requires_the_actual_triangle_test_budget(self):
        with tempfile.TemporaryDirectory() as directory:
            source=tetra()
            self.assertTrue(validate_cells({'cell':source},max_pairs=1000,max_tests=100,checkpoint_directory=directory)['passed'])
            self.assertIsNone(replay_tissue_cell_audit(source,directory,max_pairs=1000))
            self.assertTrue(replay_tissue_cell_audit(source,directory,max_pairs=1000,max_tests=100)['passed'])

    def test_discarded_buffer_entries_need_literal_audited_corner_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            base=tetra();vertices=np.vstack([base.vertices,[[7,8,9]]]);uv=np.arange(10,dtype='<f4').reshape(5,2)
            original=Mesh(vertices,base.faces,base.frame,{'uvs':uv})
            validate_cells({'cell':original},max_pairs=1000,checkpoint_directory=directory)
            bank=vertices.astype('<f4').tobytes()+uv.tobytes()+base.faces.astype('<u4').tobytes()
            row=dict(vertexCount=5,indexCount=12,positions=0,uvs=60,indices=100)
            layout=dict(parts={'cell':row},frame=base.frame,origin=[0,0,0],_bank=bank,
                _path='audited-layout.json',_sha256='a'*64,bank='audited.bin',bank_sha256='b'*64)
            compact=Mesh(base.vertices,base.faces,base.frame,{'uvs':uv[:4]},source={'record':{'id':'cell'}})
            self.assertIsNone(replay_tissue_cell_audit(compact,directory,max_pairs=1000))
            replay=replay_tissue_cell_audit(compact,directory,max_pairs=1000,audited_layouts=[layout])
            self.assertTrue(replay['passed']);self.assertEqual(replay['attributes_sha256'],compact.attributes_sha)
            self.assertTrue(replay['buffer_equivalence']['all_emitted_triangles_and_corner_fields_literal'])
            changed_uv=uv[:4].copy();changed_uv[0,0]+=.1
            changed=Mesh(base.vertices,base.faces,base.frame,{'uvs':changed_uv},source=compact.source)
            with self.assertRaisesRegex(ValueError,'corner fields'):
                replay_tissue_cell_audit(changed,directory,max_pairs=1000,audited_layouts=[layout])


if __name__=='__main__':unittest.main()
