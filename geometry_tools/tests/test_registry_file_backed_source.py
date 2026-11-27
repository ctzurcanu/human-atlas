import json,tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import AtlasReader
from geometry_tools.registry import file_backed_source


class RegistryFileBackedSourceTests(unittest.TestCase):
 def test_sparse_source_maps_and_attributes_reconstruct_from_bound_bank(self):
  with tempfile.TemporaryDirectory(prefix='atlas-registry-source-')as directory:
   root=Path(directory);bank=bytearray()
   def add(value):
    while len(bank)%4:bank.append(0)
    at=len(bank);bank.extend(value.tobytes());return at
   p=dict(id='fixture',vertexCount=6,indexCount=6,chunk=0,geometryOrigin=[3.1,4.2,5.3],geometryScale=2.,sourceOffset=[.01,.02,.03]);p['positions']=add(np.array([[0,0,0],[99,99,99],[1,0,0],[88,88,88],[0,1,0],[1,1,0]],dtype='<f4'));p['normals']=add(np.tile(np.array([0,0,32767],dtype='<i2'),(6,1)));p['normalSigns']=add(np.ones(6,dtype='i1'));p['uvs']=add(np.arange(12,dtype='<f4').reshape(6,2));p['colors']=add(np.arange(18,dtype='u1').reshape(6,3));p['indices']=add(np.array([[0,2,4],[4,2,5]],dtype='<u4'));path=root/'public/models/bank.bin';path.parent.mkdir(parents=True);path.write_bytes(bank);a=dict(chunks=[dict(url='/models/bank.bin')]);reader=AtlasReader(root);m=reader.read_part(a,p,frame='fixture-metres');source=file_backed_source(m.source);self.assertNotIn('source_vertex_ids',source);self.assertNotIn('source_face_ids',source);self.assertEqual(m.source['source_vertex_ids'],[0,2,4,5]);reloaded=reader.reload(json.loads(json.dumps(source)),m.frame);self.assertEqual(reloaded.geometry_sha,m.geometry_sha);self.assertEqual(reloaded.attributes_sha,m.attributes_sha);self.assertEqual(reloaded.source['source_vertex_ids'],[0,2,4,5]);self.assertEqual(reloaded.source['source_face_ids'],[0,1]);np.testing.assert_array_equal(reloaded.vertices,m.vertices)
   altered=bytearray(bank);altered[0]^=1;path.write_bytes(altered)
   with self.assertRaises(ValueError):reader.reload(source,m.frame)
 def test_non_atlas_provenance_is_retained(self):
  source=dict(kind='npz',path='/example.npz',source_vertex_ids=[4,2],source_face_ids=[7]);self.assertEqual(file_backed_source(source),source)

if __name__=='__main__':unittest.main()
