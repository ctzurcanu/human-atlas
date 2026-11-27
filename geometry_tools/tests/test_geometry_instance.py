import copy,json,tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.tissue import validate_cells

class GeometryInstanceTests(unittest.TestCase):
 def setUp(self):
  self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name);(self.root/'.local-models').mkdir();v=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]],'<f4');n=np.tile([0,0,32767],(5,1)).astype('<i2');a=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]],'<u4');b=np.array([[0,1,2],[0,4,1],[1,4,2],[2,4,0]],'<u4');data=bytearray();offsets=[]
  for array in [v,n,a,b]:data.extend(b'\0'*((-len(data))%4));offsets.append(len(data));data.extend(array.tobytes())
  (self.root/'.local-models/bank.bin').write_bytes(data);q=[0.,0.,float(np.sin(.31)),float(np.cos(.31))];common=dict(positions=offsets[0],normals=offsets[1],vertexCount=5,indexCount=12,chunk=0,geometryOrigin=[.2,1.3,-.1],geometryInstance=dict(rotation=q,scale=.91));self.parts=[dict(common,id='a',indices=offsets[2]),dict(common,id='b',indices=offsets[3])];self.atlas=dict(chunks=[dict(url='/local-models/bank.bin')]);self.reader=AtlasReader(root=self.root)
 def tearDown(self):self.temp.cleanup()
 def read(self,p):return self.reader.read_part(self.atlas,p,frame='model:test')
 def test_positive_similarity_keeps_literal_geometry_and_scales_metrics(self):
  m=self.read(self.parts[0]);result=audit(m,require_closed=True);self.assertTrue(result['passed']);self.assertTrue(result['geometry_instance']['all_world_positions_and_normal_fields_replayed']);self.assertAlmostEqual(result['topology']['oriented_triangle_integral_m3'],.91**3/6)
 def test_mutation_cannot_borrow_instance_proof(self):
  m=self.read(self.parts[0]);bad=Mesh(m.vertices+[1e-5,0,0],m.faces,m.frame,m.attributes,m.source)
  with self.assertRaises(ValueError):audit(bad,require_closed=True)
 def test_invalid_or_reflected_placement_rejected(self):
  for change in [dict(rotation=[0,0,0,0],scale=1),dict(rotation=[0,0,0,1],scale=-1)]:
   p=copy.deepcopy(self.parts[0]);p['geometryInstance']=change
   with self.assertRaises(ValueError):self.read(p)
 def test_occupied_volume_checks_use_the_same_verified_preimage(self):
  result=validate_cells({p['id']:self.read(p)for p in self.parts});self.assertTrue(result['passed'],result['failures']);self.assertEqual(len(result['pair_occupied_volume_checks']),1);self.assertAlmostEqual(result['cell_volumes_m3']['a'],.91**3/6)

if __name__=='__main__':unittest.main()
