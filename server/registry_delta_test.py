"""A geometry cache must reject changed frames, banks and every decoded field."""
import copy,sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.registry_delta import reusable_geometry

class RegistryDeltaTest(unittest.TestCase):
 def setUp(self):
  self.record=dict(vertexCount=3,indexCount=3,positions=0,indices=36,normals=48,
                   normalSigns=66,uvs=72,colors=96,geometryOrigin=[0,1.25,0],sourceOffset=[0,0,0])
  self.old=dict(record=self.record,error=None,frame_id='male',bank_path='local.bin',bank_sha256='a'*64)
 def reusable(self,old=None,record=None,**changes):
  args=dict(frame='male',bank_path='local.bin',bank_sha256='a'*64);args.update(changes)
  return reusable_geometry(self.old if old is None else old,self.record if record is None else record,**args)
 def test_metadata_changes_reuse_only_the_literal_geometry(self):
  self.assertTrue(self.reusable(record={**self.record,'name':'Renamed','suppressed':True,'defaultHidden':True}))
 def test_every_decoded_buffer_and_transform_is_part_of_the_guard(self):
  for key in self.record:
   with self.subTest(field=key):
    record=copy.deepcopy(self.record);record[key]=[1,2,3]if isinstance(record[key],list)else record[key]+1
    self.assertFalse(self.reusable(record=record))
 def test_bank_hash_path_and_frame_cannot_be_inherited(self):
  for changes in [dict(frame='female'),dict(bank_path='other.bin'),dict(bank_sha256='b'*64)]:
   self.assertFalse(self.reusable(**changes))
  self.assertFalse(self.reusable(old={**self.old,'error':'failed source'}))
  self.assertFalse(reusable_geometry(None,self.record,'male','local.bin','a'*64))
 def test_instance_pose_and_uniform_scale_cannot_reuse_world_geometry(self):
  for field,value in [('geometryScale',1.1),('geometryInstance',dict(rotation=[0,0,0,1],scale=.9))]:
   self.assertFalse(self.reusable(record={**self.record,field:value}))
  instanced={**self.record,'geometryInstance':dict(rotation=[0,0,0,1],scale=.9)}
  old={**self.old,'record':instanced}
  self.assertTrue(self.reusable(old=old,record=copy.deepcopy(instanced)))
  for rotation,scale in [([0,0,0,1],.8),([0,0,1,0],.9)]:
   self.assertFalse(self.reusable(old=old,record={**instanced,'geometryInstance':dict(rotation=rotation,scale=scale)}))

if __name__=='__main__':unittest.main()
