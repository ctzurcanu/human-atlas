"""A positive similarity places a literal local geometry bank without rebaking.

Topology, self-intersection and occupancy are invariant under this placement.
Numerical proofs can therefore use the literal preimage, after verifying every
world position and normal against the declared proper pose and source bank.
"""
import copy
import numpy as np

def instance_pose(record):
 value=record.get('geometryInstance')
 if value is None:return None
 q=np.asarray(value.get('rotation'),float);scale=float(value.get('scale',0))
 if q.shape!=(4,)or not np.isfinite(q).all()or not np.isfinite(scale)or scale<=0 or abs(q@q-1)>1e-8:raise ValueError('Only a finite unit quaternion and positive uniform instance scale are allowed')
 q=q/np.linalg.norm(q);x,y,z,w=q;r=np.array([[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],[2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],[2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]])
 return r,scale

def instance_preimage(mesh):
 source=mesh.source;record=source.get('record',{});pose=instance_pose(record)
 if pose is None:return None
 from .model import AtlasReader
 from pathlib import Path
 path=Path(source['bank_path']);root=Path(source['project_root']);url='/local-models/'+str(path.relative_to(root/'.local-models'))if path.is_relative_to(root/'.local-models')else'/'+str(path.relative_to(root/'public'));local=copy.deepcopy(record);local.pop('geometryInstance');local['geometryOrigin']=[0.,0.,0.];local['chunk']=0;reader=AtlasReader(root=root);canonical=reader.read_part({'chunks':[{'url':url}]},local,frame=mesh.frame);r,scale=pose;origin=np.asarray(record.get('geometryOrigin',[0,0,0]),float);expected=(canonical.vertices@r.T)*scale+origin
 if canonical.source['bank_sha256']!=source['bank_sha256']or not np.array_equal(expected,mesh.vertices)or not np.array_equal(canonical.faces,mesh.faces):raise ValueError('Instance does not replay the literal stored geometry and declared placement')
 for key,value in canonical.attributes.items():
  actual=value.astype(float)@r.T if key=='normals'else value
  if key not in mesh.attributes or not np.array_equal(actual,mesh.attributes[key]):raise ValueError('Instance corner fields do not replay their declared normal basis')
 return canonical,dict(rotation=record['geometryInstance']['rotation'],scale=scale,origin=origin.tolist(),proper_positive_similarity=True,stored_geometry_sha256=canonical.geometry_sha,stored_attributes_sha256=canonical.attributes_sha,bank_sha256=source['bank_sha256'],all_world_positions_and_normal_fields_replayed=True)
