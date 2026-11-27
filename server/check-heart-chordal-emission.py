"""Replay actual stored chordal exteriors against the immutable native sources."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.distance import barycentric
ap=argparse.ArgumentParser();ap.add_argument('--bank',required=True);ap.add_argument('--candidate',required=True);args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';bank=h/args.bank;stage=h/args.candidate;meta=json.loads((bank/'construction.json').read_text());assert meta['passed'];data=(bank/'cells.bin').read_bytes();assert digest_file(bank/'cells.bin')==meta['bank_sha256'];a=json.loads((stage/'atlas.json').read_text());reader=AtlasReader();rows=[]
def view(record,key,dtype,width):return np.frombuffer(data,dtype=dtype,count=record['vertexCount']*width,offset=record[key]).reshape(-1,width)
for row in meta['parts']:
 i=row['source'];arr=np.load(bank/f'part-{i}/arrangement.npz');p=arr['source_positions'];q=arr['positions'];f=arr['source_faces'];original=Mesh.load(h/'apparatus-source-v203'/f'part-{i}.npz');r0=original.source['record'];old=reader._bank(Path(original.source['bank_path']));used=np.array(original.source['source_vertex_ids'],int);w0=reader.view(old,r0['normals'],r0['vertexCount'],'<i2',3)[used][original.faces];uv0=original.attributes.get('uvs');uv0=uv0[original.faces]if uv0 is not None else None
 tri=p[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);bad=np.einsum('tki,ti->tk',w0.astype(float),fn)<=0;assert int(bad.sum())==row['source_normal_corners_repaired'];w0=w0.copy();fixed=np.rint(fn/np.linalg.norm(fn,axis=1)[:,None]*32767).astype('<i2');w0[bad]=np.broadcast_to(fixed[:,None],w0.shape)[bad]
 source_area=np.linalg.norm(fn,axis=1)*.5;covered=np.zeros(len(f));exteriors=[];caps=[];max_plane=0.;literal=0;max_new_error=0.;world_parts=[]
 for role,sid in [('leaflet','leaflet'),('chordae','chord')]:
  r=row['parts'][role];v=view(r,'positions','<f4',3);w=view(r,'normals','<i2',3);uv=view(r,'uvs','<f4',2)if'uvs'in r else None;idx=np.frombuffer(data,dtype='<u4',count=r['indexCount'],offset=r['indices']).reshape(-1,3);signs=np.frombuffer(data,dtype='i1',count=r['vertexCount'],offset=r['normalSigns']);ext=r['tissueDrawGroups'][0]['count']//3;parent=arr[sid+'_parent_faces'];assert len(parent)==ext;actual=v[idx[:ext]];ideal=q[arr[sid+'_exterior']];assert np.array_equal(actual,ideal.astype('<f4'));max_new_error=max(max_new_error,float(np.max(np.linalg.norm(actual.astype(float)-ideal,axis=2))));exteriors.append(actual);caps.append(v[idx[ext:]])
  fn_actual=np.cross(v[idx[:,1]].astype(float)-v[idx[:,0]],v[idx[:,2]].astype(float)-v[idx[:,0]]);assert np.all(np.einsum('tki,ti->tk',(w[idx]*signs[idx,None]).astype(float),fn_actual)>0)
  source_tri=p[f[parent]];normal=fn[parent]/np.linalg.norm(fn[parent],axis=1)[:,None];plane=abs(np.einsum('tki,ti->tk',actual.astype(float)-source_tri[:,None,0].reshape(-1,1,3),normal));max_plane=max(max_plane,float(plane.max()));np.add.at(covered,parent,np.linalg.norm(np.cross(actual[:,1].astype(float)-actual[:,0],actual[:,2].astype(float)-actual[:,0]),axis=1)*.5)
  for corner in range(3):
   weights,valid=barycentric(ideal[:,corner],source_tri);assert valid.all()
   if uv is not None:assert np.array_equal(uv[idx[:ext,corner]],np.einsum('tk,tkj->tj',weights,uv0[parent]).astype('<f4'))or np.max(abs(uv[idx[:ext,corner]]-np.einsum('tk,tkj->tj',weights,uv0[parent])))<1e-7
   for k in range(3):
    mask=arr[sid+'_exterior'][:,corner]==f[parent,k];assert np.array_equal(w[idx[:ext,corner]][mask],w0[parent[mask],k]);literal+=int(mask.sum())
  pid=r['id'];candidate_id=pid if a['finiteCardiacIntegration']['sex']=='male'else next(p['id']for p in a['parts']if p.get('cardiacTemplateSourcePartId')==pid)if role=='leaflet'else pid.replace(':male:',':female:');candidate_record=next(p for p in a['parts']if p['id']==candidate_id);mesh=reader.read_part(a,candidate_record,frame='local-ta98'if a['finiteCardiacIntegration']['sex']=='male'else'local-female-ta98');world_parts.append(dict(id=candidate_id,geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,bounds=mesh.bounds.tolist()))
 assert np.array_equal(caps[0],caps[1][:,[0,2,1]]);points=np.unique(np.concatenate(exteriors).reshape(-1,3),axis=0);point_set=set(map(tuple,points));assert all(tuple(x)in point_set for x in p.astype('<f4'));assert max_plane<=1.2e-7 and max_new_error<=1.2e-7
 area_error=abs(covered-source_area);perimeter=np.linalg.norm(tri[:,1]-tri[:,0],axis=1)+np.linalg.norm(tri[:,2]-tri[:,1],axis=1)+np.linalg.norm(tri[:,0]-tri[:,2],axis=1);assert np.all(area_error<=perimeter*max_new_error*3+1e-16)
 rows.append(dict(source=i,passed=True,original_controls_literal=True,unchanged_source_normal_corners_checked=literal,source_normal_corners_repaired=int(bad.sum()),reciprocal_caps_literal=True,maximum_new_control_rounding_m=max_new_error,maximum_source_plane_residual_m=max_plane,maximum_parent_face_area_error_m2=float(area_error.max()),source_face_coverage=True,parts=world_parts))
result=dict(passed=True,candidate_catalogue_sha256=digest_file(stage/'atlas.json'),bank_sha256=meta['bank_sha256'],rows=rows,elapsed_seconds=time.monotonic()-start);atomic_json(stage/'emission-replay.json',result);print(json.dumps(result),flush=True)
