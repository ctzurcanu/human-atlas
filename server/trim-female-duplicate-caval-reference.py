"""Retain the cardiac cuff above the already present native caval segment.

This clips an outer reference only, with an open distal mouth and no cap.
It does not manufacture a finite vessel wall or certify an unbuilt lumen join.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface
ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic()
h=ROOT/'.local-models/geometry-workbench/heart';source=h/args.candidate;out=h/args.stage;out.mkdir(exist_ok=False);path=source/'atlas.json';a=json.loads(path.read_text());reader=AtlasReader();parts={p['id']:p for p in a['parts']};pid='LOCAL:female:heart-template:za-inferior-vena-cava-thoracic-part';p=parts[pid];m=reader.read_part(a,p,frame='model:local-female-ta98:assembled-metres',catalogue_path=path);cut=1.2605;clip=clip_scalar_surface(m.vertices,m.faces,m.vertices[:,1],cut);v=clip['positions'];f=clip['positive_faces'];edges=clip['new_edge_parents'];mu=clip['fractions'];attributes={}
for k,value in m.attributes.items():
 extra=value[edges[:,0]]*(1-mu[:,None])+value[edges[:,1]]*mu[:,None]if value.ndim==2 else value[edges[:,0]]*(1-mu)+value[edges[:,1]]*mu
 if k=='normalSigns':assert np.array_equal(value[edges[:,0]],value[edges[:,1]]);extra=value[edges[:,0]]
 attributes[k]=np.r_[value,extra]
used,inv=np.unique(f,return_inverse=True);v=v[used];f=inv.reshape(-1,3);attributes={k:x[used]for k,x in attributes.items()};origin=np.array([0.,1.25,0.]);buffers=[];offset=0;layout={}
def add(k,array):
 global offset
 pad=(-offset)%4
 if pad:buffers.append(bytes(pad));offset+=pad
 layout[k]=offset;b=array.tobytes();buffers.append(b);offset+=len(b)
add('positions',(v-origin).astype('<f4'))
for k,arr in attributes.items():
 if k=='normals':
  arr=arr.astype(float);arr/=np.linalg.norm(arr,axis=1)[:,None];arr=np.rint(arr*32767).astype('<i2')
 elif k=='uvs':arr=arr.astype('<f4')
 elif k=='normalSigns':arr=arr.astype('i1')
 elif k=='colors':arr=np.rint(arr).astype('u1')
 else:raise ValueError(k)
 add(k,arr)
add('indices',f.astype('<u4'));bank=out/'caval-reference.bin';atomic_bytes(bank,b''.join(buffers));chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
for key in ['positions','normals','normalSigns','uvs','colors','indices','geometryInstance','geometryScale','sourceOffset']:p.pop(key,None)
p.update(layout,chunk=chunk,geometryOrigin=origin.tolist(),vertexCount=len(v),indexCount=int(f.size),triangles=len(f));p['bounds']=reader.read_part(a,p,frame=m.frame,catalogue_path=path).bounds.tolist();p['cardiacRootReferenceTrim']=dict(schema='retained-cardiac-caval-cuff/v1',sourceCatalogueSha256=digest_file(path),nativeCavalPart='TA98HEPATIC:female:region:A12.3.09.001',retainedAboveYMetres=cut,capsAdded=0,finiteWallEstablished=False,exactDistalJoinEstablished=False)
fit=json.loads((source/'packing.json').read_text());fit['additional_modified_parts']=[pid];fit['source_pose_candidate_sha256']=digest_file(path);fit['reference_trim']=p['cardiacRootReferenceTrim'];fit['banks_modified']=False;fit['replacement_banks']=[str(bank)];atomic_json(out/'packing.json',fit);atomic_json(out/'atlas.json',a);atomic_json(out/'trim.json',dict(source_catalogue_sha256=digest_file(path),part_id=pid,triangles_before=len(m.faces),triangles_after=len(f),cut_y_m=cut,caps_added=0,native_caval_geometry_unchanged=True,upper_mouth_preserved=True,elapsed_seconds=time.monotonic()-start,installed=False));print('Retained the upper cardiac cuff; native caval geometry unchanged',flush=True)
