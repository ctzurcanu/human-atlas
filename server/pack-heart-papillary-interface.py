"""Emit shared finite native leaflet/chordal/papillary cells with source charts."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_npz,atomic_bytes,digest_file
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.tissue import decode_packed,validate_cells
ap=argparse.ArgumentParser();ap.add_argument('--construction',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';stage=h/args.construction;meta=json.loads((stage/'construction.json').read_text());assert meta['passed'];z=np.load(stage/'arrangement.npz');out=h/args.stage;out.mkdir(exist_ok=False);a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());reader=AtlasReader();m={}
for name,source in meta['source_meshes'].items():
 mesh=reader.read_part(a,source['record'],frame='model:local-ta98:assembled-metres');assert mesh.geometry_sha==source['geometry_sha256']and mesh.attributes_sha==source['attributes_sha256'];m[name]=mesh
origin=np.array(meta['origin_m']);p=z['source_positions'];pm=z['papillary_source_positions'];q=z['positions'];lp=meta['source_meshes']['leaflet']['record'];cp=meta['source_meshes']['chordae']['record'];lc=lp['tissueDrawGroups'][0]['count']//3;cc=cp['tissueDrawGroups'][0]['count']//3
def fields(name,sl):
 mesh=m[name];f=mesh.faces[sl];normals=mesh.attributes['normals'][f].copy();sign=mesh.attributes.get('normalSigns');
 if sign is not None:normals*=sign[f,None]
 return normals,mesh.attributes['uvs'][f].copy()if'uvs'in mesh.attributes else np.zeros((len(f),3,2),dtype='<f4')
ln,lu=fields('leaflet',slice(0,lc));cn,cu=fields('chordae',slice(0,cc));nn,nu=fields('chordae',slice(cc,None));pn,pu=fields('papillary',slice(None));source_defs={'leaflet-exterior':(p,z['source_leaflet_exterior'],ln,lu),'chord-exterior':(p,z['source_chord_exterior'],cn,cu),'neck':(p,z['source_neck'],nn,nu),'occupied-outer':(pm,z['papillary_source_faces'],pn,pu)};repairs={}
for sid,(vertices,faces,words,uv)in source_defs.items():
 tri=vertices[faces];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);bad=np.einsum('tki,ti->tk',words.astype(float),fn)<=0;repairs[sid]=int(bad.sum())
 if bad.any():flat=np.rint(fn/np.linalg.norm(fn,axis=1)[:,None]*32767).astype('<i2');words[bad]=np.broadcast_to(flat[:,None],words.shape)[bad]
surfaces={};words={};uvs={};rows={r['id']:r for r in meta['surface_rows']};proof=[]
for definition in meta['surfaces']:
 sid=definition['id'];fs=z[definition['array']];
 if not len(fs):continue
 row=rows[sid];sp,sf,sn,su=source_defs[row.get('parentSurface',sid)];parent=z[definition['parent_faces']];bary=z[definition['barycentrics']];source_tri=sp[sf[parent]];field=np.einsum('tck,tkj->tcj',bary,sn[parent].astype(float));field/=np.linalg.norm(field,axis=2)[:,:,None];field=np.rint(field*32767).astype('<i2');texture=np.einsum('tck,tkj->tcj',bary,su[parent]).astype('<f4');literal=0
 for k in range(3):
  for j in range(3):
   take=np.all(q[fs[:,k]]==source_tri[:,j],axis=1);field[take,k]=sn[parent[take],j];texture[take,k]=su[parent[take],j];literal+=int(take.sum())
 field*=row.get('parentNormalSign',1);surfaces[sid]=fs;words[sid]=field;uvs[sid]=texture;source_normal=np.cross(source_tri[:,1]-source_tri[:,0],source_tri[:,2]-source_tri[:,0])*row.get('parentNormalSign',1);actual=q.astype('<f4')[fs].astype(float);fn=np.cross(actual[:,1]-actual[:,0],actual[:,2]-actual[:,0]);dot=np.einsum('ti,ti->t',fn,source_normal);bad=np.flatnonzero(dot<=0);proof.append(dict(surface=sid,literal_source_corners=literal,reversed_float32_fragments=int(len(bad)),minimum_parent_orientation_dot=float(dot.min()),source_parent_faces=parent.tolist()if len(bad)else[]))
atomic_json(out/'parent-orientation.json',dict(passed=all(r['reversed_float32_fragments']==0 for r in proof),surfaces=proof));assert all(r['reversed_float32_fragments']==0 for r in proof),'Reversed Float32 fragments; no bank emitted'
packed=pack_finite_cells(q,surfaces,words,meta['cells'],meta['part_ids'],corner_uvs=uvs);atomic_bytes(out/'cells.bin',packed.data);decoded=decode_packed(packed,'heart-native-local-metres',np.zeros(3));checks=validate_cells(decoded,max_pairs=1000000,max_tests=10000000,checkpoint_directory=out/'checkpoints');atomic_json(out/'checked-geometry.json',checks);assert checks['passed'];atomic_npz(out/'source-arrangement.npz',positions=q,**{f'surface_{i}':faces for i,faces in enumerate(surfaces.values())});emission=dict(passed=True,source=meta['source'],construction_stage=args.construction,construction_sha256=digest_file(stage/'construction.json'),source_catalogue_sha256=meta['source_catalogue_sha256'],source_meshes=meta['source_meshes'],parts=packed.cells,bank_sha256=digest_file(out/'cells.bin'),bytes=len(packed.data),geometry_origin=origin.tolist(),source_normal_corners_repaired=repairs,maximum_control_rounding_m=float(np.max(np.linalg.norm(q-q.astype('<f4').astype(float),axis=1))),elapsed_seconds=time.monotonic()-start,installed=False);atomic_json(out/'emission.json',emission);print(json.dumps({k:v for k,v in emission.items()if k not in ['parts','source_meshes']}),flush=True)
