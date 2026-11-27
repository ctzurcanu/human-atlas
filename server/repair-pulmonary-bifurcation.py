"""Union overlapping native branch envelopes; classifier caps never draw."""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,copy,json,sys,time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit,topology
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.communicating_vascular_network import boolean_union
from geometry_tools.algorithms.finite_tissue_boolean import barycentric
ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();started=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';source=h/args.candidate;path=source/'atlas.json';a=json.loads(path.read_text());out=h/args.stage;out.mkdir(exist_ok=False);pid='ZA:Bifurcation of pulmonary trunk';p=next(p for p in a['parts']if p['id']==pid);m=AtlasReader().read_part(a,p,frame='assembled-metres',catalogue_path=path);edges=m.faces[:,[[0,1],[1,2],[2,0]]].reshape(-1,2);graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(m.vertices),len(m.vertices)));count,component=connected_components(graph,directed=False);assert count==2;tubes=[];source_faces={};cap_ids=[];offset=0
for ci in range(count):
 fis=np.flatnonzero(component[m.faces[:,0]]==ci);faces=m.faces[fis];used,inv=np.unique(faces,return_inverse=True);v=m.vertices[used].copy();f=inv.reshape(-1,3);all_edges=[tuple(map(int,e))for t in f for e in zip(t,np.roll(t,-1))];counts=Counter(tuple(sorted(e))for e in all_edges);boundary=[e for e in all_edges if counts[tuple(sorted(e))]==1];by=defaultdict(set)
 for i,j in boundary:by[i].add(j);by[j].add(i)
 remaining=set(by);caps=[];closed=f.tolist()
 while remaining:
  stack=[next(iter(remaining))];ring=[]
  while stack:
   i=stack.pop()
   if i not in remaining:continue
   remaining.remove(i);ring.append(i);stack.extend(by[i])
  assert len(ring)==12;center=len(v);v=np.r_[v,v[ring].mean(0)[None]];ids=[];ring_set=set(ring)
  for i,j in boundary:
   if i in ring_set:ids.append(len(closed));closed.append([j,i,center])
  caps.append(np.array(ids))
 for i,fi in enumerate(fis):source_faces[offset+i]=int(fi)
 cap_ids.extend((np.concatenate(caps)+offset).tolist());tubes.append((str(ci),v,np.array(closed,np.uint32),caps));offset+=len(closed)
v,f,parents,provenance,guide=boolean_union(tubes,m.vertices.mean(0));take=~np.isin(parents,cap_ids);assert take.any();f=f[take];parents=parents[take];assert all(int(i)in source_faces for i in parents);source_ids=np.array([source_faces[int(i)]for i in parents]);tree=cKDTree(m.vertices);dist,ix=tree.query(v);v[dist<1e-12]=m.vertices[ix[dist<1e-12]];origin=np.array(p.get('geometryOrigin',[0,0,0]));pos=(v-origin).astype('<f4');actual=pos.astype(float)+origin;w=np.array([barycentric(actual[t],m.triangles[fi])for t,fi in zip(f,source_ids)]);corners=pos[f].reshape(-1,3);attrs={}
for key,value in m.attributes.items():
 field=np.einsum('fij,fjk->fik',w,value[m.faces[source_ids]][...,None])[...,0]if value.ndim==1 else np.einsum('fij,fjk->fik',w,value[m.faces[source_ids]])
 if key=='normals':field/=np.linalg.norm(field,axis=2)[:,:,None];field=np.rint(field*32767).astype('<i2')
 elif key=='uvs':field=field.astype('<f4')
 elif key=='normalSigns':field=np.rint(field).astype('i1')
 else:field=np.rint(field).astype('u1')
 attrs[key]=field.reshape((-1,field.shape[-1]))if value.ndim>1 else field.reshape(-1)
buffers=[];at=0
def put(array):
 global at
 pad=(-at)%4
 if pad:buffers.append(bytes(pad));at+=pad
 result=at;data=array.tobytes();buffers.append(data);at+=len(data);return result
layout=dict(positions=put(corners.astype('<f4')))
for key,value in attrs.items():layout[key]=put(value)
indices=np.arange(len(corners),dtype='<u4').reshape(-1,3);layout['indices']=put(indices);bank=out/'bifurcation.bin';atomic_bytes(bank,b''.join(buffers));p.update(layout,chunk=len(a['chunks']),vertexCount=len(corners),indexCount=int(indices.size),triangles=len(indices),bounds=np.array([actual[f].min((0,1)),actual[f].max((0,1))]).tolist());a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size));p['pulmonaryBifurcationRepair']=dict(schema='native-branch-envelope-union/v1',sourceCatalogueSha256=digest_file(path),calculationCapsExcluded=True,finiteVesselWallEstablished=False,originalMouthsRetained=True);candidate_mesh=AtlasReader().read_part(a,p,frame=m.frame,catalogue_path=path);quality=audit(candidate_mesh,max_pairs=500000,max_tests=4000000);report=dict(source_catalogue_sha256=digest_file(path),part_id=pid,original=topology(m),emitted=quality,classifier_caps_excluded=True,retained_faces=len(indices),construction_guide_components=int(guide.num_vert()>0),elapsed_seconds=time.monotonic()-started,passed=quality['passed']);atomic_json(out/'repair.json',report);print('Bifurcation quality',quality['passed'],quality['failures'],quality['self_intersections'].get('defect_pairs'),quality['topology']['components'],flush=True);assert quality['passed'];atomic_json(out/'atlas.json',a);packing=json.loads((source/'packing.json').read_text());packing.update(source_candidate_sha256=digest_file(path),bifurcation_repair=str(out/'repair.json'),bifurcation_finite_wall_established=False);atomic_json(out/'packing.json',packing)
