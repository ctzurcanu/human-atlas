"""Bounded source-surface cleanup; no vascular lumen disks are constructed.

Remove only exact zero-area faces and internal branch partitions incident to
nonmanifold edges. Lock every coordinate shared with another active arterial
source. A complete global crossing guard accepts only improving local moves.
"""
from pathlib import Path
from collections import Counter
import argparse, json, sys, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_json,atomic_npz,atomic_bytes,digest_file
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.vertex_clearance_lp import repair_vertex_clearance

def repair_local_diagonals(value,faces,remaining):
 """Retriangulate an existing quad only, preserving its four-edge boundary."""
 guard=LocalCrossingGuard(value,faces,remaining,maximum_displacement_m=.02,maximum_pair_tests=4000000);moves=[];start=time.monotonic()
 for step in range(12):
  if not guard.remaining or time.monotonic()-start>5.:break
  by={}
  for fi,t in enumerate(guard.faces):
   for i,j in zip(t,np.roll(t,-1)):by.setdefault(tuple(sorted((int(i),int(j)))),[]).append(fi)
  near=set(np.asarray(sorted(guard.remaining)).ravel());accepted=False
  for (i,j),fs in by.items():
   if len(fs)!=2 or not near.intersection(fs):continue
   old=guard.faces[fs];x=next(int(k)for k in old[0]if k not in [i,j]);y=next(int(k)for k in old[1]if k not in [i,j])
   if x==y or tuple(sorted((x,y)))in by:continue
   new=np.array([[x,y,i],[y,x,j]])
   old_tri=value[old];old_n=np.cross(old_tri[:,1]-old_tri[:,0],old_tri[:,2]-old_tri[:,0]);new_tri=value[new];new_n=np.cross(new_tri[:,1]-new_tri[:,0],new_tri[:,2]-new_tri[:,0])
   if np.dot(new_n.sum(0),old_n.sum(0))<0:new=new[:,[0,2,1]];new_n=-new_n
   if np.any(np.einsum('ij,ij->i',old_n,new_n)<=0):continue
   if np.linalg.norm(value[new]-value[old],axis=2).max()>.02:continue
   token=guard.check_faces(fs,new,np.ones(2,bool))
   if token is None or len(token.remaining)>=len(guard.remaining):continue
   guard.commit(token);moves.append(dict(faces=fs,remaining=len(guard.remaining)));accepted=True;break
  if not accepted:break
 return guard.faces.copy(),dict(passed=not guard.remaining,moves=moves,guard=guard.evidence(),elapsed_seconds=time.monotonic()-start)

ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--maximum-displacement-mm',type=float,default=1.);ap.add_argument('--right-branch-checkpoint');args=ap.parse_args();assert 0<args.maximum_displacement_mm<=2.;maximum=args.maximum_displacement_mm/1000
h=ROOT/'.local-models/geometry-workbench/heart';source=h/args.candidate;path=source/'atlas.json';a=json.loads(path.read_text());out=h/args.stage;out.mkdir(exist_ok=False);started=time.monotonic()
parts={p['id']:p for p in a['parts']};reader=AtlasReader();hidden=set(a.get('defaultHidden',[]));frame='assembled-metres'
ids=['LOCAL:female:l_intercostal_arteries:1743','LOCAL:female:r_intercostal_arteries:1736','LOCAL:female:l_transverse_cervical_artery:1731','LOCAL:female:aortic_arch:1513']
meshes={p['id']:reader.read_part(a,p,frame=frame,catalogue_path=path)for p in a['parts']if p['system']=='arterial'and p['id']not in hidden and not p.get('suppressed')and not p.get('defaultHidden')}
assert all(pid in meshes for pid in ids)
owners={}
for pid,m in meshes.items():
 for x in np.unique(m.vertices,axis=0):owners.setdefault(tuple(x),set()).add(pid)
report=dict(source_catalogue_sha256=digest_file(path),parts=[],complete=False,passed=False,all_shared_source_controls_locked=True)
buffers=[];offset=0;chunk=len(a['chunks'])
def put(value):
 global offset
 pad=(-offset)%4
 if pad:buffers.append(bytes(pad));offset+=pad
 start=offset;data=value.tobytes();buffers.append(data);offset+=len(data);return start
for pid in ids:
 p=parts[pid];m=meshes[pid];v,iv=np.unique(m.vertices,axis=0,return_inverse=True);all_faces=iv[m.faces]
 area=np.linalg.norm(np.cross(v[all_faces[:,1]]-v[all_faces[:,0]],v[all_faces[:,2]]-v[all_faces[:,0]]),axis=1);live=area>1e-20
 counts=Counter(tuple(sorted(e))for t in all_faces[live]for e in zip(t,np.roll(t,-1)))
 internal=np.array([sum(counts[tuple(sorted(e))]>2 for e in zip(t,np.roll(t,-1)))>=2 for t in all_faces]);take=live&~internal;faces=all_faces[take]
 clean=Mesh(v,faces,frame);c,s=pair_contacts(clean,clean,0.,selfcheck=True,max_pairs=500000,max_tests=4000000);assert s['complete_within_pair_budget']
 bad=c['pairs'][np.isin(c['relation'],['crossing','coplanar-overlap'])];locked=[i for i,x in enumerate(v)if len(owners[tuple(x)])>1]
 tri=v[faces];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];origin=np.array(p['geometryOrigin']);guard=LocalCrossingGuard(v,faces,bad,maximum_displacement_m=maximum,maximum_pair_tests=4000000)
 if pid=='LOCAL:female:r_intercostal_arteries:1736'and args.right_branch_checkpoint:
  checkpoint=Path(args.right_branch_checkpoint);z=np.load(checkpoint);assert np.array_equal(z['source'],v)and np.array_equal(z['faces'],faces);assert np.linalg.norm(z['positions']-v,axis=1).max()<=maximum;assert np.array_equal(z['positions'][locked],v[locked]);result=dict(positions=z['positions'],evidence=dict(passed=True,checkpoint=str(checkpoint),sha256=digest_file(checkpoint),prior_bounded_stages=True))
 else:result=repair_vertex_clearance(v,v,faces,np.repeat(normal[:,None],3,axis=1),origin=origin,global_guard=guard,maximum_displacement_m=maximum,locked=locked,wall_seconds=16.,maximum_lp_solves=300)
 value=result['positions'];diagonals=None
 if guard.remaining and result['evidence'].get('budget_or_stop_reason')=='no improving bounded linear clearance proposal':faces,diagonals=repair_local_diagonals(value,faces,sorted(guard.remaining))
 quality=audit(Mesh(value,faces,frame),max_pairs=500000,max_tests=4000000)
 row=dict(id=pid,removed_zero_area_faces=np.flatnonzero(~live).tolist(),removed_internal_branch_partition_faces=np.flatnonzero(internal&live).tolist(),locked_shared_controls=len(locked),repair=result['evidence'],diagonals=diagonals,audit=quality,maximum_displacement_m=float(np.linalg.norm(value-v,axis=1).max()))
 report['parts'].append(row);atomic_json(out/'repair.json',report,overwrite=True);print(pid,'repair',result['evidence']['passed'],'quality',quality['passed'],quality['failures'],flush=True)
 atomic_npz(out/('solver-'+str(len(report['parts']))+'.npz'),source=v,positions=value,faces=faces,source_vertex_map=iv,source_face_ids=np.flatnonzero(take))
 assert (result['evidence']['passed']or diagonals and diagonals['passed'])and quality['passed'],'Incomplete source repair cannot be emitted'
 assert not diagonals or not diagonals['moves'],'Changed quad diagonals need a distinct source-chart interpolation step before emission'
 assert np.array_equal(value[locked],v[locked])
 assert not p.get('tissueDrawGroups'),'Partitioned attributes need an explicit remap'
 source_faces=np.flatnonzero(take);corner_ids=m.faces[take].ravel();positions=(value[iv[corner_ids]]-origin).astype('<f4');layout=dict(positions=put(positions));attrs={}
 # Preserve each source chart. Transport normal charts with the proper local
 # least-squares rotation of the retained incident edge star.
 affected=set(np.flatnonzero(np.any(value!=v,axis=1)));neighbors=[set()for _ in v]
 for t in faces:
  for i in t:neighbors[i].update(int(j)for j in t if j!=i)
 rotations={}
 for i,near in enumerate(neighbors):
  if i not in affected and not near.intersection(affected):continue
  ns=np.array(sorted(near),int)
  if len(ns)<2:continue
  u,_,vt=np.linalg.svd((v[ns]-v[i]).T@(value[ns]-value[i]));r=vt.T@u.T
  if np.linalg.det(r)<0:vt[-1]*=-1;r=vt.T@u.T
  rotations[i]=r
 for key,field in m.attributes.items():
  data=field[corner_ids].copy()
  if key=='normals':
   for k,i in enumerate(iv[corner_ids]):
    if i in rotations:data[k]=np.rint(rotations[i]@data[k])
  dtype={'normals':'<i2','normalSigns':'i1','uvs':'<f4','colors':'u1'}[key];attrs[key]=data.astype(dtype);layout[key]=put(attrs[key])
 indices=np.arange(len(positions),dtype='<u4').reshape(-1,3);layout['indices']=put(indices);p.update(layout,chunk=chunk,vertexCount=len(positions),indexCount=int(indices.size),triangles=len(indices),bounds=np.array([value[faces].min((0,1)),value[faces].max((0,1))]).tolist());p['sourceSurfaceCleanup']=dict(schema='locked-interface-bounded-source-fold-repair/v1',sourceCatalogueSha256=digest_file(path),maximumDisplacementMetres=row['maximum_displacement_m'],sharedCoordinatesLiteral=True,uvChartsLiteral=True,noLumenCapsAdded=True)
 atomic_npz(out/('source-chart-'+str(len(report['parts']))+'.npz'),source_face_ids=source_faces,source_unique_vertices=v,repaired_unique_vertices=value,source_corner_ids=corner_ids)
bank=out/'source-fold-repairs.bin';atomic_bytes(bank,b''.join(buffers));a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
# Audit literal re-decoded storage and normal fields, rather than the solver's
# temporary arrays alone.
for row in report['parts']:
 m=AtlasReader().read_part(a,parts[row['id']],frame=frame,catalogue_path=path);row['emitted_audit']=audit(m,max_pairs=500000,max_tests=4000000);assert row['emitted_audit']['passed']
report.update(complete=True,passed=True,elapsed_seconds=time.monotonic()-started);atomic_json(out/'repair.json',report,overwrite=True);atomic_json(out/'atlas.json',a)
packing=json.loads((source/'packing.json').read_text());packing.update(source_candidate_sha256=digest_file(path),source_surface_repair=str(out/'repair.json'));atomic_json(out/'packing.json',packing)
