"""Repair shared pulmonary charts on one common pool, pinning the exterior.

This changes no tissue ownership, faces, UVs, airway course or outer contour.
Every bounded proposal must reduce the complete crossing list without adding
a crossing anywhere in the common surface complex. Installation is separate.
"""
from pathlib import Path
import argparse, copy, json, sys, time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.common import atomic_bytes, atomic_json, digest_file
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.vertex_clearance_lp import repair_vertex_clearance
ap=argparse.ArgumentParser();ap.add_argument('--catalogue',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();started=time.monotonic()
path=ROOT/args.catalogue;a=json.loads(path.read_text());parts={p['id']:p for p in a['parts']};out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False);reader=AtlasReader();ids=a['pulmonarySegmentIntegration']['partIds'];frame='model:local-male-ta98:assembled-metres';meshes={i:reader.read_part(a,parts[i],frame=frame,catalogue_path=path)for i in ids}
offset=np.r_[0,np.cumsum([len(meshes[i].vertices)for i in ids])];v,inv=np.unique(np.concatenate([meshes[i].vertices for i in ids]),axis=0,return_inverse=True);maps={i:inv[offset[j]:offset[j+1]]for j,i in enumerate(ids)};cell_faces={i:maps[i][meshes[i].faces]for i in ids};all_faces=np.concatenate(list(cell_faces.values()));_,which=np.unique(np.sort(all_faces,axis=1),axis=0,return_index=True);faces=all_faces[np.sort(which)];exterior=[]
for pid in ids:
 for g in parts[pid]['tissueDrawGroups']:
  if not g.get('neighbor'):exterior.append(cell_faces[pid][g['start']//3:(g['start']+g['count'])//3])
outer=np.unique(np.concatenate(exterior));locked=np.union1d(outer,np.flatnonzero(v[:,1]>=1.475));m=Mesh(v,faces,frame);c,s=pair_contacts(m,m,0.,selfcheck=True,max_pairs=1500000,max_tests=16000000);bad=c['pairs'][np.isin(c['relation'],['crossing','coplanar-overlap'])];atomic_json(out/'initial.json',dict(catalogue_sha256=digest_file(path),contacts=s,bad_faces=bad.tolist(),fixed_outer_controls=len(outer),unique_faces=len(faces)));assert s['complete_within_pair_budget'];print('Initial common-complex crossing pairs',len(bad),flush=True)
tri=v[faces];n=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];words=np.repeat(n[:,None],3,axis=1);origin=np.array([0.,1.3521,0.]);guard=LocalCrossingGuard(v,faces,bad,maximum_displacement_m=.0001,maximum_pair_tests=16000000)
result=repair_vertex_clearance(v,v,faces,words,origin=origin,global_guard=guard,maximum_displacement_m=.0001,locked=locked,wall_seconds=18.,maximum_lp_solves=150);e=result['evidence'];atomic_json(out/'single-vertex-repair.json',e);print('Single-vertex repair leaves',len(guard.remaining),'pairs',flush=True)
# A boundary fan can need two interior corners to move simultaneously: moving
# either alone temporarily introduces a crossing. Test the entire free fan
# against the same immutable global guard and original face hemispheres.
joint=[]
if guard.remaining and e['budget_or_stop_reason']=='no improving bounded linear clearance proposal':
 fixed_set=set(map(int,locked));bad_faces=np.unique(np.asarray(sorted(guard.remaining)).ravel());free=np.array(sorted(set(map(int,faces[bad_faces].ravel()))-fixed_set),int);planes=[fi for fi in bad_faces if all(int(k)in fixed_set for k in faces[fi])]
 if 0<len(free)<=32 and planes:
  direction=-n[planes].sum(0);direction/=np.linalg.norm(direction)
  for microns in [5.,10.,20.,30.,40.,60.,80.]:
   value=(guard.positions[free]+direction*microns*1e-6-origin).astype(np.float32).astype(float)+origin
   if np.linalg.norm(value-v[free],axis=1).max()>.0001:continue
   proposed=guard.positions.copy();proposed[free]=value;t=proposed[faces];outward=bool(np.all(np.einsum('fi,fi->f',n,np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]))>0))
   token=guard.check_move(free,value)if outward else None;joint.append(dict(microns=microns,free_controls=free.tolist(),outward=outward,remaining_pairs=None if token is None else len(token.remaining)))
   if token is not None and len(token.remaining)<len(guard.remaining):guard.commit(token);break
fixed=guard.positions;e=dict(passed=not guard.remaining,single_vertex=result['evidence'],joint_fan_proposals=joint,global_crossing_guard=guard.evidence());atomic_json(out/'repair.json',e);print('Final repair leaves',len(guard.remaining),'pairs; joint tests',joint,flush=True)
assert e['passed'],'Bounded pinned-exterior chart repair did not pass';assert np.array_equal(fixed[locked],v[locked]);moving=np.flatnonzero(np.any(fixed!=v,axis=1));candidate=copy.deepcopy(a);cp={p['id']:p for p in candidate['parts']};buffers=[];byte_offset=0;chunk=len(a['chunks']);changed=[]
def put(value):
 global byte_offset
 pad=(-byte_offset)%4
 if pad:buffers.append(bytes(pad));byte_offset+=pad
 start=byte_offset;data=value.tobytes();buffers.append(data);byte_offset+=len(data);return start
for pid in ids:
 old=meshes[pid];actual=fixed[maps[pid]]
 if np.array_equal(actual,old.vertices):continue
 p=cp[pid];assert p['geometryOrigin']==origin.tolist();positions=(actual-origin).astype('<f4');assert np.array_equal(positions.astype(float)+origin,actual);layout=dict(positions=put(positions))
 for key,value in old.attributes.items():
  dtype={'normals':'<i2','normalSigns':'i1','uvs':'<f4','colors':'u1'}[key];layout[key]=put(value.astype(dtype))
 layout['indices']=put(old.faces.astype('<u4'));p.update(layout,chunk=chunk);p['sharedChartRepair']=dict(schema='bounded-common-chart-repair/v1',sourceCatalogueSha256=digest_file(path),outerSurfaceLiteral=True,maximumDisplacementMetres=float(np.linalg.norm(actual-old.vertices,axis=1).max()),facesAndAttributesLiteral=True);changed.append(pid)
bank=out/'shared-chart-repair.bin';atomic_bytes(bank,b''.join(buffers));candidate['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size));atomic_json(out/'atlas.json',candidate);atomic_json(out/'packing.json',dict(source_candidate_sha256=digest_file(path),changed_chart_parts=changed,shared_controls_moved=len(moving),maximum_displacement_m=float(np.linalg.norm(fixed-v,axis=1).max()),outer_surface_literal=True,all_nonmembers_literal=True,all_faces_and_attributes_literal=True,bank_sha256=digest_file(bank),elapsed_seconds=time.monotonic()-started,installed=False,anatomical_acceptance=False));print('Prepared',len(changed),'chart parts; maximum move um',float(np.linalg.norm(fixed-v,axis=1).max()*1e6),flush=True)
