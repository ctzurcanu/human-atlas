"""Apply explicit bounded local impressions to the complete pulmonary cohort.

The same positive-Jacobian fields carry every tissue chart and bronchial
source. Original vessels, muscle, skeleton and cardiac assembly stay fixed.
This only emits a proposal; discrete geometry and whole-context checks follow.
"""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.local_translation_cage import LocalTranslationCage
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--plan',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();started=time.monotonic()
path=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{args.sex}/atlas.json';original=json.loads(path.read_text());a=copy.deepcopy(original);parts={p['id']:p for p in a['parts']};old={p['id']:p for p in original['parts']};plan_path=ROOT/args.plan;plan=json.loads(plan_path.read_text());assert plan['sex']==args.sex and plan['source_catalogue_sha256']==digest_file(path);fields=[LocalTranslationCage(np.array(f['windows_m']),np.array(f['translation_m']))for f in plan['fields']];ids=a['pulmonarySegmentIntegration']['partIds']+a['pulmonaryBranchIntegration']['partIds'];assert len(ids)==46;assert not set(ids)&set(plan['fixed_parts']);out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False);reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres';buffers=[];offset=0;chunk=len(a['chunks']);changed=[];records=[]
def put(value):
 global offset
 pad=(-offset)%4
 if pad:buffers.append(bytes(pad));offset+=pad
 at=offset;value=value.tobytes();buffers.append(value);offset+=len(value);return at
for pid in ids:
 p=parts[pid];assert not p.get('geometryInstance')and p.get('geometryScale',1)==1 and not p.get('sourceOffset');mesh=reader.read_part(original,old[pid],frame=frame,catalogue_path=path);v=mesh.vertices.copy();j=np.broadcast_to(np.eye(3),(len(v),3,3)).copy()
 for field in fields:
  next_j=field.jacobians(v);v=field.map(v);j=next_j@j
 active=np.any(v!=mesh.vertices,axis=1)
 if not np.any(active):continue
 origin=np.array(p.get('geometryOrigin',[0,0,0]));pos=(v-origin).astype('<f4');actual=pos.astype(float)+origin;row=dict(positions=put(pos))
 for key,value in mesh.attributes.items():
  if key=='normals':
   n=np.linalg.solve(j.transpose(0,2,1),value.astype(float)[...,None])[...,0];length=np.linalg.norm(n,axis=1);assert np.all(length>0);n=np.rint(n/length[:,None]*32767).astype('<i2');n[~active]=value[~active].astype('<i2');value=n
  else:value=value.astype({'uvs':'<f4','normalSigns':'i1','colors':'u1'}[key])
  row[key]=put(value)
 row['indices']=put(mesh.faces.astype('<u4'));p.update(row,chunk=chunk,vertexCount=len(mesh.vertices),indexCount=int(mesh.faces.size),triangles=len(mesh.faces),bounds=np.array([actual.min(0),actual.max(0)]).tolist());p['pulmonaryContactImpressions']=dict(schema='bounded-local-pulmonary-impressions/v1',sourceCatalogueSha256=digest_file(path),planSha256=digest_file(plan_path),fixedParts=plan['fixed_parts'],commonPulmonaryCohort=ids,minimumContinuousJacobian=float(np.prod([f.determinant_lower_bound for f in fields])),noTissueRemoved=True);changed.append(pid);records.append(dict(id=pid,maximum_displacement_m=float(np.linalg.norm(actual-mesh.vertices,axis=1).max()),moved_controls=int(active.sum()),upper_airway_controls_changed=int(np.sum(active&(mesh.vertices[:,1]>=1.475)))if pid in a['pulmonaryBranchIntegration']['partIds']else None))
bank=out/'pulmonary-contact-impressions.bin';atomic_bytes(bank,b''.join(buffers));a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size));assert all(p==old[pid]for pid,p in parts.items()if pid not in changed)
for key in ['pulmonarySegmentIntegration','pulmonaryBranchIntegration']:a[key]['contactImpressions']=dict(schema='bounded-local-pulmonary-impressions/v1',planSha256=digest_file(plan_path),fixedParts=plan['fixed_parts'],bankSha256=digest_file(bank),anatomicalAcceptance=False)
atomic_json(out/'atlas.json',a);atomic_json(out/'packing.json',dict(source_catalogue_sha256=digest_file(path),plan=str(plan_path),plan_sha256=digest_file(plan_path),moving_parts=changed,field_applied_to_cohort=ids,all_nonmembers_literal=True,source_banks_literal=True,replacement_bank=str(bank),replacement_bank_sha256=digest_file(bank),records=records,skeleton_unchanged=True,reference_vessels_unchanged=True,neck_muscle_unchanged=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-started));print('Prepared',args.sex,len(changed),'coordinated pulmonary parts; upper airway changes',sum(r['upper_airway_controls_changed']or 0 for r in records),flush=True)
