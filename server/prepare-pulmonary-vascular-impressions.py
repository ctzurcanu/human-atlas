"""Local medial lung impressions; reference vessels and skeleton remain fixed."""
from pathlib import Path
import argparse,copy,json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.local_translation_cage import LocalTranslationCage
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--catalogue',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--amplitude-mm',type=float,required=True);args=ap.parse_args();start=time.monotonic();assert 0<args.amplitude_mm<=12
path=ROOT/args.catalogue;original=json.loads(path.read_text());a=copy.deepcopy(original);parts={p['id']:p for p in a['parts']};old={p['id']:p for p in original['parts']};reader=AtlasReader();frame=f'model:local-{args.sex}-ta98:assembled-metres';out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False)
refs=['LOCAL:female:aortic_arch:1513','LOCAL:female:heart-template:za-pulmonary-trunk']if args.sex=='female'else['ZA:Pulmonary trunk'];bounds=np.array([parts[i]['bounds']for i in refs]);lo=bounds[:,0].min(0);hi=bounds[:,1].max(0);windows=np.array([[0.,.015,.065,.095],[lo[1]-.015,lo[1],hi[1],hi[1]+.015],[lo[2]-.015,lo[2],hi[2],hi[2]+.015]]);cage=LocalTranslationCage(windows,np.array([args.amplitude_mm/1000,0,0]));ids=a['pulmonarySegmentIntegration']['partIds']+a['pulmonaryBranchIntegration']['partIds'];buffers=[];offset=0;chunk=len(a['chunks']);changed=[];records=[]
def put(array):
 global offset
 pad=(-offset)%4
 if pad:buffers.append(bytes(pad));offset+=pad
 result=offset;data=array.tobytes();buffers.append(data);offset+=len(data);return result
for pid in ids:
 p=parts[pid];assert not p.get('geometryInstance')and p.get('geometryScale',1)==1 and not p.get('sourceOffset');m=reader.read_part(original,old[pid],frame=frame,catalogue_path=path);weight,_=cage.weights(m.vertices)
 if not np.any(weight):continue
 v=cage.map(m.vertices);j=cage.jacobians(m.vertices);origin=np.array(p.get('geometryOrigin',[0,0,0]));pos=(v-origin).astype('<f4');row=dict(positions=put(pos));active=weight>0
 for k,value in m.attributes.items():
  if k=='normals':
   n=np.linalg.solve(j.transpose(0,2,1),value.astype(float)[...,None])[...,0];length=np.linalg.norm(n,axis=1);assert np.all(length>0);n=np.rint(n/length[:,None]*32767).astype('<i2');n[~active]=value[~active].astype('<i2');value=n
  elif k=='uvs':value=value.astype('<f4')
  elif k=='normalSigns':value=value.astype('i1')
  elif k=='colors':value=value.astype('u1')
  else:raise ValueError(k)
  row[k]=put(value)
 row['indices']=put(m.faces.astype('<u4'));p.update(row,chunk=chunk,vertexCount=len(v),indexCount=int(m.faces.size),triangles=len(m.faces),bounds=np.array([(pos.astype(float)+origin).min(0),(pos.astype(float)+origin).max(0)]).tolist());p['pulmonaryVascularImpression']=dict(schema='local-medial-pulmonary-impression/v1',sourceCatalogueSha256=digest_file(path),fixedVesselParts=refs,windowsMetres=windows.tolist(),translationMetres=cage.translation.tolist(),minimumContinuousJacobian=cage.determinant_lower_bound,allBronchopulmonaryPartsUseOneField=True);changed.append(pid);records.append(dict(id=pid,maximum_displacement_m=float(np.linalg.norm(v-m.vertices,axis=1).max()),unmodified_controls=int((~active).sum()),minimum_sample_jacobian=float(np.linalg.det(j).min())))
assert changed;bank=out/'pulmonary-impressions.bin';atomic_bytes(bank,b''.join(buffers));a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
for p in a['parts']:
 if p['id']not in changed:assert p==old[p['id']]
for k in ['pulmonarySegmentIntegration','pulmonaryBranchIntegration']:a[k]['vascularImpressions']=dict(schema='local-medial-pulmonary-impression/v1',bankSha256=digest_file(bank),fixedVesselParts=refs,sourceCatalogueSha256=digest_file(path),anatomicalAcceptance=False)
atomic_json(out/'atlas.json',a);atomic_json(out/'packing.json',dict(source_catalogue_sha256=digest_file(path),moving_parts=changed,field_applied_to_cohort=ids,all_nonmembers_literal=True,source_banks_literal=True,replacement_bank=str(bank),replacement_bank_sha256=digest_file(bank),records=records,windows_m=windows.tolist(),translation_m=cage.translation.tolist(),minimum_continuous_jacobian=cage.determinant_lower_bound,reference_vessels_unchanged=True,skeleton_unchanged=True,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start));print('Prepared local impressions in',len(changed),'pulmonary parts; vessels and skeleton literal',flush=True)
