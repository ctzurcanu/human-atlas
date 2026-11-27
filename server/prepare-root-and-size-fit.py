"""Coherent cardiac size correction and source-bound great-root placement.

Keep the entire skeleton, chest muscles and all finite tissue ownership fixed.
Only declared cardiac/pulmonary cohorts and the local aortic arterial field move.
No new solid, cap, vascular wall or certified recipient join is invented.
"""
from pathlib import Path
from collections import Counter,defaultdict
import argparse,copy,json,sys,time
import numpy as np
from scipy.spatial.transform import Rotation
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader
from geometry_tools.geometry_instance import instance_pose
from geometry_tools.common import atomic_json,atomic_bytes,digest_file
from geometry_tools.algorithms.local_translation_cage import LocalTranslationCage
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--stage',required=True);ap.add_argument('--arch-dx-mm',type=float,default=10.);ap.add_argument('--arch-dz-mm',type=float,default=0.);ap.add_argument('--posterior-arch',action='store_true');args=ap.parse_args();start=time.monotonic();sex=args.sex;path=ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json';original=json.loads(path.read_text());a=copy.deepcopy(original);parts={p['id']:p for p in a['parts']};bp={p['id']:p for p in original['parts']};reader=AtlasReader();frame=f'model:local-{sex}-ta98:assembled-metres';out=ROOT/'.local-models/geometry-workbench/heart'/args.stage;out.mkdir(exist_ok=False);forest_path=ROOT/f'reports/structural-dependencies/current/{sex}-forest.json';forest=json.loads(forest_path.read_text());assert forest['sourceSha256']==digest_file(path);heart=[g for root in forest['forest']for g in root['children']if g['id']=='assembly:heart'];assert len(heart)==1;heart_ids=heart[0]['activePartIds']
def read(pid,new=True):return reader.read_part(a if new else original,parts[pid]if new else bp[pid],frame=frame,catalogue_path=path)
scale=.95 if sex=='male'else .9;pivot=np.mean([read(f'TA98HEART:{sex}:fibrous:{s}-annulus').vertices.mean(0)for s in ['left','right']],axis=0);mat=np.eye(4);mat[:3,:3]*=scale;mat[:3,3]=pivot*(1-scale)
def point(v):return (np.asarray(v)*scale+mat[:3,3]).tolist()
for pid in heart_ids:
 p=parts[pid];prior=instance_pose(p);rot,s=prior if prior is not None else(np.eye(3),1.);p['geometryOrigin']=point(p.get('geometryOrigin',[0,0,0]));p['geometryInstance']=dict(rotation=Rotation.from_matrix(rot).as_quat().tolist(),scale=s*scale);p['cardiacThoracicFit']=dict(schema='cardiac-size-and-root-packing/v3',sourceCatalogueSha256=digest_file(path),matrix=mat.tolist(),relativeScale=scale);p['bounds']=read(pid).bounds.tolist()
for group in a.get('cardiacFeatureIntegration',{}).get('termGroups',[]):group['anchor']['point']=point(group['anchor']['point'])
a['finiteCardiacIntegration']['geometryInstance']=copy.deepcopy(parts[a['finiteCardiacIntegration']['wallPartIds'][0]]['geometryInstance']);a['finiteCardiacIntegration']['thoracicFit']=dict(schema='cardiac-size-and-root-packing/v3',sourceCatalogueSha256=digest_file(path),matrix=mat.tolist(),movingPartIds=heart_ids,clinicalSizeEstablished=False,independentRecipientJoinsVerified=False)
# One field carries every pulmonary chart at the male apex; the cervical
# muscle and all upper airway controls remain literal.
apex=LocalTranslationCage(np.array([[.006,.01,.022,.047],[1.49,1.5,1.525,1.54],[-.04,-.028,-.006,.007]]),np.array([.004,0,0]))
assert 0<=args.arch_dx_mm<=12 and 0<=args.arch_dz_mm<=20
arch_windows=[[0,.03,.045,.08],[1.31,1.335,1.402,1.42],[-.10,-.02,.042,.06]]if args.posterior_arch else[[0,.02,.045,.08],[1.345,1.364,1.402,1.42],[-.035,-.02,.042,.06]]
arch=LocalTranslationCage(np.array(arch_windows),np.array([-args.arch_dx_mm/1000,0,-args.arch_dz_mm/1000]),minimum_jacobian=.1)
pulmonary=a['pulmonarySegmentIntegration']['partIds']+a['pulmonaryBranchIntegration']['partIds'];hidden=set(a.get('defaultHidden',[]));arterial=[p['id']for p in a['parts']if p['system']=='arterial'and p['id']not in hidden and p['id']not in heart_ids and not p.get('defaultHidden')and not p.get('suppressed')]
if sex=='female':arterial=sorted(set(arterial+['LOCAL:female:aortic_arch:1513']))
pool=pulmonary if sex=='male'else arterial;field=apex if sex=='male'else arch;buffers=[];offset=0;chunk=len(a['chunks']);baked=[];origin=np.array([0.,1.3521 if sex=='male'else 1.25,0.])
def put(value):
 global offset
 pad=(-offset)%4
 if pad:buffers.append(bytes(pad));offset+=pad
 at=offset;value=value.tobytes();buffers.append(value);offset+=len(value);return at
for pid in pool:
 m=read(pid);w,_=field.weights(m.vertices)
 if not np.any(w):continue
 v=field.map(m.vertices);j=field.jacobians(m.vertices);pos=(v-origin).astype('<f4');actual=pos.astype(float)+origin;p=parts[pid];row=dict(positions=put(pos))
 for key,value in m.attributes.items():
  if key=='normals':
   n=np.linalg.solve(j.transpose(0,2,1),value.astype(float)[...,None])[...,0];n/=np.linalg.norm(n,axis=1)[:,None];words=np.rint(n*32767).astype('<i2');words[w==0]=value[w==0].astype('<i2');value=words
  else:value=value.astype({'normalSigns':'i1','uvs':'<f4','colors':'u1'}[key])
  row[key]=put(value)
 row['indices']=put(m.faces.astype('<u4'))
 for k in ['geometryScale','sourceOffset','geometryInstance']:p.pop(k,None)
 p.update(row,chunk=chunk,geometryOrigin=origin.tolist(),vertexCount=len(m.vertices),indexCount=int(m.faces.size),triangles=len(m.faces),bounds=np.array([actual.min(0),actual.max(0)]).tolist());p['localRootPacking']=dict(schema='source-bound-local-packing/v1',sourceCatalogueSha256=digest_file(path),field='left-apex'if sex=='male'else'aortic-arch-arterial-cohort',windowsMetres=field.windows.tolist(),translationMetres=field.translation.tolist(),minimumContinuousJacobian=field.determinant_lower_bound);baked.append(pid)
if baked:
 bank=out/'local-packing.bin';atomic_bytes(bank,b''.join(buffers));a['chunks'].append(dict(url='/local-models/'+str(bank.relative_to(ROOT/'.local-models')),bytes=bank.stat().st_size))
# The retained open pulmonary-trunk reference is registered at the active
# valve, rather than left at the old unscaled source heart's position.
pid='ZA:Pulmonary trunk'if sex=='male'else'LOCAL:female:heart-template:za-pulmonary-trunk';m=read(pid);edges=Counter(tuple(sorted(map(int,e)))for f in m.faces for e in zip(f,np.roll(f,-1)));graph=defaultdict(set)
for(i,j),n in edges.items():
 if n==1:graph[i].add(j);graph[j].add(i)
remaining=set(graph);rings=[]
while remaining:
 stack=[next(iter(remaining))];ids=[]
 while stack:
  k=stack.pop()
  if k not in remaining:continue
  remaining.remove(k);ids.append(k);stack.extend(graph[k])
 rings.append(np.array(ids))
assert len(rings)==2 and all(len(x)==23 for x in rings);valves=[i for i in a['finiteCardiacIntegration']['valvePartIds']if 'pulmonary'in parts[i]['name'].lower()];valve=np.concatenate([read(i).vertices for i in valves]).mean(0);prox=min(rings,key=lambda ids:np.linalg.norm(m.vertices[ids].mean(0)-valve));source_centre=m.vertices[prox].mean(0);delta=valve-source_centre;root_ids=[pid]
if sex=='male':root_ids.append('ZA:Bifurcation of pulmonary trunk')
for rid in root_ids:
 p=parts[rid];p['geometryOrigin']=(np.array(p.get('geometryOrigin',[0,0,0]))+delta).tolist();p['bounds']=read(rid).bounds.tolist();p['cardiacRootRegistration']=dict(schema='pulmonary-root-valve-centre-registration/v1',sourceCatalogueSha256=digest_file(path),translationMetres=delta.tolist(),sourceMouthCentreMetres=source_centre.tolist(),valveCentreMetres=valve.tolist(),sourceFacesAndMouthShapeLiteral=True,exactAnnularJoinEstablished=False,distalRecipientJoinEstablished=False)
changed=sorted(pid for pid,p in parts.items()if p!=bp[pid]);assert all(p==bp[pid]for pid,p in parts.items()if p['system']=='skeletal');assert all(p==bp[pid]for pid,p in parts.items()if p['system']=='muscular');assert not set(baked)&set(heart_ids)
atomic_json(out/'atlas.json',a);atomic_json(out/'packing.json',dict(source_catalogue_sha256=digest_file(path),moving_parts=changed,heart_closure=heart_ids,heart_relative_scale=scale,heart_matrix=mat.tolist(),baked_field_parts=baked,pulmonary_field_cohort=pulmonary if sex=='male'else[],arterial_field_cohort=arterial if sex=='female'else[],pulmonary_root_parts=root_ids,root_translation_m=delta.tolist(),root_mouth_shape_and_faces_literal=True,skeleton_unchanged=True,muscles_unchanged=True,all_other_records_literal=True,remaining_limits='Root-centre registration is not an exact finite annular/recipient wall join; source vascular surfaces remain open references.',installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-start));print(sex,'prepared',len(changed),'changes; root translation mm',(delta*1000).tolist(),'field parts',len(baked),flush=True)
