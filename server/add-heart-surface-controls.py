"""Bind named cardiac exterior regions to real facets, without duplicate solids."""
from pathlib import Path
import argparse,json,sys
import numpy as np
import igl
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--candidate',required=True);ap.add_argument('--install',required=True);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.candidate;stage=h/args.install;path=out/'atlas.json';a=json.loads(path.read_text());sex=args.sex;marker=a['finiteCardiacIntegration'];parts={p['id']:p for p in a['parts']};host=parts[f'TA98HEART:{sex}:epicardium'];reader=AtlasReader();mesh=reader.read_part(a,host,frame='current-body-metres',catalogue_path=path);features=[];classification={};term_keys=[('base','A12.1.00.002','Base of heart'),('anterior-surface','A12.1.00.003','Anterior surface of heart'),('diaphragmatic-surface','A12.1.00.004','Diaphragmatic surface of heart'),('right-pulmonary-surface','A12.1.00.005','Right pulmonary surface of heart'),('left-pulmonary-surface','A12.1.00.031','Left pulmonary surface of heart'),('apex','A12.1.00.007','Apex of heart')]
# Rebind our own existing controls after a changed host bank; never duplicate
# an earlier control or retain a stale face ordinal from the previous bank.
replaced={g['conceptId']for g in marker['termGroups']if g['key'].startswith('feature:')or g['key']in ['atria','ventricles']}
a['concepts']=[c for c in a['concepts']if c['id']not in replaced]
marker['termGroups']=[g for g in marker['termGroups']if g['conceptId']not in replaced]
if sex=='male':
 free=np.concatenate([np.arange(g['start']//3,(g['start']+g['count'])//3)for g in host['tissueDrawGroups']if not g.get('neighbor')]);tri=mesh.triangles[free];centers=tri.mean(1);normals=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None];names=['right-atrium','right-ventricle','left-atrium','left-ventricle'];distance=[]
 for name in names:
  p=parts[f'TA98HEART:male:myocardium:{name}'];m=reader.read_part(a,p,frame=mesh.frame,catalogue_path=path);distance.append(igl.point_mesh_squared_distance(centers,m.vertices,m.faces)[0])
 tissue=np.argmin(distance,axis=0);directions=np.array([[0.,.25,-1.],[0.,.05,1.],[0.,-.8,-.6],[-1.,0.,0.],[1.,0.,0.]]);directions/=np.linalg.norm(directions,axis=1)[:,None];scores=normals@directions.T;scores[(tissue!=0)&(tissue!=2),0]=-np.inf;scores[tissue!=0,3]=-np.inf;labels=scores.argmax(1)
 valve_ids=next(g['partIds']for g in marker['termGroups']if g['key']=='mitral-valve');base=np.mean([reader.read_part(a,parts[pid],frame=mesh.frame,catalogue_path=path).triangles.mean(axis=(0,1))for pid in valve_ids],axis=0);eligible=free[tissue==3];apical_tri=mesh.triangles[eligible];which=np.unravel_index(np.linalg.norm(apical_tri-base,axis=2).argmax(),apical_tri.shape[:2]);apex_face=int(eligible[which[0]]);apex=mesh.triangles[apex_face,which[1]];apical=free[np.linalg.norm(centers-apex,axis=1)<=.003]
 if apex_face not in apical:apical=np.r_[apical,apex_face]
 for k,(key,term,name)in enumerate(term_keys):
  faces=free[labels==k]if k<5 else apical;assert len(faces);points=mesh.triangles[faces].mean(1);face=int(faces[np.argmin(np.linalg.norm(points-points.mean(0),axis=1))])if k<5 else apex_face;bary=[1/3]*3 if k<5 else [float(i==which[1])for i in range(3)];point=np.array(bary)@mesh.triangles[face];features.append(dict(term=term,conceptId=f'TA98HEART:{sex}:feature:{key}',partIds=[host['id']],representation='surface-region',surface=dict(partId=host['id'],faces=faces.tolist(),sourceBankSha256=marker['coupledBankSha256']),anchor=dict(partId=host['id'],face=face,barycentric=bary,point=point.tolist()),boundaryProvisional=True));classification[term]=dict(facets=len(faces),area_m2=float(np.linalg.norm(np.cross(*(np.moveaxis(mesh.triangles[faces][:,1:]-mesh.triangles[faces][:,:1],1,0))),axis=1).sum()/2))
else:
 male=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text());old=male['cardiacFeatureIntegration'];assert old['bankSha256']==marker['coupledBankSha256'];r,scale=instance_pose(host);morigin=np.array(next(p for p in male['parts']if p['id']==old['hostPartId'])['geometryOrigin']);forigin=np.array(host['geometryOrigin'])
 for g in old['termGroups']:
  g=json.loads(json.dumps(g));g['conceptId']=g['conceptId'].replace(':male:',':female:');g['partIds']=[host['id']];g['surface']['partId']=host['id'];g['anchor']['partId']=host['id'];g['anchor']['point']=((np.array(g['anchor']['point'])-morigin)@r.T*scale+forigin).tolist();features.append(g)
for g in features:
 point=np.array(g['anchor']['barycentric'])@mesh.triangles[g['anchor']['face']];assert np.max(abs(point-g['anchor']['point']))<2e-15;key,term,name=next(row for row in term_keys if row[1]==g['term']);concept=dict(id=g['conceptId'],name=name,elements=[host['id']],ta98Term=term,ta98Kind='exact',partOf=marker['wholeConceptId'],coverageLimitation='Named exterior surface footprint on the stored cardiac coat. Region boundaries are illustrative; this control adds no duplicate tissue solid.');a['concepts'].append(concept);marker['termGroups'].append(dict(key='feature:'+key,name=name,partIds=[host['id']],term=term,conceptId=g['conceptId'],parent='heart',kind='exact'))
for key,term,name,chambers in [('atria','A12.1.00.017','Atria',['right-atrium','left-atrium']),('ventricles','A12.1.00.012','Ventricles',['right-ventricle','left-ventricle'])]:
 ids=[pid for chamber in chambers for pid in next(g['partIds']for g in marker['termGroups']if g['key']==chamber)];cid=f'TA98HEART:{sex}:{key}';a['concepts'].append(dict(id=cid,name=name,elements=ids,ta98Term=term,ta98Kind='exact',partOf=marker['wholeConceptId'],coverageLimitation='Collective control for the separately represented named chambers; no extra tissue is added.'));marker['termGroups'].append(dict(key=key,name=name,partIds=ids,term=term,conceptId=cid,parent='heart',kind='exact'))
 for g in marker['termGroups']:
  if g['key']in chambers:g['parent']=key
for pid,p in parts.items():
 if p.get('cardiacCoupledAttachment'):p['cardiacCoupledAttachment']['interfaces']=[dict(peer=r['second']if r['first']==pid else r['first'],kind=r['kind'],continuous=r['fusion_inferred'])for r in marker['coupledAttachmentSemantics']if pid in [r['first'],r['second']]]
a['cardiacFeatureIntegration']=dict(schema='cardiac-surface-features/v1',hostPartId=host['id'],bankSha256=marker['coupledBankSha256'],termGroups=features);atomic_json(path,a,overwrite=True);atomic_json(out/'surface-features.json',dict(passed=True,sex=sex,bank_sha256=marker['coupledBankSha256'],source_facets=classification,features=features,additional_solid_meshes=0,new_image_assets=0,classification_boundaries_provisional=True,distinct_exact_term_controls=len({g['term']for g in marker['termGroups']if g['kind']=='exact'})));inventory=json.loads((stage/'inventory.json').read_text());termmap={g['term']:g['partIds']for g in marker['termGroups']if g['kind']=='exact'}
for row in inventory['terms']:
 row['candidate_part_ids']=termmap.get(row['authority']['id'],[])
 if row['authority']['id']in {g['term']for g in features}:row['candidate_representation']='source-bound exterior surface footprint; provisional boundary, no extra solid'
atomic_json(stage/'inventory.json',inventory,overwrite=True);spec=json.loads((stage/'install-spec.json').read_text());spec['evidence_files'].append(str(out/'surface-features.json'));atomic_json(stage/'install-spec.json',spec,overwrite=True);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec),overwrite=True);print(sex,len({g['term']for g in marker['termGroups']if g['kind']=='exact'}),'distinct exact TA98 controls, six source-bound exterior features, zero additional solids',flush=True)
