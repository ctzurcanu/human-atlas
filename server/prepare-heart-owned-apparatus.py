"""Stage a whole checked apparatus and its exact reciprocal interface controls."""
from pathlib import Path
import argparse,copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--bank',required=True);ap.add_argument('--candidate',required=True);ap.add_argument('--install',required=True);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';bank=h/args.bank;e=json.loads((bank/'emission.json').read_text());assert e['passed']and digest_file(bank/'cells.bin')==e['bank_sha256']and json.loads((bank/'checked-geometry.json').read_text())['passed'];sex=args.sex;model='local-ta98'if sex=='male'else'local-female-ta98';current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';a=json.loads(current.read_text());marker=a['finiteCardiacIntegration'];parts={p['id']:p for p in a['parts']};reader=AtlasReader();origin=np.array(e['geometry_origin']);chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=e['bytes']));male_source={v['record']['id']:v['record']for s in e['sources']for v in s['source_meshes'].values()};mapping={}
for pid in e['parts']:
 if sex=='male':mapping[pid]=pid
 elif pid.startswith('TA98HEART:male:chordae:'):mapping[pid]=pid.replace(':male:',':female:')
 else:
  found=[p['id']for p in a['parts']if p.get('cardiacTemplateSourcePartId')==pid];assert len(found)==1,(pid,found);mapping[pid]=found[0]
common_origin=origin if sex=='male'else np.array(parts[marker['wallPartIds'][0]]['geometryOrigin']);transfers=[]
for pid,layout in e['parts'].items():
 p=parts[mapping[pid]];old=copy.deepcopy(p);source=male_source[pid];assert not p.get('sourceOffset')and float(source.get('geometryScale',1.))==1.;pose=instance_pose(source);assert pose is None or(np.array_equal(pose[0],np.eye(3))and pose[1]==1.)
 if sex=='female':
  r,scale=instance_pose(p);declared=np.array(p.get('geometryOrigin',[0.,0.,0.]))+(origin-np.array(source.get('geometryOrigin',[0.,0.,0.])))@r.T*scale;error=float(np.linalg.norm(declared-common_origin));assert error<2e-15,(pid,error);transfers.append(dict(part_id=p['id'],common_origin_residual_m=error,scale=scale,rotation=p['geometryInstance']['rotation']))
 for key in ['positions','normals','indices','normalSigns','uvs','vertexCount','indexCount','triangles','tissueDrawGroups']:p.pop(key,None)
 p.update({k:copy.deepcopy(v)for k,v in layout.items()if k not in ['id','bounds','geometryOrigin','geometryInstance','geometryScale','sourceOffset']});p.update(chunk=chunk,geometryOrigin=common_origin.tolist(),pairedTissue=True,sectionSolid=True,tissueWallAssembly=f'heart:{sex}',tissueWallLayer='chordae'if p['id']in marker['chordalPartIds']else'papillary'if p['id']in marker['papillaryPartIds']else'valve')
 if sex=='male':p['geometryInstance']=dict(rotation=[0.,0.,0.,1.],scale=1.)
 for group in p['tissueDrawGroups']:
  if group.get('neighbor'):group['neighbor']=mapping[group['neighbor']]
 if p.get('cardiacChordalPartition'):
  c=p['cardiacChordalPartition'];c.update(papillaryInterfaceVerified=True,sharedPapillaryInterface=True,bankSha256=e['bank_sha256']);p['coverageLimitation']='Native leaflet and chordal source tissue with shared finite interfaces at the leaflet neck and papillary insertion. Computational ownership boundaries, complete fine chordal distribution and annular anatomy remain unverified.'
 if p['id']in marker['papillaryPartIds']:p['cardiacPapillaryInterface']=dict(schema='native-papillary-interface/v1',bankSha256=e['bank_sha256'],illustrativeBasalOwnership=True)
 p['cardiacApparatusOwnership']=dict(schema='native-apparatus-ownership/v1',sourcePartId=pid,bankSha256=e['bank_sha256'],sourceUnionPreserved=True,maximumControlRoundingMeters=e['maximum_control_rounding_m']);mesh=reader.read_part(a,p,frame='model:'+model+':assembled-metres');p['bounds']=mesh.bounds.tolist();a['materials'][p['material']]['pairedNormalSigns']=True
required=[]
for s in e['sources']:
 ids={role:mapping[r['id']]for role,r in s['parts'].items()};required.extend([dict(first=ids['leaflet'],second=ids['chordae'],kind='native-leaflet-neck'),dict(first=ids['chordae'],second=ids['papillary'],kind='native-papillary-insertion')])
interfaces=[];semantics=[]
for pid in mapping.values():
 p=parts[pid]
 for group in p['tissueDrawGroups']:
  neighbor=group.get('neighbor')
  if not neighbor or pid>=neighbor:continue
  peers=[g for g in parts[neighbor]['tissueDrawGroups']if g.get('neighbor')==pid];assert peers;interfaces.append(dict(kind='shared-interface',first=pid,second=neighbor,parameters=dict(first_faces=list(range(group['start']//3,(group['start']+group['count'])//3)),second_faces=[i for g in peers for i in range(g['start']//3,(g['start']+g['count'])//3)])));declared=next((r['kind']for r in required if{r['first'],r['second']}=={pid,neighbor}),None);semantics.append(dict(first=pid,second=neighbor,kind=declared or('coaptation-sliding'if pid in marker['valvePartIds']and neighbor in marker['valvePartIds']else'regional-muscle-boundary'if pid in marker['papillaryPartIds']and neighbor in marker['papillaryPartIds']else'provisional-source-contact'),fusion_inferred=bool(declared)))
assert len([s for s in semantics if s['kind']=='native-papillary-insertion'])==5
remaining='Finite chambers, septal regions, wall coats and native valves with five conforming leaflet–chordal–papillary insertions. Septal, basal-muscle and strand-neck ownership boundaries and lining thickness are illustrative. Complete fine chordal distribution, annular attachment, conduction tissues, pericardial detail and distal recipient vessel joins remain unfinished.';marker.update(remainingLimits=remaining,apparatusBankSha256=e['bank_sha256'],apparatusInterfaceSemantics=semantics);next(c for c in a['concepts']if c['id']==marker['wholeConceptId'])['coverageLimitation']=remaining;out=h/args.candidate;out.mkdir(exist_ok=False);atomic_json(out/'atlas.json',a);atomic_json(out/'transfer.json',dict(passed=True,shared_bank_sha256=e['bank_sha256'],declared_positive_instances=transfers,new_geometry_banks=0 if sex=='female'else 1,all_interfaces_in_one_common_frame=True));atomic_json(out/'candidate.json',dict(passed=True,sex=sex,source_catalogue_sha256=digest_file(current),catalogue_sha256=digest_file(out/'atlas.json'),changed_parts=list(mapping.values()),interfaces=semantics,remaining_limits=remaining,installed=False));stage=h/args.install;stage.mkdir(exist_ok=False);inventory=copy.deepcopy(json.loads((h/f'source-v1-{sex}/inventory.json').read_text()));inventory['source_catalogue_sha256']=digest_file(current);terms={}
for g in marker['termGroups']:
 if g['kind']=='exact':terms.setdefault(g['term'],[]).extend(g['partIds'])
for row in inventory['terms']:row['candidate_part_ids']=list(dict.fromkeys(terms.get(row['authority']['id'],[])));row['candidate_representation']='stored finite tissue or collective control'if row['candidate_part_ids']else None
atomic_json(stage/'inventory.json',inventory);spec=dict(model=model,proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(current),changed_parts=list(mapping.values()),closed_parts=list(mapping.values()),require_closed=False,max_pairs=1500000,max_tests=16000000,inventory=str(stage/'inventory.json'),tissue_audit_layouts=[],interface_checks=interfaces,assembly=dict(whole_organ=marker['wholeConceptId'],part_ids=next(c['elements']for c in a['concepts']if c['id']==marker['wholeConceptId']),fixed=['skeleton','lungs','diaphragm','cardiac chamber walls','vessel roots','source apparatus union and course'],adjoining_interfaces=semantics,remaining_limits=remaining),evidence_files=[str(bank/'emission.json'),str(bank/'checked-geometry.json'),str(bank/'storage-preflight.json'),str(h/e['source_stage']/'ownership.json'),str(out/'candidate.json'),str(out/'transfer.json')]);atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec));print(sex,len(mapping),'finite parts,',len(interfaces),'reciprocal charts; five papillary insertions',flush=True)
