"""Stage checked coupled heart tissues in one body frame, with explicit contacts."""
from pathlib import Path
import argparse,copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.geometry_instance import instance_pose
ap=argparse.ArgumentParser();ap.add_argument('--sex',required=True,choices=['male','female']);ap.add_argument('--bank',required=True);ap.add_argument('--candidate',required=True);ap.add_argument('--install',required=True);ap.add_argument('--transfer-source');ap.add_argument('--audit-checkpoints');args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';bank=h/args.bank;e=json.loads((bank/'emission.json').read_text());checks=json.loads((bank/'checked-geometry.json').read_text());assert checks['passed']and checks['complete_bank']and checks['bank_sha256']==e['bank_sha256']and digest_file(bank/'cells.bin')==e['bank_sha256'];source=h/e['source_stage'];meta=json.loads((source/'ownership.json').read_text());sex=args.sex;model='local-ta98'if sex=='male'else'local-female-ta98';current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';a=json.loads(current.read_text());marker=a['finiteCardiacIntegration'];male_current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';male_baseline=h/args.transfer_source if args.transfer_source else male_current;male=json.loads(male_baseline.read_text());mp={p['id']:p for p in male['parts']};mm=male['finiteCardiacIntegration'];parts={p['id']:p for p in a['parts']};mapping={}
for key in ['wallPartIds','rootWallPartIds','rootPartIds','valvePartIds','chordalPartIds','papillaryPartIds']:
 assert len(mm[key])==len(marker[key]);mapping.update(zip(mm[key],marker[key]))
new_regions=set()
for region in meta.get('annular_regions',[]):
 pid=region['region_id'];target=pid.replace(':male:',':'+sex+':');mapping[pid]=target
 if target in parts:continue
 template=copy.deepcopy(parts[f'TA98HEART:{sex}:fibrous:interventricular-membranous-septum'])
 name=region['side'].capitalize()+' fibrous ring';template.update(id=target,conceptId=target,name=name,ta98Term=region['term'],tissueWallLayer='fibrous-septum',sourceId='Heart:measured-annular-partition',coverageLimitation='Finite connective support partitioned within the existing myocardium at measured peripheral leaflet attachments. Histological extent and discontinuities remain provisional.',cardiacAnnularPartition=dict(schema='owned-annular-region/v1',term=region['term'],sourceRegion=region,sourceBankSha256=meta['annular_base_bank_sha256']))
 template['provenance']=dict(label='Local finite annular construction',url='https://github.com/Z-Anatomy/Models-of-human-anatomy',detail='Derived myocardial host partition at measured leaflet hinges. Supplied local model consulted for shape only; no copied reference mesh or texture.')
 a['parts'].append(template);parts[target]=template;new_regions.add(target);marker['wallPartIds'].append(target)
 group=dict(key='fibrous-'+region['side']+'-ring',name=name,partIds=[target],term=region['term'],conceptId=target,parent='heart',kind='exact');marker['termGroups'].append(group)
 a['concepts'].append(dict(id=target,name=name,elements=[target],ta98Term=region['term'],ta98Kind='exact',partOf=marker['wholeConceptId'],coverageLimitation=template['coverageLimitation']))
 whole=next(c for c in a['concepts']if c['id']==marker['wholeConceptId']);whole['elements']=list(dict.fromkeys(whole['elements']+[target]));next(g for g in marker['termGroups']if g['key']=='heart')['partIds']=whole['elements']
assert set(e['parts'])<=set(mapping);origin=np.array(e['geometry_origin']);common=origin if sex=='male'else np.array(parts[marker['wallPartIds'][0]]['geometryOrigin']);reader=AtlasReader();old_male=AtlasReader();transfers=[];chunk=len(a['chunks']);a['chunks'].append(dict(url='/local-models/'+str((bank/'cells.bin').relative_to(ROOT/'.local-models')),bytes=e['bytes']))
for pid,layout in e['parts'].items():
 target=mapping[pid];p=parts[target];old=copy.deepcopy(p)
 if sex=='female' and target not in new_regions:
  r,scale=instance_pose(old);actual=reader.read_part(a,old,frame='recipient-source',catalogue_path=current)
  if args.transfer_source:src=old_male.read_part(male,mp[pid],frame='male-source',catalogue_path=male_baseline)
  else:source_records=json.loads((h/'coupled-attachments-v249/source.json').read_text())['sources'];src=old_male.reload(source_records[pid],frame='male-source')
  expected=(src.triangles-origin)@r.T*scale+common;error=float(np.abs(actual.triangles-expected).max());assert actual.triangles.shape==expected.shape and error<2e-15,(pid,error);transfers.append(dict(part=target,maximum_source_transfer_error_m=error,scale=scale))
 elif sex=='female':
  r,scale=instance_pose(old);transfers.append(dict(part=target,new_owned_region=True,literal_shared_bank=True,scale=scale,proper_rotation=bool(np.linalg.det(r)>0)))
 for key in ['positions','normals','normalSigns','uvs','indices','vertexCount','indexCount','triangles','tissueDrawGroups']:p.pop(key,None)
 p.update({k:copy.deepcopy(v)for k,v in layout.items()if k not in ['id','bounds','geometryOrigin','geometryInstance','sourceOffset','geometryScale']});p.update(chunk=chunk,geometryOrigin=common.tolist(),pairedTissue=True,sectionSolid=True)
 if sex=='male':p['geometryInstance']=dict(rotation=[0.,0.,0.,1.],scale=1.)
 for group in p['tissueDrawGroups']:
  if group.get('neighbor'):group['neighbor']=mapping[group['neighbor']]
 p['cardiacCoupledAttachment']=dict(schema='coupled-cardiac-attachment/v1',bankSha256=e['bank_sha256'],sourcePartId=pid,sourceUnionPreservedBeforeBoundedStorage=True)
 if p.get('cardiacApparatusOwnership'):p['cardiacApparatusOwnership']['bankSha256']=e['bank_sha256']
 if p.get('cardiacChordalPartition'):p['cardiacChordalPartition']['bankSha256']=e['bank_sha256']
 if p.get('cardiacPapillaryInterface'):p['cardiacPapillaryInterface'].update(bankSha256=e['bank_sha256'],wallInterfaceVerified=True)
 if pid=='ZA:Superior vena cava':
  assert json.loads((bank/'caval-path.json').read_text())['passed'];material=f'HEART:{sex}:caval-wall';a['materials'][material]=copy.deepcopy(a['materials'][p['material']]);p.update(material=material,tissueWallAssembly=f'heart:{sex}',tissueWallLayer='caval-wall',cardiacRootWall=dict(schema='finite-caval-root-wall/v1',bankSha256=e['bank_sha256'],sourceMalePartId=pid,distalRecipientJoinVerified=False,histologicalThicknessVerified=False,illustrativeWallThicknessMeters=.0006,atriumEntryTracksVerified=True));marker['rootWallPartIds']=list(dict.fromkeys(marker['rootWallPartIds']+[target]))
 p['bounds']=reader.read_part(a,p,frame='model:'+model+':assembled-metres').bounds.tolist();a['materials'][p['material']]['pairedNormalSigns']=True
interfaces=[];semantic=[];old_semantics=marker.get('coupledAttachmentSemantics',marker['apparatusInterfaceSemantics']);new_ids={mapping[pid]for pid in e['parts']};seen_pairs=set()
for pid in new_ids:
 for group in parts[pid]['tissueDrawGroups']:
  peer=group.get('neighbor')
  if not peer or pid>=peer or (pid,peer)in seen_pairs:continue
  seen_pairs.add((pid,peer))
  matching=[g for g in parts[peer]['tissueDrawGroups']if g.get('neighbor')==pid];assert matching
  first_matching=[g for g in parts[pid]['tissueDrawGroups']if g.get('neighbor')==peer]
  interfaces.append(dict(kind='shared-interface',first=pid,second=peer,parameters=dict(first_faces=[i for g in first_matching for i in range(g['start']//3,(g['start']+g['count'])//3)],second_faces=[i for g in matching for i in range(g['start']//3,(g['start']+g['count'])//3)])))
  prior=next((r for r in old_semantics if{r['first'],r['second']}=={pid,peer}),None)
  if prior:semantic.append(prior);continue
  first,second=parts[pid],parts[peer];is_wall=lambda p:p['id']in marker['wallPartIds'];host=first if is_wall(first)else second;member=second if is_wall(first)else first
  kind='cardiac-source-contact';fusion=False
  if is_wall(host)and member['id']in marker['papillaryPartIds']and host.get('tissueWallLayer')in ['myocardium','endocardium']:kind='papillary-wall-insertion';fusion=True
  elif is_wall(host)and member['id']in marker['valvePartIds']and host.get('tissueWallLayer')in ['myocardium','endocardium']:kind='leaflet-host-attachment';fusion=True
  elif is_wall(host)and member.get('cardiacRootWall',{}).get('schema')=='finite-caval-root-wall/v1':kind='vessel-host-attachment';fusion=True
  if any(p.get('cardiacAnnularPartition')for p in [first,second]):
   if first['id']in marker['valvePartIds']or second['id']in marker['valvePartIds']:kind='annular-leaflet-insertion';fusion=True
   elif is_wall(first)and is_wall(second):kind='annular-myocardial-insertion';fusion=True
  semantic.append(dict(first=pid,second=peer,kind=kind,fusion_inferred=fusion,regional_boundary_provisional=True))
for pid in new_ids:parts[pid]['cardiacCoupledAttachment']['interfaces']=[dict(peer=r['second']if r['first']==pid else r['first'],kind=r['kind'],continuous=r['fusion_inferred'])for r in semantic if pid in [r['first'],r['second']]]
remaining='Finite chambers, wall coats, source leaflet and papillary attachments with reciprocal tissue interfaces. Regional septal/annular boundaries, fine chordal distribution, full LV papillary homology, fibrous-ring histology and lining thickness remain provisional. Conduction, pericardial detail and distal vessel joins remain unfinished.';marker.update(coupledBankSha256=e['bank_sha256'],wallBankSha256=e['bank_sha256'],apparatusBankSha256=e['bank_sha256'],coupledAttachmentSemantics=semantic,remainingLimits=remaining);next(c for c in a['concepts']if c['id']==marker['wholeConceptId'])['coverageLimitation']=remaining
out=h/args.candidate;out.mkdir(exist_ok=False);atomic_json(out/'atlas.json',a)
if sex=='male':atomic_json(out/'male-source.json',male)
atomic_json(out/'transfer.json',dict(passed=True,shared_bank_sha256=e['bank_sha256'],common_origin=common.tolist(),source_transfers=transfers));stage=h/args.install;stage.mkdir(exist_ok=False);inventory=copy.deepcopy(json.loads((h/f'source-v1-{sex}/inventory.json').read_text()));inventory['source_catalogue_sha256']=digest_file(current);terms={}
for g in marker['termGroups']:
 if g['kind']=='exact':terms.setdefault(g['term'],[]).extend(g['partIds'])
for row in inventory['terms']:row['candidate_part_ids']=list(dict.fromkeys(terms.get(row['authority']['id'],[])))
atomic_json(stage/'inventory.json',inventory);spec=dict(model=model,proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(current),changed_parts=sorted(new_ids),closed_parts=sorted(new_ids),require_closed=False,max_pairs=1500000,max_tests=16000000,inventory=str(stage/'inventory.json'),tissue_audit_layouts=[],tissue_audit_checkpoints=str(h/args.audit_checkpoints) if args.audit_checkpoints else str(bank/'checkpoints'),tissue_audit_canonical_frame='heart-native-local-metres',interface_checks=interfaces,assembly=dict(whole_organ=marker['wholeConceptId'],part_ids=next(c['elements']for c in a['concepts']if c['id']==marker['wholeConceptId']),fixed=['skeleton','lungs','diaphragm','original source controls and overall cardiac placement','unmodified vessel courses'],adjoining_interfaces=semantic,remaining_limits=remaining),evidence_files=[str(bank/'emission.json'),str(bank/'checked-geometry.json'),str(bank/'storage-preflight.json'),str(source/'ownership.json'),str(out/'transfer.json')]);atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec));print(sex,len(new_ids),'finite cells;',len(interfaces),'reciprocal tissue charts',flush=True)
