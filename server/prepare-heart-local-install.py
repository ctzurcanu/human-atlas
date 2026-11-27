"""Bind the reviewed cardiac bank, inventory and reciprocal seams for local install."""
from pathlib import Path
import argparse,copy,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],default='male');ap.add_argument('--candidate',required=True);ap.add_argument('--wall',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();b=ROOT/'.local-models/geometry-workbench/heart';candidate=b/args.candidate/'atlas.json';wall=b/args.wall;checks=json.loads((wall/'checked-geometry-v3.json').read_text());assert checks['passed'];a=json.loads(candidate.read_text());sex=args.sex;model='local-ta98'if sex=='male'else'local-female-ta98';source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';before=json.loads(source.read_text());bp={p['id']:p for p in before['parts']};marker=a['finiteCardiacIntegration'];out=b/args.stage;out.mkdir(exist_ok=False);inventory=copy.deepcopy(json.loads((b/('source-v1-'+sex)/'inventory.json').read_text()));inventory['source_catalogue_sha256']=digest_file(source);byterm={}
for group in marker['termGroups']:
 if group['kind']=='exact':byterm.setdefault(group['term'],[]).extend(group['partIds'])
for row in inventory['terms']:
 row['candidate_part_ids']=list(dict.fromkeys(byterm.get(row['authority']['id'],[])));row['candidate_representation']='stored finite tissue, source valve tissue or collective control'if row['candidate_part_ids']else None;row['regional_boundary_provisional']=True
atomic_json(out/'inventory.json',inventory);e=json.loads((wall/'emission.json').read_text());layout=dict(e,parts={e['part_ids'][name]:record for name,record in e['parts'].items()});atomic_json(out/'audited-layout.json',layout);interfaces=[]
for pid,record in layout['parts'].items():
 for group in record['tissueDrawGroups']:
  neighbor=group.get('neighbor')
  if not neighbor or pid>=neighbor:continue
  peer=layout['parts'][neighbor];other=[g for g in peer['tissueDrawGroups']if g.get('neighbor')==pid];interfaces.append(dict(kind='shared-interface',first=pid,second=neighbor,parameters=dict(first_faces=list(range(group['start']//3,(group['start']+group['count'])//3)),second_faces=[i for g in other for i in range(g['start']//3,(g['start']+g['count'])//3)])))
changed=[p['id']for p in a['parts']if p!=bp.get(p['id'])];closed=[p['id']for p in a['parts']if p['id']in changed and not p.get('suppressed')and p.get('sectionSolid')];spec=dict(model=model,proposed_catalogue=str(candidate),expected_catalogue_sha256=digest_file(source),changed_parts=changed,closed_parts=closed,require_closed=False,max_pairs=1300000,max_tests=14000000,inventory=str(out/'inventory.json'),tissue_audit_layouts=[] if any(p.get('geometryInstance')for p in a['parts']if p['id']in marker['wallPartIds']) else [str(out/'audited-layout.json')],tissue_audit_checkpoints=str(wall/'checkpoints'),interface_checks=interfaces,assembly=dict(whole_organ=marker['wholeConceptId'],part_ids=a['concepts'][next(i for i,c in enumerate(a['concepts'])if c['id']==marker['wholeConceptId'])]['elements'],fixed=['skeleton','lungs','diaphragm'],adjoining_interfaces=[dict(kind='shared-interface',scope='all reciprocal cardiac tissue charts',count=len(interfaces)),dict(kind='open-conduit',scope='native valvular outflows, caval inflows and pulmonary venous openings',source=str(b/'native-pv-mouths-v69/construction.json'))],remaining_limits=marker.get('remainingLimits','')),evidence_files=[str(wall/'checked-geometry-v3.json'),str(b/'native-pv-mouths-v69/construction.json'),str(b/'anterior-valves-v52/construction.json')]);
if (wall/'transfer.json').exists():spec['evidence_files'].append(str(wall/'transfer.json'))
if marker.get('rootWallPartIds'):
 flow=wall/'root-fluid-paths.json'
 if not flow.exists() and e.get('bank'):flow=Path(e['bank']).parent/'root-fluid-paths.json'
 proof=json.loads(flow.read_text());assert proof['passed'] and proof['bank_sha256']==e['bank_sha256']
 interface=Path(e['root_interface_evidence']);assert digest_file(interface)==e['root_interface_evidence_sha256']
 spec['evidence_files'].extend([str(flow),str(interface),str(b/e['source_root_stage']/'construction.json')])
 spec['assembly']['adjoining_interfaces'].append(dict(kind='shared-interface',scope='four finite pulmonary-vein roots joined reciprocally to atrial tissue; eight complete fluid centerlines checked',source=str(flow)))
atomic_json(out/'install-spec.json',spec);atomic_json(out/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec));print(sex,len(changed),'changed parts;',len(closed),'closed tissues;',len(interfaces),'reciprocal interfaces',flush=True)
