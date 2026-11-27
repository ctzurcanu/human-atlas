"""Prepare a guarded, rollback-capable local liver install from passed evidence."""
from pathlib import Path
import argparse,copy,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--sex',choices=['male','female'],required=True);ap.add_argument('--candidate');ap.add_argument('--stage');args=ap.parse_args();sex=args.sex;model='local-female-ta98'if sex=='female'else'local-ta98';b=ROOT/'.local-models/geometry-workbench/liver-completion';candidate=Path(args.candidate)if args.candidate else b/'final-ready-v265'/(sex+'.json');a=json.loads(candidate.read_text());source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json';before=json.loads(source.read_text());bp={p['id']:p for p in before['parts']};stage=b/(args.stage or'final-install-v266-'+sex);stage.mkdir(exist_ok=False)
assert json.loads((candidate.parent/(sex+'-context.json')).read_text())['passed'];assert json.loads((b/'final-ready-v265'/(sex+'-interfaces.json')).read_text())['passed']
geometry_keys=['positions','normals','indices','normalSigns','uvs','colors','gastricLumen','tissueChart','vertexCount','indexCount','triangles','geometryOrigin','geometryScale','sourceOffset','chunk','bounds','tissueDrawGroups']
for p in a['parts']:
 if p.get('suppressed')and p['id']in bp:
  for key in geometry_keys:
   if key in bp[p['id']]:p[key]=copy.deepcopy(bp[p['id']][key])
   else:p.pop(key,None)
proposed=stage/'candidate.json';atomic_json(proposed,a);changed=[p['id']for p in a['parts']if p!=bp.get(p['id'])];closed=[p['id']for p in a['parts']if p['id']in changed and not p.get('suppressed')and p.get('pairedTissue')and p.get('sectionSolid')]
inv=copy.deepcopy(json.loads((b/'final-repair-plan-v209/plan.json').read_text())['authority_inventory']);inv.update(model=model,model_sex=sex,source_catalogue_sha256=digest_file(source));groups={g['term']:g for key in ['liverIntegration','liverFeatureIntegration','liverCompletionIntegration','hepaticConnectionIntegration']for g in a[key]['termGroups']}
for t in inv['terms']:
 g=groups.get(t['authority']['id']);t['candidate_part_ids']=g.get('parts',g.get('partIds',[]))if g else[];t['candidate_representation']=g.get('representation','finite-tissue-or-owned-region')if g else None;t['regional_boundary_provisional']=True
inventory=stage/'inventory.json';atomic_json(inventory,inv)
interfaces=[dict(kind='shared-interface',first=q['first'],second=q['second'],parameters=dict(first_faces=q['first_faces'],second_faces=q['second_faces']))for q in a['liverFinalIntegration']['interfaces']]
evidence=[candidate.parent/(sex+'-context.json'),b/'final-ready-v265'/(sex+'-interfaces.json'),b/'final-processes-v223/processes.json',b/'right-caudate-duct-v226/network.json',b/'right-caudate-duct-v226/adjoining-checks.json']
if(candidate.parent/'field-evidence.json').exists():evidence.append(candidate.parent/'field-evidence.json')
spec=dict(model=model,proposed_catalogue=str(proposed),expected_catalogue_sha256=digest_file(source),changed_parts=changed,inventory=str(inventory),require_closed=False,closed_parts=closed,max_pairs=1000000,max_tests=12000000,interface_checks=interfaces,assembly=dict(whole_organ=a['liverIntegration']['wholeConceptId'],part_ids=a['liverIntegration']['partIds'],fixed=['skeleton','diaphragm','oesophagus','gastric mouths'],adjoining_interfaces=[dict(kind='disjoint',scope='stomach, fixed oesophagus, diaphragm and active skeleton',evidence=str(candidate.parent/(sex+'-context.json'))),dict(kind='shared-tissue-interface',scope='three independently matched annular walls',evidence=str(b/'final-ready-v265'/(sex+'-interfaces.json')))],remaining_limits=a['liverFinalIntegration']['remainingConnectionLimits']),evidence_files=list(map(str,evidence)))
atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec));print(sex,len(changed),'changed parts;',len(closed),'finite cells')
