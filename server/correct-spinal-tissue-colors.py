"""Correct RGB units only; all stored geometry and anatomical identities stay literal."""
from pathlib import Path
import argparse,copy,json
from geometry_tools.common import ROOT,atomic_json,digest_file
ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';before=json.loads(source.read_text());a=copy.deepcopy(before)
ids=a['finiteSpinalIntegration']['partIds'];colors={'grey':[158,122,115],'white':[232,214,185]}
for p in a['parts']:
 if p['id']not in ids:continue
 p['tissueAppearance']['color']=[v/255 for v in colors[p['tissueWallLayer']]]
for family,color in colors.items():
 m=a['materials']['spinal:female:'+family];m['color']=color;m['tissueAppearance']['color']=[v/255 for v in color]
for p,q in zip(before['parts'],a['parts']):
 old=copy.deepcopy(p);new=copy.deepcopy(q)
 if p['id']in ids:old['tissueAppearance'].pop('color');new['tissueAppearance'].pop('color')
 assert old==new
assert a['chunks']==before['chunks'] and a['concepts']==before['concepts'] and a['defaultHidden']==before['defaultHidden']
atomic_json(out/'atlas.json',a)
base=ROOT/'.local-models/geometry-workbench/spinal-cord-completion/female-install-v2'
spec=json.loads((base/'job.json').read_text())['parameters'];inventory=json.loads((base/'inventory.json').read_text())
inventory.update(source_catalogue=str(source),source_catalogue_sha256=digest_file(source));atomic_json(out/'inventory.json',inventory)
spec.update(proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(source),changed_parts=ids,closed_parts=ids,inventory=str(out/'inventory.json'))
spec['assembly']['remaining_limits'].append('This revision changes only RGB units; geometry and coverage are unchanged.')
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec));print(json.dumps(dict(changed='RGB units only',parts=len(ids),geometry_and_identifiers_literal=True)))
