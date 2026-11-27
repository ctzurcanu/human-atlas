from pathlib import Path
import copy
import hashlib
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.db import Database
from geometry_tools.model import AtlasReader
from geometry_tools.registry import inventory

base=ROOT/'.local-models/geometry-workbench/liver-completion'
output=base/'male-liver-section-install-v142'
output.mkdir(exist_ok=False)
path=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
raw=path.read_bytes();before=json.loads(raw);after=copy.deepcopy(before)
preview=json.loads((base/'liver-section-volumes-v141/candidate.json').read_bytes())
assert hashlib.sha256(raw).hexdigest()==preview['liverSectionIntegration']['sourceCatalogueSha256']
old_chunk=next(p['chunk']for p in preview['parts']if p['id']=='TA98PART:ZA:Liver:right')
new_chunk=len(after['chunks']);after['chunks'].append(preview['chunks'][old_chunk])
changed=before['liverIntegration']['partIds']
for p in after['parts']:
    if p['id'] in changed:p['chunk']=new_chunk
after['liverIntegration']=preview['liverIntegration']
after['liverSectionIntegration']=preview['liverSectionIntegration']
reader=AtlasReader();old_parts={p['id']:p for p in before['parts']};literal={}
for p in after['parts']:
    if p['id']not in changed:continue
    old=reader.read_part(before,old_parts[p['id']],frame='model:local-ta98:assembled-metres')
    new=reader.read_part(after,p,frame=old.frame)
    assert old.geometry_sha==new.geometry_sha and old.attributes_sha==new.attributes_sha
    literal[p['id']]=dict(geometry_sha256=old.geometry_sha,attributes_sha256=old.attributes_sha)
old_inventory=json.loads((base/'current-inventory-v88.json').read_text())
classes={t['authority']['id']:t['representation_class']for t in old_inventory['terms']}
for n in [*range(2,11),12,13,14,*range(18,24)]:classes[f'A05.8.01.{n:03d}']='surface-region'
current_inventory=inventory(Database(),'A05.8.01.001','local-ta98',classifications=classes)
atomic_json(output/'inventory.json',current_inventory)
atomic_json(output/'atlas.json',after)
atomic_json(output/'literal-display-evidence.json',dict(passed=True,parts=literal,
    sourceCatalogueSha256=hashlib.sha256(raw).hexdigest(),otherPartsLiteral=True,
    skeletonLiteral=True,addedTa98Terms=0))
spec=dict(model='local-ta98',proposed_catalogue=str(output/'atlas.json'),
    expected_catalogue_sha256=hashlib.sha256(raw).hexdigest(),changed_parts=changed,
    inventory=str(output/'inventory.json'),require_closed=False,
    closed_parts=['TA98PART:ZA:Liver:capsule'],max_pairs=2000000,max_tests=16000000,
    assembly=dict(whole_organ='ZA:Liver',part_ids=changed,
        adjoining_interfaces=[dict(role='parenchymal-section',source='Four audited closed source volumes minus the union of five audited conduit envelopes. Computational disks classify pixels only; none are displayed.'),
            dict(role='vascular-and-biliary-hila',source='Current external surfaces, actual vessel walls and all adjoining interfaces remain byte-identical.')],
        fixed=['complete skeleton','diaphragm','all five literal liver tissue geometries','all vascular and biliary geometries','all unrelated anatomy']),
    evidence_files=[str(base/'liver-section-volumes-v141/evidence.json'),str(base/'liver-section-volumes-v141/checks.json'),str(output/'literal-display-evidence.json')])
atomic_json(output/'install-spec.json',spec)
atomic_json(output/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(candidate=str(output/'atlas.json'),terms=len(current_inventory['terms']),changed_parts=changed)))
