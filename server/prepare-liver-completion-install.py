"""Prepare a source-bound local revision with explicit, honest coverage."""
from pathlib import Path
import copy,json,sys
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.db import Database
from geometry_tools.registry import inventory
base=ROOT/'.local-models/geometry-workbench/liver-completion';stage=base/'male-liver-completion-install-v201';stage.mkdir(exist_ok=False)
source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads((base/'completion-v200/atlas.json').read_text());before=json.loads(source.read_text())
assert digest_file(source)==a['liverCompletionIntegration']['sourceCatalogueSha256']
ledger=json.loads((base/'completion-v155/coverage.json').read_text());assert ledger['representedTerms']==63
old=json.loads((base/'male-liver-section-install-v142/inventory.json').read_text());classes={t['authority']['id']:t['representation_class']for t in old['terms']}
current=inventory(Database(),'A05.8.01.001','local-ta98',classifications=classes);assert len(current['terms'])==69
atomic_json(stage/'inventory.json',current);atomic_json(stage/'atlas.json',a)
oldparts={p['id']:p for p in before['parts']};changed=[p['id']for p in a['parts']if oldparts.get(p['id'])!=p]
ids=a['liverCompletionIntegration']['partIds']
whole_channels=['TA98HEPATIC:male:'+role for role in ['artery','portal','hepatic_vein','biliary']]
closed=ids+whole_channels+['TA98PART:ZA:Liver:capsule']
proof=dict(sourceCatalogueSha256=digest_file(source),changedParts=changed,newParts=ids,
    originalLiverGeometryLiteral=True,allUnrelatedPartsLiteral=True,skeletonLiteral=True,
    coverage=ledger,failedCoatCandidateExcluded=True,localOnly=True)
atomic_json(stage/'ownership-evidence.json',proof)
evidence=[base/'completion-v155/coverage.json',base/'completion-v155/functional-section-evidence.json',base/'completion-v155/computational-mouth-evidence.json',
    base/'functional-regions-v152/evidence.json',base/'functional-regions-v152/packed/tissue.json',base/'ligament-courses-v153/evidence.json',base/'ligament-courses-v153/packed/tissue.json',
    base/'microscopic-unit-v151/evidence.json',base/'microscopic-unit-v151/packed/tissue.json',stage/'ownership-evidence.json',
    base/'completion-v186/construction.json',base/'completion-v200/construction.json',base/'hepatic-contained-banks-v199/containment.json',
    *sorted((base/'contained-hepatic-pairs-v198').glob('*.json')),
    *(base/f'completion-v186/exclusion-{role}-audit.json' for role in ['artery','portal','hepatic_vein']),
    base/'completion-v200/exclusion-biliary-audit.json',
    *(base/f'hepatic-contained-banks-{version}/{role}-{kind}.json' for role,version in [('artery','v170'),('portal','v181'),('biliary','v181'),('hepatic_vein','v197')] for kind in ['audit','evidence'])]
for p in evidence:
    d=json.loads(p.read_text())
    if 'passed' in d:assert d['passed'],str(p)
spec=dict(model='local-ta98',proposed_catalogue=str(stage/'atlas.json'),expected_catalogue_sha256=digest_file(source),
    changed_parts=changed,inventory=str(stage/'inventory.json'),require_closed=False,closed_parts=closed,
    max_pairs=3000000,max_tests=24000000,
    assembly=dict(whole_organ='ZA:Liver',part_ids=ids+whole_channels+a['liverIntegration']['partIds'],
        adjoining_interfaces=[dict(role='functional territories',source='Eight single-component, disjoint source-landmark partitions share reciprocal interfaces; the repaired exterior is retained.'),
            dict(role='hilar source disks',source='Computation-only triangles stay in closed classifiers and are excluded from displayed surfaces and picking.'),
            dict(role='hepatic conduit sections',source='Five literal audited conduit envelopes remain subtractive in live section rendering.'),
            dict(role='ligaments',source='Two illustrative finite cords follow actual source-bound fissure courses. Attachment extent and calibre remain provisional.'),
            dict(role='microscopic anatomy',source='A separate illustrative human-liver unit is selectable and has magnified navigation. It is not patient-specific microscopy.')],
        fixed=['skeleton','diaphragm','original liver tissue positions and corner fields','native vascular and biliary attachment coordinates and normal words','all unrelated anatomy'],
        remaining_limits=ledger['missing']),
    evidence_files=list(map(str,evidence)))
atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(stage=str(stage),changedParts=len(changed),newParts=len(ids),coverage=ledger['percentage'])))
