"""Prepare a guarded existing-bank heart pose update after packing review."""
from pathlib import Path
import argparse,copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH
ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True);ap.add_argument('--stage',required=True);ap.add_argument('--context-stages',nargs='*',default=[]);ap.add_argument('--volume-stage');ap.add_argument('--chest-stage');ap.add_argument('--bronchi-stage');ap.add_argument('--lungs-stage');args=ap.parse_args()
h=ROOT/'.local-models/geometry-workbench/heart';candidate=h/args.candidate;stage=h/args.stage;stage.mkdir(exist_ok=False)
current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json';a=json.loads((candidate/'atlas.json').read_text());before=json.loads(current.read_text());fit=json.loads((candidate/'packing.json').read_text());assert fit['source_catalogue_sha256']==digest_file(current)
moving=fit['moving_parts'];parts={p['id']:p for p in a['parts']};bp={p['id']:p for p in before['parts']};reader=AtlasReader();frame='model:local-female-ta98:assembled-metres'
for pid in moving:
 if pid in fit.get('additional_modified_parts',[]):
  trim=parts[pid].get('cardiacRootReferenceTrim',{});assert trim.get('capsAdded')==0 and trim.get('finiteWallEstablished')is False and trim.get('sourceCatalogueSha256')==fit['source_pose_candidate_sha256'];continue
 old,new=copy.deepcopy(bp[pid]),copy.deepcopy(parts[pid])
 for p in [old,new]:
  for key in ['bounds','geometryOrigin','geometryInstance','cardiacThoracicFit']:p.pop(key,None)
 assert old==new,('An instance refit changed a source record',pid)
assert a['chunks'][:len(before['chunks'])]==before['chunks']
if args.volume_stage:
 assert all([args.chest_stage,args.bronchi_stage,args.lungs_stage]);sha=digest_file(candidate/'atlas.json');paths=[h/name/'report.json'for name in [args.volume_stage,args.chest_stage,args.bronchi_stage,args.lungs_stage]];volume,chest,bronchi,lungs=[json.loads(p.read_text())for p in paths]
 assert volume['source_catalogue_sha256']==sha and volume['no_heart_lung_tissue_overlap']
 for r in [chest,bronchi,lungs]:assert r['catalogue_sha256']==sha and r['complete_within_pair_budget']
 for r in [chest,bronchi]:assert r['crossing_face_pairs']==0 and r['coplanar_overlap_face_pairs']==0
 assert all(row['heart']in a['finiteCardiacIntegration']['rootPartIds']and not parts[row['heart']].get('pairedTissue')for row in lungs['part_pairs'])
 remaining='Female heart body has zero measured occupied-tissue overlap with either lung, zero wall/bronchus crossings, and no full-assembly crossings with the native sternum, ribs, chest-wall muscles or diaphragm. Separate aortic-arch and pulmonary-trunk reference geometry retains '+str(lungs['crossing_face_pairs'])+' lung-surface crossing pairs. Exact distal vessel joins, including the retained open caval cuff, fine anatomy and clinical size remain unverified.'
 extra_evidence=[str(p)for p in paths]
else:
 context=[json.loads((h/name/'context.json').read_text())for name in args.context_stages]
 assert {r['context']for r in context}=={'lungs','diaphragm','bones','oesophagus'}
 for r in context:
  assert r['source_catalogue_sha256']==digest_file(current)and r['candidate_catalogue_sha256']==digest_file(candidate/'atlas.json')
  assert all(p['complete_within_pair_budget']for p in r['contacts'])
  if r['context']in ['bones','diaphragm']:assert r['candidate_no_crossings']
  assert r['contacts'][1]['crossing_face_pairs']<=r['contacts'][0]['crossing_face_pairs'],r['context']
 remaining='Female common-instance thoracic correction; sternum, native rib cage and diaphragm clear. Residual pulmonary packing contacts, fine anatomy, arterial entry repair and distal recipient vessel joins remain unfinished. No clinical size or complete anatomical acceptance is claimed.'
 extra_evidence=[str(h/name/'context.json')for name in args.context_stages]
spec=copy.deepcopy(json.loads((h/'annular-female-install-v384/install-spec.json').read_text()));spec.update(proposed_catalogue=str(candidate/'atlas.json'),expected_catalogue_sha256=digest_file(current),changed_parts=sorted(moving),closed_parts=sorted(pid for pid in moving if parts[pid].get('pairedTissue')and parts[pid].get('sectionSolid')))
inventory=json.loads((h/'annular-female-install-v384/inventory.json').read_text());inventory['source_catalogue_sha256']=digest_file(current);atomic_json(stage/'inventory.json',inventory);spec['inventory']=str(stage/'inventory.json')
spec['assembly']['remaining_limits']=remaining;spec['assembly']['fixed']=['every nonmember remains literal','skeleton','lungs','diaphragm','oesophagus','native complete aortic arch and its branches']
spec['evidence_files']=[p for p in spec['evidence_files']if not p.endswith('/transfer.json')and not p.endswith('/surface-features.json')]
forest=candidate/'source-forest.json'
if not forest.exists():forest=h/'female-thoracic-clear-candidate-v456/source-forest.json'
spec['evidence_files'] += [str(candidate/'packing.json'),str(forest)]+extra_evidence
if (candidate/'trim.json').exists():spec['evidence_files'].append(str(candidate/'trim.json'))
heart=reader.read_part(a,parts['TA98HEART:female:epicardium'],frame=frame,catalogue_path=candidate/'atlas.json');sternum=reader.read_part(a,parts['LOCAL:female:sternum:342'],frame=frame,catalogue_path=candidate/'atlas.json');nearest=BVH(heart).closest_mesh(BVH(sternum));assert nearest['distance_m']>.0005
spec['interface_checks'].append(dict(kind='clearance',first='TA98HEART:female:epicardium',second='LOCAL:female:sternum:342',parameters=dict(minimum_m=.0005)))
atomic_json(stage/'sternal-minimum-distance.json',nearest);spec['evidence_files'].append(str(stage/'sternal-minimum-distance.json'))
atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec))
print('Prepared',len(moving),'common-instance part updates; actual minimum sternum gap',nearest['distance_m']*1000,'mm',flush=True)
