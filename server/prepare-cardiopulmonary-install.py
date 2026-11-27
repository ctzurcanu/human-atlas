"""Bind the reviewed male heart/lung packing closure to guarded installation."""
from pathlib import Path
import argparse,copy,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH
ap=argparse.ArgumentParser();ap.add_argument('--candidate',required=True);ap.add_argument('--stage',required=True);args=ap.parse_args();h=ROOT/'.local-models/geometry-workbench/heart';candidate=h/args.candidate;stage=h/args.stage;stage.mkdir(exist_ok=False);path=candidate/'atlas.json';current=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';a=json.loads(path.read_text());before=json.loads(current.read_text());sha=digest_file(path);parts={p['id']:p for p in a['parts']};bp={p['id']:p for p in before['parts']};changed=sorted(i for i,p in parts.items()if p!=bp[i]);assert set(parts)==set(bp)
heart=json.loads((h/'male-cardiopulmonary-final-candidate-v500/packing.json').read_text());pulmonary=json.loads((h/'male-pulmonary-raised-cage-v494/packing.json').read_text());assert pulmonary['source_catalogue_sha256']==digest_file(current);assert set(changed)==set(heart['moving_parts'])|set(pulmonary['moving_parts']);assert len(changed)==84;assert a['chunks'][:len(before['chunks'])]==before['chunks']
names=['male-final-heart-lung-volume-v520','male-final-cardiac-chest-v521','male-final-cardiac-bronchi-v522','male-final-cardiac-roots-v523','male-final-pulmonary-chest-v524','male-pulmonary-cell-quality-v518','male-pulmonary-owned-check-v519'];reports=[json.loads((h/n/'report.json').read_text())for n in names];volume,chest,bronchi,roots,lung_chest,quality,owned=reports
assert volume['source_catalogue_sha256']==sha and volume['no_heart_lung_tissue_overlap']
for r in reports[1:]:assert r['catalogue_sha256']==sha
for r in [chest,bronchi]:assert r['complete_within_pair_budget']and r['crossing_face_pairs']==r['coplanar_overlap_face_pairs']==0
assert quality['complete']and quality['passed']and len(quality['cells'])==19
assert owned['passed']and len(owned['occupied_volume_pairs'])==171 and len(owned['named_destinations'])==20
assert roots['complete_within_pair_budget']and all(r['heart']=='ZA:Pulmonary trunk'for r in roots['part_pairs'])and not parts['ZA:Pulmonary trunk'].get('pairedTissue')
baseline=json.loads((h/'male-current-pulmonary-chest-v465/report.json').read_text());old={(r['heart'],r['context']):r['crossings']for r in baseline['part_pairs']};assert lung_chest['complete_within_pair_budget']and all(r['crossings']<=old.get((r['heart'],r['context']),0)for r in lung_chest['part_pairs'])
reader=AtlasReader();frame='model:local-male-ta98:assembled-metres'
def read(pid):return reader.read_part(a,parts[pid],frame=frame,catalogue_path=path)
for pid in pulmonary['moving_parts']:
 old_mesh=reader.read_part(before,bp[pid],frame=frame,catalogue_path=current);new=read(pid);assert np.array_equal(old_mesh.faces,new.faces);assert np.array_equal(old_mesh.attributes.get('uvs'),new.attributes.get('uvs'));assert np.array_equal(old_mesh.attributes.get('normalSigns'),new.attributes.get('normalSigns'))
spec=copy.deepcopy(json.loads((h/'annular-male-install-v382/install-spec.json').read_text()));spec.update(proposed_catalogue=str(path),expected_catalogue_sha256=digest_file(current),changed_parts=changed,closed_parts=sorted(pid for pid in changed if parts[pid].get('pairedTissue')and parts[pid].get('sectionSolid')))
inventory_paths=[h/'annular-male-install-v382/inventory.json',ROOT/'.local-models/geometry-workbench/pulmonary-completion/male-lung-install-v596/inventory.json'];inventories=[json.loads(p.read_text())for p in inventory_paths];inventory=copy.deepcopy(inventories[0]);terms={t['authority']['id']:t for inv in inventories for t in inv['terms']};inventory.update(terms=list(terms.values()),term_count=len(terms),root_terms=sorted(set(sum([inv.get('root_terms',inv.get('roots',[]))for inv in inventories],[]))),source_catalogue_sha256=digest_file(current),source_catalogue=dict(path=str(current),sha256=digest_file(current)),source_inventories=[dict(path=str(p),sha256=digest_file(p))for p in inventory_paths]);atomic_json(stage/'inventory.json',inventory);spec['inventory']=str(stage/'inventory.json')
lung_interfaces=0
for pid in a['pulmonarySegmentIntegration']['partIds']:
 for group in parts[pid]['tissueDrawGroups']:
  peer=group.get('neighbor')
  if not peer or pid>=peer:continue
  other=[g for g in parts[peer]['tissueDrawGroups']if g.get('neighbor')==pid];assert len(other)==1;other=other[0]
  params=dict(first_faces=list(range(group['start']//3,(group['start']+group['count'])//3)),second_faces=list(range(other['start']//3,(other['start']+other['count'])//3)))
  spec['interface_checks'].append(dict(kind='shared-interface',first=pid,second=peer,parameters=params));spec['assembly']['adjoining_interfaces'].append(dict(kind='literal-shared-pulmonary-chart',first=pid,second=peer,anatomical_boundary_provisional=True));lung_interfaces+=1
assert lung_interfaces==42
spec['assembly'].update(part_ids=sorted(set(spec['assembly']['part_ids'])|set(pulmonary['moving_parts'])),assembly_ids=['TA98HEART:male:whole','TA98PULMONARY:male:whole'],fixed=['every nonmember of the explicit heart/pulmonary closure','entire skeleton','diaphragm','chest-wall muscles','native vascular references','literal upper airway controls at and above Y 1.475 m'],remaining_limits=f"Heart-body/lung occupied overlap is zero bilaterally; the whole cardiac assembly clears the native sternum, ribs, chest muscles and diaphragm, and cardiac walls clear bronchi. All 19 closed pulmonary cells and 171 disjointness pairs pass; 20 named airway destinations persist. The unchanged coarse pulmonary-trunk reference retains {roots['crossing_face_pairs']} lung-surface crossing pairs, and the unchanged left apical source retains {lung_chest['crossing_face_pairs']} crossing pairs with left longus colli. Exact distal vessel joins, microscopic anatomy, finite airway walls, pulmonary territory anatomy and clinical size remain unverified. This is local packing improvement, not complete anatomical acceptance.")
sternal_ids=['ZA:Body of sternum:Bone-7','ZA:Manubrium of sternum:Bone-6','ZA:Xiphoid process:Bone-8'];rows=[];epi=BVH(read('TA98HEART:male:epicardium'))
for pid in sternal_ids:
 nearest=epi.closest_mesh(BVH(read(pid)));assert nearest['distance_m']>.0005;rows.append(dict(part_id=pid,**nearest));spec['interface_checks'].append(dict(kind='clearance',first='TA98HEART:male:epicardium',second=pid,parameters=dict(minimum_m=.0005)))
atomic_json(stage/'sternal-minimum-distances.json',rows)
forest=ROOT/'reports/structural-dependencies/current/male-forest.json';f=json.loads(forest.read_text());assert f['sourceSha256']==digest_file(current)
spec['evidence_files']=[p for p in spec['evidence_files']if not p.endswith('/transfer.json')and not p.endswith('/surface-features.json')]+[str(h/n/'report.json')for n in names]+[str(h/'male-pulmonary-raised-cage-v494/packing.json'),str(h/'male-cardiopulmonary-final-candidate-v500/packing.json'),str(candidate/'packing.json'),str(candidate/'repair.json'),str(forest),str(stage/'sternal-minimum-distances.json')]+[str(p)for p in inventory_paths]
atomic_json(stage/'install-spec.json',spec);atomic_json(stage/'job-spec.json',dict(kind='catalogue-candidate',parameters=spec));print('Prepared',len(changed),'parts,',len(spec['closed_parts']),'finite cells,',len(spec['interface_checks']),'interface/clearance checks; sternum min mm',min(r['distance_m']for r in rows)*1000,flush=True)
