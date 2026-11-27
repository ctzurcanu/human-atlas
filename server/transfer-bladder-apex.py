"""Transfer the existing male apex onto its literal female bladder carrier.

No mesh, texture, attachment, filling state or clinical validation is invented.
The saved core similarity supplies correspondence; the actual encoded target
facet supplies the selectable anatomical point.
"""
from pathlib import Path
import argparse,copy,json,time
import numpy as np
from geometry_tools.common import ROOT,atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.spatial import BVH

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();started=time.monotonic()
out=Path(args.output).resolve();out.mkdir(exist_ok=False)
base=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'
paths={s:base/s/'atlas.json' for s in ['male','female']}
atlases={s:json.loads(p.read_text()) for s,p in paths.items()}
male,female=atlases['male'],atlases['female'];before=copy.deepcopy(female)
term='A08.3.01.002';source=next(g for g in male['finiteUrinaryIntegration']['termGroups'] if g['term']==term)
ref=source['pointReference'];source_id=ref['partId'];target_id=source_id.replace(':male:',':female:')
records={s:{p['id']:p for p in a['parts']} for s,a in atlases.items()};reader=AtlasReader()
source_mesh=reader.read_part(male,records['male'][source_id],frame='model:local-ta98:assembled-metres',catalogue_path=paths['male'])
point=np.asarray(ref['point']);assert np.linalg.norm(np.asarray(ref['barycentric'])@source_mesh.triangles[ref['face']]-point)<1e-12
plan_path=ROOT/'.local-models/geometry-workbench/bladder-completion/female-transfer-v1/transfer-plan.json'
plan=json.loads(plan_path.read_text());matrix=np.asarray(plan['matrix']);scale=plan['uniform_scale']
assert np.allclose(matrix[:3,:3].T@matrix[:3,:3],np.eye(3)*scale**2,rtol=1e-12,atol=1e-14) and np.linalg.det(matrix[:3,:3])>0
mapped=matrix[:3,:3]@point+matrix[:3,3];target_part=records['female'][target_id]
target=reader.read_part(female,target_part,frame='model:local-female-ta98:assembled-metres',catalogue_path=paths['female'])
query,proof=BVH(target).closest_points(mapped[None,:]);face=int(query['face_ids'][0]);weights=query['barycentric'][0]
actual=weights@target.triangles[face];distance=float(np.linalg.norm(actual-mapped))
assert distance<1e-7,'The saved core similarity no longer establishes this correspondence'
assert any(not g.get('neighbor') and g['start']//3<=face<(g['start']+g['count'])//3 for g in target_part['tissueDrawGroups'])
marker=female['finiteUrinaryIntegration'];bank=next(g['bankSha256'] for g in marker['charts'] if target_id in g['partIds'])
target_ref=dict(representation='point-landmark',partId=target_id,point=actual.tolist(),face=face,barycentric=weights.tolist(),sourceBankSha256=bank,sourcePlanSha256=digest_file(plan_path))
transfer=dict(schema='bladder-apex-similarity/v1',sourcePartId=source_id,sourceCatalogueSha256=digest_file(paths['male']),sourceBankSha256=ref['sourceBankSha256'],sourcePoint=point.tolist(),sourcePartRecord=copy.deepcopy(records['male'][source_id]),sourceChunk=copy.deepcopy(male['chunks'][records['male'][source_id]['chunk']]),sourceReferenceFace=ref['face'],sourceBarycentric=ref['barycentric'],matrix=matrix.tolist(),transferPlanPath=str(plan_path),transferPlanSha256=digest_file(plan_path),maximumSourceOffsetM=distance,targetPartId=target_id)
concept_id='URINARY:female:bladder:'+term
assert not any(g['term']==term for g in marker['termGroups'])
assert not any(c['id']==concept_id for c in female['concepts'])
marker['apexTransfer']=transfer
marker['termGroups'].append(dict(term=term,conceptId=concept_id,partIds=[target_id],pointReference=target_ref))
female['concepts'].append(dict(id=concept_id,name='Apex of bladder',elements=[target_id],componentElements=[],ta98Term=term,ta98Kind='exact',partOf='URINARY:female:bladder:A08.3.01.001',coverageLimitation='Source-corresponding point on the transferred bladder exterior. The female median umbilical ligament and its tissue attachment remain unbuilt.'))
assert all(female[k]==before[k] for k in ['parts','materials','chunks','defaultHidden'])
inventory_path=ROOT/'.local-models/geometry-workbench/bladder-completion/female-transfer-v1/inventory.json'
inventory=copy.deepcopy(json.loads(inventory_path.read_text()));inventory.update(source_catalogue=str(paths['female']),source_catalogue_sha256=digest_file(paths['female']),source_male_catalogue=str(paths['male']),source_male_catalogue_sha256=digest_file(paths['male']))
for row in inventory['terms']:
 code=row['authority']['id'];row['source_meshes']=[dict(part_id=p['id'],suppressed=bool(p.get('suppressed')),chunk=p['chunk']) for p in female['parts'] if p.get('ta98Term')==code]
atomic_json(out/'inventory.json',inventory);atomic_json(out/'atlas.json',female)
evidence=dict(schema='transferred-bladder-point/v1',source_catalogue_sha256=digest_file(paths['male']),target_catalogue_sha256=digest_file(paths['female']),candidate_sha256=digest_file(out/'atlas.json'),source_reference=ref,transfer=transfer,target_reference=target_ref,source_facet=source_mesh.triangles[ref['face']].tolist(),target_facet=target.triangles[face].tolist(),source_query=proof,maximum_source_offset_m=distance,parts_materials_chunks_and_visibility_literal=True,new_geometry_bytes=0,new_texture_bytes=0,anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-started)
atomic_json(out/'reference.json',evidence)
spec=dict(model='local-female-ta98',proposed_catalogue=str(out/'atlas.json'),expected_catalogue_sha256=digest_file(paths['female']),changed_parts=[],inventory=str(out/'inventory.json'),assembly=dict(whole_organ='Female bladder with male-derived source-corresponding apex point',adjoining_interfaces='Unchanged bilateral ureteric and female urethral finite walls; all bones and pelvic tissues remain literal. No ligament is implied by this point.'),evidence=dict(reference=str(out/'reference.json')))
atomic_json(out/'job.json',dict(kind='catalogue-candidate',parameters=spec))
print(json.dumps(dict(candidate=str(out/'atlas.json'),source_offset_nm=distance*1e9,new_geometry_bytes=0,new_texture_bytes=0,elapsed_seconds=time.monotonic()-started)),flush=True)
