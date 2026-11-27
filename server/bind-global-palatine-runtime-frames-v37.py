"""Bind the extra v33 domain sources to actual loader participation and exact decoded geometry."""
from pathlib import Path
import json,sys,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from collections import Counter
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.db import Database

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17/global-palatine-graph-v33';P=D/'additional-runtime-v36.json';runtime=json.loads(P.read_bytes());expected=runtime['source_catalogue_sha256'];revision=runtime['registry_revision'];source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();ids=[r['id']for r in runtime['rows']];assert len(ids)==117
    db=Database();db.assert_current_catalogues('local-ta98')
    with db.connect()as connection:
        current=connection.execute('SELECT a.catalogue_sha256,c.metadata FROM atlas_geom.current_catalogues a JOIN atlas_geom.catalogues c ON c.sha256=a.catalogue_sha256 WHERE a.model_id=%s',('local-ta98',)).fetchone();assert current['catalogue_sha256']==revision and current['metadata']['json_sha256']==expected
        records=connection.execute('SELECT id,part_id,frame_id,geometry_sha256,attributes_sha256,metadata,source,error FROM atlas_geom.meshes WHERE catalogue_sha256=%s AND part_id=ANY(%s)',(revision,ids)).fetchall()
    registered={r['part_id']:r for r in records};assert set(registered)==set(ids)
    frame='model:local-ta98:assembled-metres';height=atlas['heightAdjustment'];previous_path=reader.asset(height['previousCatalogue']);assert digest_file(previous_path)==height['previousAtlasSha256'];before=json.loads(previous_path.read_bytes());before_parts={p['id']:p for p in before['parts']};rows=[];bindings={str(source):expected,str(P):digest_file(P),str(previous_path):height['previousAtlasSha256']}
    for state in runtime['rows']:
        if time.monotonic()-began>47:raise RuntimeError('Frame decode internal deadline reached; incomplete evidence cannot be used')
        pid=state['id'];p=parts[pid];r=registered[pid];assert not r['error'] and r['metadata']['record']==p and r['frame_id']==frame
        mesh=reader.read_part(atlas,p,frame=frame,catalogue_path=source,catalogue_sha=expected);assert mesh.geometry_sha==r['geometry_sha256'] and mesh.attributes_sha==r['attributes_sha256'];assert mesh.source['bank_sha256']==r['source']['bank_sha256'];bindings[mesh.source['bank_path']]=mesh.source['bank_sha256']
        if state['drawing_kind']=='instanced-metric-hair-shafts-not-carrier-triangles':kind='instanced-field-draw';reason='Actual scene drawing chooses the explicit hairField prototype and roots.'
        elif state.get('tissueOwnershipClaim')is False or 'reference' in (state.get('representationMode')or''):
            kind='explicit-nonowning-anatomical-reference-envelope';reason='The record explicitly declares non-owning reference representation or denies tissue ownership. It is not an occupied-tissue classifier.'
        elif state['finite_tissue_metadata'] and state.get('tissueDrawGroups'):
            kind='owned-derivative-finite-tissue-cell';reason='Explicit pairedTissue, sectionSolid and tissueDrawGroups metadata; this identifies intended cell ownership but does not independently accept anatomy.'
        elif state['system']=='attachments':
            kind='declared-attachment-reference-surface';reason='Explicit attachments system and source material record describe attachment drawings. Finite tissue volume is not established by those triangles.'
        elif state['declared_point_or_line_reference']:
            kind='declared-point-or-curve-reference';reason='Explicit point/line/curve metadata; never infer solid occupancy from a locator.'
        else:kind='physical-or-anatomical-indexed-surface-source';reason='Actual runtime draws the indexed triangles; no reference-envelope/field override or finite-owned-cell declaration establishes another class. Volume and anatomical registration remain to be checked.'
        row=dict(state,representation_class=kind,representation_class_basis=reason,retained_as_possible_obstacle=not state['normal_runtime_suppressed'],
          registry_geometry_binding=dict(registry_mesh_id=r['id'],registry_revision=revision,record_matches_registry=True,geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,bank_path=mesh.source['bank_path'],bank_sha256=mesh.source['bank_sha256'],world_frame=frame,coordinate_rule=mesh.source['coordinates'],world_bounds_m=mesh.bounds.tolist(),triangles=len(mesh.faces),decoded_actual_record_once=True),
          independent_anatomical_registration_not_established=True,provider_prefix_not_used_to_exclude=True)
        if pid in before_parts:
            old=reader.read_part(before,before_parts[pid],frame='historical-pre-height-source-metres',catalogue_path=previous_path,catalogue_sha=height['previousAtlasSha256']);same=old.vertices.shape==mesh.vertices.shape and np.array_equal(old.faces,mesh.faces)and old.source['source_vertex_ids']==mesh.source['source_vertex_ids'];error=float(np.max(np.linalg.norm(old.vertices*height['scale']+np.array(height['translationMeters'])-mesh.vertices,axis=1)))if same else None;tolerance=float(height['maximumPositionRoundingMeters'])*np.sqrt(3)+1e-10
            row['saved_height_transform_replay']=dict(verified=bool(same and error<=tolerance),same_source_controls_and_faces=same,scale=height['scale'],translation_m=height['translationMeters'],maximum_vertex_residual_m=error,tolerance_m=tolerance,source_geometry_sha256=old.geometry_sha,anatomical_landmark_registration_not_inferred=True)
        else:row['saved_height_transform_replay']=dict(verified=False,reason='No literal same-ID pre-height source; explicit current pose and bank remain verified.')
        rows.append(row)
    retained=[r for r in rows if not r['normal_runtime_suppressed']];visible0=[r for r in retained if r['normal_default_visible_skin0']];visible1=[r for r in retained if r['normal_visible_skin1']]
    result=dict(schema='global-palatine-additional-runtime-frame-binding/v37',source_catalogue_sha256=expected,registry_revision=revision,ordinary_runtime_input_sha256=runtime['ordinary_body_sha256'],runtime_evidence=str(P.relative_to(ROOT)),runtime_evidence_sha256=digest_file(P),code_bindings=runtime['code_bindings'],file_bindings=bindings,source_bounds_only_no_contact_inference=True,
      counts=dict(total=len(rows),suppressed_by_actual_default_loader=len(rows)-len(retained),retained=len(retained),visible_skin0=len(visible0),visible_skin1=len(visible1),retained_classes=dict(Counter(r['representation_class']for r in retained)),suppression_basis=dict(Counter(r['exclusion_basis']for r in rows if r['normal_runtime_suppressed'])),height_transform_replays=sum(r['saved_height_transform_replay']['verified']for r in rows)),
      retained_runtime_ids=[r['id']for r in retained],retained_indexed_surface_or_tissue_ids=[r['id']for r in retained if r['representation_class']in['physical-or-anatomical-indexed-surface-source','owned-derivative-finite-tissue-cell']],rows=rows,
      no_geometry_visibility_or_catalogue_mutation=True,no_source_prefix_exclusions=True,no_solid_occupancy_or_anatomical_acceptance=True,elapsed_seconds=time.monotonic()-began)
    target=D/'additional-participation-frame-v37.json';atomic_json(target,result);assert digest_file(source)==expected;print(json.dumps(dict(output=str(target.relative_to(ROOT)),counts=result['counts'],elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
