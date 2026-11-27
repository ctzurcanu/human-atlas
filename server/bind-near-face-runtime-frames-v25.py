"""Bind runtime participation to literal registered source coordinates only."""
from pathlib import Path
import json
import sys

ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json,digest_file
from geometry_tools.model import AtlasReader
from geometry_tools.db import Database

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17'
SOURCE=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
REVISION='31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2'


def main():
    assert digest_file(SOURCE)==EXPECTED
    atlas=json.loads(SOURCE.read_bytes());parts={p['id']:p for p in atlas['parts']}
    runtime_path=BASE/'nearby-unresolved-runtime-v25.json';runtime=json.loads(runtime_path.read_bytes());assert runtime['source_catalogue_sha256']==EXPECTED
    old_context_path=BASE/'nearby-derivative-context-inventory-v24.json';old_context=json.loads(old_context_path.read_bytes());registered={r['part_id']:r for r in old_context['all_registered_overlap_candidates']}
    db=Database();db.assert_current_catalogues('local-ta98')
    with db.connect() as c:
        current=c.execute('SELECT a.catalogue_sha256,c.metadata FROM atlas_geom.current_catalogues a JOIN atlas_geom.catalogues c ON c.sha256=a.catalogue_sha256 WHERE a.model_id=%s',('local-ta98',)).fetchone()
        assert current['catalogue_sha256']==REVISION and current['metadata']['json_sha256']==EXPECTED
    reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    height=atlas['heightAdjustment'];before_path=reader.asset(height['previousCatalogue']);before=None
    if before_path.exists():
        assert digest_file(before_path)==height['previousAtlasSha256'];before=json.loads(before_path.read_bytes());before_parts={p['id']:p for p in before['parts']}
    rows=[];bindings={str(SOURCE):EXPECTED,str(runtime_path):digest_file(runtime_path),str(old_context_path):digest_file(old_context_path)}
    if before is not None:bindings[str(before_path)]=digest_file(before_path)
    for status in runtime['rows']:
        pid=status['id'];p=parts[pid];reference=registered[pid]
        mesh=reader.read_part(atlas,p,frame=frame,catalogue_path=SOURCE,catalogue_sha=EXPECTED)
        assert mesh.geometry_sha==reference['geometry_sha256'] and mesh.attributes_sha==reference['attributes_sha256']
        assert mesh.source['bank_sha256']==reference['source_bank_sha256'];bindings[mesh.source['bank_path']]=mesh.source['bank_sha256']
        row=dict(status,registered_geometry_binding=dict(registry_revision=REVISION,registry_mesh_id=reference['registry_mesh_id'],
            geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,bank_path=mesh.source['bank_path'],bank_sha256=mesh.source['bank_sha256'],
            world_frame=frame,actual_record_decoded_once=True,coordinate_rule=mesh.source['coordinates']),
            clinical_or_native_landmark_registration_established=False,
            source_frame_resolution='Literal actual installed world-metre coordinates are verified. Provider name or generic nativePartFrame does not certify anatomical registration.',
            source_surface_is_not_automatically_a_point_or_nonphysical_reference=True)
        if before is None or pid not in before_parts:
            row['explicit_global_height_binding']=dict(verified=False,reason='No literal same-ID source exists in the saved pre-height catalogue')
        else:
            old=reader.read_part(before,before_parts[pid],frame='historical-pre-height-source-metres',catalogue_path=before_path,catalogue_sha=height['previousAtlasSha256'])
            same=old.vertices.shape==mesh.vertices.shape and np.array_equal(old.faces,mesh.faces) and old.source['source_vertex_ids']==mesh.source['source_vertex_ids']
            error=float(np.max(np.linalg.norm(old.vertices*height['scale']+np.asarray(height['translationMeters'])-mesh.vertices,axis=1))) if same else None
            tolerance=float(height['maximumPositionRoundingMeters'])*np.sqrt(3)+1e-10
            row['explicit_global_height_binding']=dict(verified=bool(same and error<=tolerance),source_catalogue=str(before_path),source_catalogue_sha256=height['previousAtlasSha256'],
                source_geometry_sha256=old.geometry_sha,same_source_triangles_and_controls=same,uniform_scale=height['scale'],translation_m=height['translationMeters'],
                maximum_world_vertex_residual_m=error,roundoff_tolerance_m=tolerance,
                scope='Replay of a saved whole-model coordinate transform, not independent anatomical registration to ZA landmarks.')
        if p.get('hairField'):
            hair=p['hairField'];hp=reader.asset(hair['url']);assert digest_file(hp)==hair['sha256'];bindings[str(hp)]=hair['sha256']
            raw=np.memmap(hp,dtype='u1',mode='r');roots=reader.view(raw,hair['positions'],hair['count'],'<f4',3).astype(float)
            row['hair_drawing_binding']=dict(review_source=p.get('reviewSource'),field_path=str(hp),field_sha256=hair['sha256'],
                shaft_count=hair['count'],diameter_m=hair['diameter'],length_m=hair['length'],spacing_m=hair['spacing'],actual_root_bounds_m=[roots.min(0).tolist(),roots.max(0).tolist()],
                carrier_mesh_is_not_rendered_occupied_volume=True,source='app/metric-hair-field.ts: metricHairGeometry constructs a closed cylinder prototype at each literal root; scene.tsx chooses this branch when hairField exists.',
                skin0_hidden_but_not_retired=status['normal_default_visible_skin0'] is False and not status['normal_runtime_suppressed'])
        rows.append(row)
    retained=[r for r in rows if not r['normal_runtime_suppressed']]
    result=dict(schema='near-face-runtime-source-frame-evidence/v1',source_catalogue_sha256=EXPECTED,registry_revision=REVISION,
        actual_runtime_input_sha256=runtime['ordinary_body_sha256'],source_inventory=str(old_context_path),runtime_evidence=str(runtime_path),
        code_bindings=runtime['code_bindings'],source_file_bindings=bindings,rows=rows,
        counts=dict(total=len(rows),genuinely_suppressed_in_normal_loader=sum(r['normal_runtime_suppressed'] for r in rows),
            specific_preferred_replacement=sum(r['exclusion_basis']=='specific-preferred-source-retirement' for r in rows),
            ordinary_loader_default_hidden_alternatives=sum(r['exclusion_basis']=='default-review-alternative-suppression' for r in rows),
            actual_retained_runtime_items=len(retained),actual_default_visible_skin0=sum(r['normal_default_visible_skin0'] for r in rows),
            explicitly_typed_point_or_curve_references=sum(r['declared_point_or_line_reference'] for r in rows),
            verified_literal_height_transforms=sum(r['explicit_global_height_binding']['verified'] for r in rows)),
        retained_runtime_ids=[r['id'] for r in retained],
        conclusion=['Default source exclusions are established by executing the actual preferred-source and hidden-alternative loader rules, not inferred from provenance.',
            'Fifteen retained alternatives can be explicitly viewed with review=alternatives. FJ1578 remains retired by its specific preferred native digastric replacement.',
            'Gingiva, soft palate and the O3M trapezius are actual default physical source meshes; unregistered provenance does not exclude them from current context.',
            'Metric scalp hair remains an active appearance assembly, hidden at the local default skin=0; use its hairField drawing rather than treating the carrier as a solid organ.',
            'No item among these20 is explicitly typed as a point/curve annotation. Finite occupied volume cannot be inferred from an ordinary source mesh alone.'],
        no_collision_or_occupied_volume_acceptance=True,no_visibility_catalogue_source_or_relation_edits=True)
    target=BASE/'nearby-unresolved-participation-frame-v25.json';atomic_json(target,result);assert digest_file(SOURCE)==EXPECTED
    print(json.dumps(dict(output=str(target.relative_to(ROOT)),counts=result['counts'],retained_runtime_ids=result['retained_runtime_ids'],
        height_transform_not_replayed=[r['id'] for r in rows if not r['explicit_global_height_binding']['verified']])))


if __name__=='__main__':main()
