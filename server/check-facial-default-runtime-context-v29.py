"""Held-v17 actual default-retained tissues; no geometry or scope replacement."""
from pathlib import Path
import sys,json,time,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.contact import pair_contacts


def main():
    began=time.monotonic();B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'default-runtime-context-v29';assert not O.exists();O.mkdir()
    runtime_file=D/'nearby-unresolved-participation-frame-v25.json';assert digest_file(runtime_file)=='55728e2322fa4d9e7c9c7c45c5866e2d8cb1c90260044579c7aa66b042b51f7d';runtime=json.loads(runtime_file.read_bytes());expected=runtime['source_catalogue_sha256'];source=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json';assert digest_file(source)==expected
    atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();wall=Mesh.load(D/'finite-wall.npz');origin=wall.bounds.mean(0);z=np.load(D/'surfaces-and-ancestry.npz');ports=json.loads((D/'ports.json').read_bytes());names=ports['component_names']
    spec=importlib.util.spec_from_file_location('default_runtime_context_v29',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper)
    plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());parent,_=helper.occupancy_classifier(Mesh.load(ROOT/plan['fixed_parent']['mesh']),True);assert parent is not None;parent_solid=helper.solid(parent,origin);components={};occupied={}
    for ci,name in enumerate(names[1:],1):
        faces=np.concatenate([z['outer_faces'][z['outer_component']==ci],z['inner_faces'][z['inner_component']==ci],z[f'annulus_{ci+1}_faces']]);components[name]=Mesh(z['vertices'],faces,wall.frame);fields=np.load(D/(name+'-fields.npz'));classifier,_=helper.occupancy_classifier(Mesh(fields['fields'][0],fields['faces'],wall.frame),True);assert classifier is not None;occupied[name]=helper.solid(classifier,origin)-parent_solid
    tested=[];crossings=[];overlaps=[];unknown=[]
    for pid in ['ZA:Gingiva','ZA:Soft palate','O3M:Trapezius muscle.l']:
        evidence=next(r for r in runtime['rows']if r['id']==pid);assert not evidence['normal_runtime_suppressed']and evidence['normal_visible_skin1'];geometry=evidence['registered_geometry_binding'];assert geometry['registry_revision']==runtime['registry_revision']
        target=reader.read_part(atlas,parts[pid],frame=wall.frame,catalogue_path=source,catalogue_sha=expected);assert target.geometry_sha==geometry['geometry_sha256']and target.source['bank_sha256']==geometry['bank_sha256'];classifier,mode=helper.occupancy_classifier(target,False);solid=helper.solid(classifier,origin)if classifier is not None else None
        for name in names[1:]:
            if time.monotonic()-began>48:raise RuntimeError('48-second default-retained context work budget reached')
            contacts,summary=pair_contacts(components[name],target,0.,max_pairs=1000000,max_tests=12000000);assert summary['complete_within_pair_budget'];intentional=pid=='ZA:Soft palate'and name=='ascending-palatine'
            row=dict(part_id=pid,branch=name,geometry_sha256=target.geometry_sha,source_bank_sha256=geometry['bank_sha256'],actual_runtime_world_coordinates_verified=True,independent_native_landmark_registration_not_inferred=True,actual_emitted_outer_inner_and_annular_wall_tested=True,summary=summary,intended_palatal_arterial_wall_tissue_entry=intentional,occupied_classifier=mode)
            if summary['crossing_face_pairs']:
                file=O/('contacts-'+str(len(tested))+'.npz');atomic_npz(file,**contacts);row['contacts']=str(file.relative_to(ROOT))
                if not intentional:crossings.append(dict(branch=name,part_id=pid,crossing_face_pairs=summary['crossing_face_pairs']))
            if solid is not None:
                overlap=occupied[name]^solid;assert overlap.status()==md.Error.NoError;row.update(new_conduit_occupied_overlap_mm3=float(overlap.volume()),occupied_intersection_empty=bool(overlap.is_empty()))
                if not intentional and not overlap.is_empty():overlaps.append(dict(branch=name,part_id=pid,occupied_overlap_mm3=float(overlap.volume())))
            elif not intentional:unknown.append(dict(branch=name,part_id=pid,reason=mode))
            tested.append(row);atomic_json(O/'progress.json',dict(completed_pairs=len(tested),rows=tested,elapsed_seconds=time.monotonic()-began),overwrite=True)
    exclusions=[dict(id=r['id'],runtime_suppressed=r['normal_runtime_suppressed'],exclusion_basis=r['exclusion_basis'],original_explicit_selection_in_normal_mode=r['original_renders_when_explicitly_selected_in_normal_mode'],preferred_replacements=r.get('preferred_replacements'),drawing_kind=r['drawing_kind'])for r in runtime['rows']if r['normal_runtime_suppressed']]
    assert len(exclusions)==16
    prior=json.loads((D/'registered-context-v25/evidence.json').read_bytes());known=prior['actual_crossings']+crossings;knownoverlap=prior['actual_occupied_overlaps']+overlaps
    result=dict(schema='held-facial-current-default-context/v29',source_catalogue_sha256=expected,registry_revision=runtime['registry_revision'],runtime_evidence=str(runtime_file.relative_to(ROOT)),runtime_evidence_sha256=digest_file(runtime_file),candidate_geometry_sha256=wall.geometry_sha,complete_pairs=len(tested),rows=tested,known_actual_crossings=known,known_occupied_overlaps=knownoverlap,new_unverified_physical_volumes=unknown,excluded_default_alternative_ids=exclusions,exclusion_basis='Exact executed ordinary-body loader proof; excludes only normally suppressed alternatives/retired source, not physical anatomy based on provider names.',scalp_hair_scope='Actual closed instanced shaft/root bounds are audited separately; the macro carrier is not an occupied solid. skin=0 hides the instanced hair; skin=1 must use the emitted instance envelope.',oral_fluid_exclusion='Intended ascending-palatine entry into actual soft-palate tissue is not evidence of absence from the oral/pharyngeal fluid space. This audit has not established a complete native closed airspace.',source_aware_gate_pass=False,held_candidate_unchanged=True,packing_performed=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result)
    routes=[r for r in plan['routes']if r['name']in{'tonsillar','glandular-0','glandular-1'}]
    for route in routes:assert digest_file(ROOT/route['guide'])==route['guide_sha256']
    scope=dict(schema='unbuilt-narrower-facial-source-binding/v29',source_catalogue=str(source.relative_to(ROOT)),source_catalogue_sha256=expected,source_snapshot_binding=str((D/'source-binding-v19/evidence.json').relative_to(ROOT)),source_snapshot_binding_sha256=digest_file(D/'source-binding-v19/evidence.json'),registry_revision=runtime['registry_revision'],runtime_input_sha256=runtime['actual_runtime_input_sha256'],runtime_evidence_sha256=digest_file(runtime_file),fixed_parent=plan['fixed_parent'],retained_literal_routes=routes,reference_only_full_held_wall_archive=str((D/'finite-wall.npz').relative_to(ROOT)),reference_only_full_held_wall_sha256=digest_file(D/'finite-wall.npz'),reference_full_wall_still_contains_unaccepted_palatine=True,new_narrower_geometry_built=False,new_owner_bank_exists=False,root_review_required_before_new_construction=True,known_scope_failures=[x for x in known if x['branch']!='ascending-palatine'],known_scope_occupied_overlaps=[x for x in knownoverlap if x['branch']!='ascending-palatine'],remaining_physical_volume_unknowns=[x for x in prior['unverified_native_occupied_volumes']+unknown if x['branch']!='ascending-palatine'],oral_pharyngeal_lumen_exclusion_unresolved=True,anatomical_acceptance=False,installed=False)
    atomic_json(O/'unbuilt-tonsillar-glandular-source-binding.json',scope);assert digest_file(source)==expected and Mesh.load(D/'finite-wall.npz').geometry_sha==wall.geometry_sha;print(json.dumps({k:v for k,v in result.items()if k not in{'rows','excluded_default_alternative_ids'}}))


if __name__=='__main__':main()
