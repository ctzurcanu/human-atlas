"""Declare fixed anatomy, real airway interfaces and unaccepted fit hypotheses."""
import json
import time
from pathlib import Path
from collections import defaultdict
import numpy as np

from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue'])
    assert digest_file(path)==family['source_catalogue_sha256']
    atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    out=BASE/'source-correspondence-contract-v7';out.mkdir();fixed_airway_ids=family['fixed_current_airway_ids']
    meshes={pid:reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres')for pid in fixed_airway_ids}
    edge_owners=defaultdict(list)
    for pid,m in meshes.items():
        for face in m.faces:
            points=[tuple(v)for v in m.vertices[face]]
            for a,b in zip(points,points[1:]+points[:1]):
                edge_owners[tuple(sorted((a,b)))].append((pid,a,b))
    arrays={};interfaces=[];unclassified=[]
    for pi,(pid,m)in enumerate(meshes.items()):
        v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces]
        try:loops=boundary_loops(f)
        except ValueError as error:
            unclassified.append(dict(part_id=pid,reason=str(error)));continue
        for li,loop in enumerate(loops):
            points=v[np.asarray(loop,int)];opposed=[];owners=[]
            for a,b in zip(points,np.roll(points,-1,axis=0)):
                aa,bb=tuple(a),tuple(b);others=[q for q in edge_owners[tuple(sorted((aa,bb)))]if q[0]!=pid]
                opposed.append(any(oa==bb and ob==aa for _,oa,ob in others));owners.extend(q[0]for q in others)
            if owners:
                key=f'airway_{pi}_boundary_{li}';arrays[key]=points
                interfaces.append(dict(part_id=pid,boundary_index=li,array_key=key,controls=len(points),
                                       geometry_sha256=m.geometry_sha,other_branch_ids=sorted(set(owners)),
                                       all_edges_have_literal_opposed_other_branch=bool(all(opposed)),
                                       opposed_edges=sum(opposed),boundary_edges=len(points),
                                       centre_of_actual_boundary_m=points.mean(0).tolist(),
                                       centre_is_anatomical_landmark_not_wall_or_hole=True,
                                       licence=parts[pid].get('license')))
    np.savez_compressed(out/'literal-fixed-airway-interface-controls.npz',**arrays)
    source_audit=json.loads((BASE/'evidence.json').read_text())
    source_ids=[r['id']for r in source_audit['sources']]
    fixed_vascular=['ZA:Pulmonary trunk','ZA:Bifurcation of pulmonary trunk',
                    'ZA:Left inferior pulmonary vein','ZA:Left superior pulmonary vein',
                    'ZA:Right inferior pulmonary vein','ZA:Right superior pulmonary vein']
    main_pa=['ZA:Left pulmonary artery','ZA:Right pulmonary artery']
    moving=[pid for pid in source_ids if pid not in fixed_vascular+main_pa]
    assert len(moving)==39
    contract=dict(schema='constrained-pulmonary-source-correspondence/v7',
        source_catalogue_sha256=family['source_catalogue_sha256'],family_contract_sha256=digest_file(BASE/'family-contract-v2.json'),
        actual_airway_interface_arrays_sha256=digest_file(out/'literal-fixed-airway-interface-controls.npz'),
        initial_fit_not_run=True,new_geometry_not_built=True,anatomy_accepted=False,catalogue_changed=False,
        fixed_whole_heart_and_native_coronaries=True,fixed_skeleton=True,fixed_airway_ids=fixed_airway_ids,
        fixed_territory_ids=family['fixed_current_territory_ids'],fixed_whole_vascular_part_ids=fixed_vascular,
        moving_distal_source_part_ids=moving,
        narrowly_conditional_main_PA_tissue_ownership=[dict(part_id=pid,proximal_and_root_profiles_literal=True,
            distal_recipient_faces_not_yet_declared=True,allowed_operation='Only source-correspondence-established distal wall/port reconstruction; no blind centroid hole or generic collar')for pid in main_pa],
        every_other_master_part_fixed=[pid for pid in parts if pid not in moving+main_pa],
        obstacle_context_policy=dict(all_current_physical_anatomy_required=True,
            current_runtime_selection_qualification_pending=True,all_relevant_skeleton_muscles_nerves_vessels_airway_and_pleural_sources_required=True,
            source_aabb_only_selects_neighbours=True,excluded_default_variants_need_actual_loader_and_provenance_evidence=True,
            unknown_open_volume_classifiers_not_waived=True,entire_displacement_envelope_must_be_covered=True),
        source_calibre_policy=dict(original_section_coordinates_and_areas_bound_to='native-curves-v3/literal-section-centres.npz',
            radii_may_not_shrink_to_make_clearance_pass=True,section_transport_requires_proper_frames=True,
            residual_storage_error_measured_separately=True,finite_inner_envelopes_and_wall_thickness_not_yet_established=True),
        airway_interface_landmarks=interfaces,unclassified_airway_interfaces=unclassified,
        correspondence_semantics=[
            'Source ZA branch takeoff and course features must be matched to the same named fixed HRA airway, with actual proximal boundary and course witnesses.',
            'A closed fixed HRA boundary with literal opposite neighbouring branch edges is a landmark chart, not permission to cut an artery or a blood-lumen certificate.',
            'Left source B I+II and the separate installed HRA apical/posterior carriers must be treated as a combined correspondence, not an invented one-to-one map.',
            'Branch origin and relation to lobar bronchi distinguish recurrent/ascending posterior branches; signed world-Y travel or nearest territory plane alone do not.',
            'Source 12-point arterial mouths remain unpaired until exact recipient patch, orientation, inner/outer profiles and unique ownership are declared.',
        ],
        proposed_arterial_assignment_hypotheses=[
            dict(terms=['A12.2.01.105','A12.2.01.106'],source='ZA:Anterior segmental artery of right lung',
                 tentative_ascending_components=[0,2,3],tentative_ascending_common_source_portion=dict(component=1,stations=[0,23]),
                 tentative_descending_components=[4,5,6,7],unassigned_common_trunk=dict(component=1,stations=[23,36]),
                 status='takeoff/course and original BIII witnesses; needs primary/render anatomical confirmation and current-HRA fit'),
            dict(terms=['A12.2.01.108','A12.2.01.109'],
                 tentative_ascending=dict(source='ZA:Posterior segmental artery of right lung',components=[0,1,2],origin='Interlobar RPA'),
                 tentative_recurrent_descending=dict(source='ZA:Apical segmental artery of right lung',components=[0,1,2,3],origin='Upper/truncus-superior path via source component9,station47'),
                 status='Two origins and original BII course support hypothesis; no corrected source label or geometry installed'),
            dict(terms=['A12.2.01.206','A12.2.01.207'],source='ZA:Anterior segmental artery of left lung',
                 tentative_two_course_groups=[[0,1],[2,3,4,5]],status='Two origins/courses visible; ascending/descending identity not accepted'),
            dict(terms=['A12.2.01.209','A12.2.01.210'],source='ZA:Posterior segmental artery of left lung',
                 tentative_two_course_groups=[[0,2],[1,3,4]],status='Two origins/courses; correct recurrent/ascending naming unaccepted'),
        ],
        prior_failed_candidate_distinction=dict(previous_private_local_female_source_part_ids=['LOCAL:female:pulmonary_trunk:1508','LOCAL:female:pulmonary_trunk:1509'],
            previous_v415_candidate_file_sha256='6876946aff90fb9b09ffe357322fa8d58f30cc6af72e965b71fe7de7ed416b77',
            previous_failed_fit_does_not_test_current_ZA_source=True,no_previous_failed_candidate_promoted=True),
        pre_fit_requirements=['Root reviews actual source branch assignments and named correspondence controls.',
                              'Exact main-PA distal patch faces/root profiles and inherited gaps/crossings declared.',
                              'All current physical neighbours and their actual frame/volume qualifications bound.',
                              'Whole communicating vascular surface repaired before final semantic tissue partition.'],
        primary_reference_urls=['https://pmc.ncbi.nlm.nih.gov/articles/PMC4454828/','https://pmc.ncbi.nlm.nih.gov/articles/PMC11461042/'],
        elapsed_seconds=time.monotonic()-began)
    (out/'contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),moving_distal_parts=len(moving),interface_landmarks=len(interfaces),
        complete_opposed_boundary_landmarks=sum(r['all_edges_have_literal_opposed_other_branch']for r in interfaces),
        unclassified_interfaces=len(unclassified),contract_sha256=digest_file(out/'contract.json'),elapsed_seconds=contract['elapsed_seconds'])))


if __name__=='__main__':main()
