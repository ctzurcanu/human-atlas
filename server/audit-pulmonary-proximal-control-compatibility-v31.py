"""Read-only qualification of v13 controls before a different construction."""
import json
import argparse
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT=BASE/'proximal-course-control-compatibility-v31'


def main():
    began=time.monotonic();OUT.mkdir()
    v12=BASE/'diverged-course-controls-v12/evidence.json';v13=BASE/'rooted-course-initialization-v13/contract.json'
    data=json.loads(v12.read_text());contract=json.loads(v13.read_text())
    assert digest_file(Path(contract['source_catalogue']))==contract['source_catalogue_sha256']
    control_arrays=BASE/'rooted-course-initialization-v13/literal-ring-transport-and-jacobians.npz'
    z=np.load(control_arrays,allow_pickle=False);entries=json.loads((BASE/'takeoff-and-root-controls-v10/evidence.json').read_text())['named_control_pairs']
    rows=data['course_controls']+data['common_parent_controls'];source=z['control_source'];target=z['control_target'];flow=z['control_flow']
    assert len(rows)==len(source)==20
    roots=np.asarray([r['centre_m']for r in contract['root_profile_protection']]);radius=np.asarray([r['radius_m']for r in contract['root_profile_protection']]);width=np.asarray([r['transition_m']for r in contract['root_profile_protection']])
    t=np.clip((np.linalg.norm(source[:,None]-roots[None],axis=2)-radius)/width,0,1)
    weights=np.prod(t**3*(10+t*(-15+6*t)),axis=1)
    qualification=[]
    for i,row in enumerate(rows):
        item=dict(index=i,kind=row['kind'],source_ids=row.get('source_ids',[row.get('source_id')]),target_id=row['target_id'],
                  source_chart_centre_m=source[i].tolist(),current_target_chart_centre_m=target[i].tolist(),
                  initial_displacement_mm=float(np.linalg.norm(source[i]-target[i])*1000),
                  v13_endpoint_residual_mm=float(np.linalg.norm(flow[i]-target[i])*1000),v13_source_root_mask_weight=float(weights[i]),
                  hard_vascular_course_position_or_tangent_qualified=False)
        if row['kind']=='named-common-parent-entry':
            item.update(compatible_role='named common-clade proximal initialization',
                        actual_source_common_node_spread_mm=row['maximum_source_common_node_representation_difference_m']*1000,
                        one_source_node_must_not_map_to_separate_child_ports=True,
                        exact_bronchial_bifurcation_or_arterial_port_not_established=True)
        else:
            entry=next(e for e in entries if e['source_id']==row['source_id'])
            common=next((r for r in data['common_parent_controls']if row['source_id']in r['source_ids']),None)
            item.update(compatible_role='same named branch course only; cross-model physical station not yet homologized',
                        source_and_target_object_entry_intrinsic_distance_mm=row['edge_graph_distance_m']*1000,
                        same_distance_from_different_part_boundaries_not_homologous_anatomical_station=True,
                        target_actual_named_entry_charts=entry['target_actual_boundary_charts'])
            if common:
                item.update(source_common_parent_to_course_chart_chord_mm=float(np.linalg.norm(source[i]-np.asarray(common['source_point_m']))*1000),
                            current_common_parent_to_course_chart_chord_mm=float(np.linalg.norm(target[i]-np.asarray(common['target_point_m']))*1000))
        qualification.append(item)
    expansions=[]
    for i in range(len(rows)):
        for j in range(i):
            ds=np.linalg.norm(source[i]-source[j]);dt=np.linalg.norm(target[i]-target[j])
            if ds<.006 and dt>2*ds:
                expansions.append(dict(first=i,second=j,source_chord_mm=float(ds*1000),current_target_chord_mm=float(dt*1000),required_local_chord_expansion=float(dt/ds)))
    result=dict(schema='pulmonary-proximal-course-control-anatomical-role-audit/v31',source_catalogue_sha256=contract['source_catalogue_sha256'],
                bound_inputs={str(p):digest_file(p)for p in [v12,v13,control_arrays,BASE/'takeoff-and-root-controls-v10/evidence.json',BASE/'source-correspondence-contract-v7/contract.json',BASE/'rooted-source-recipient-graph-v17/evidence.json',BASE/'coherent-pulmonary-root-translation-v28/contract.json',BASE/'full-context-participation-v11/evidence.json']},
                controls=qualification,nearby_source_control_expansion_witnesses=expansions,
                v13_method=dict(kind='regularized Gaussian stationary velocity,16 RK4 steps',sigma_mm=45.,regularization=.35,
                                velocity_fit_not_endpoint_position_fit=True,serial_sections='flow centres and proper polar factor of flow Jacobian; not own-vessel tangent transport',
                                nonserial_left_PA='actual vertex flow plus shared-lingular-cut and distal-profile rigid blends',nearest_bronchial_projection_used=False),
                root_mask_does_not_directly_suppress_any_of_the_20_source_controls=bool(np.all(weights==1.)),
                principal_limits=['Only proximal3mm object-entry charts, no distal parenchymal course constraints.',
                                  'Same named source/current branch charts are not homologous physical station positions.',
                                  'Source-common nodes/current separate child ports need parent-owned correspondence, not local expansion forced by labels.',
                                  'Fixed old main-PA proximal profiles are not established interfaces to moved bifurcation.',
                                  'Proper per-ring polar rotation does not align ring planes with its own candidate vessel tangent.'],
                next_construction_contract=['Preserve fixed heart/skeleton/coronaries and actual fixed current HRA tissue; assign real bifurcation outlet profiles as explicit proximal recipients with pending phase reconstruction.',
                                            'Use complete32-part arterial closure, actual source parent-station ownership and nine-piece lingular reciprocal interfaces; reconstruct main-PA course before children.',
                                            'Use actual HRA common-clade junction/course charts as proximal role-qualified guides; child source/current station correspondence must start after their respective actual shared parent-course region.',
                                            'Keep artery incidence separate from bronchial guide topology. Distal arterial extent belongs to actual parenchyma; available airway mesh ends are not targets.',
                                            'Constrain parent/child C2 position, tangent and curvature jets at saved source parameters; transport profiles by own-vessel proper frames without radius scaling.',
                                            'Actual communicating wall ownership/junction reconstruction precedes TA98 partition. First actual source-versus-candidate self/airway/territory gates, then full1524-owner context.'],
                unresolved_venous_registration_and_actual_segment_entries=True,
                no_new_candidate_or_fit=True,geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    path=OUT/'evidence.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),sha256=digest_file(path),controls=len(rows),common_clade_qualified=sum(x['kind']=='named-common-parent-entry'for x in qualification),hard_vascular_point_constraints=0,
                          root_masks_all_one=bool(np.all(weights==1.)),nearby_expansion_witnesses=len(expansions),elapsed_seconds=result['elapsed_seconds'])))


def guards():
    began=time.monotonic();p=BASE/'rooted-course-initialization-v13/contract.json';c=json.loads(p.read_text())
    z=np.load(BASE/'rooted-course-initialization-v13/literal-ring-transport-and-jacobians.npz',allow_pickle=False)
    roots=np.asarray([r['centre_m']for r in c['root_profile_protection']]);radius=np.asarray([r['radius_m']for r in c['root_profile_protection']]);width=np.asarray([r['transition_m']for r in c['root_profile_protection']])
    phases=[]
    for key in ['control_source','control_target','control_flow']:
        rr=np.linalg.norm(z[key][:,None]-roots[None],axis=2);t=np.clip((rr-radius)/width,0,1);g=np.prod(t**3*(10+t*(-15+6*t)),axis=1)
        phases.append(dict(phase=key,controls=[dict(index=int(i),guard_weight=float(g[i]),margin_to_frozen_ball_mm=((rr[i]-radius)*1000).tolist())for i in range(len(g))]))
    result=dict(schema='pulmonary-proximal-control-frozen-root-compatibility/v31',
                predecessor_control_audit_sha256=digest_file(OUT/'evidence.json'),root_contract_sha256=digest_file(p),phases=phases,
                initial_source_weight_one_does_not_exclude_root_guard_effect_along_the_flow=True,
                incompatible_target=dict(index=19,role='left superior lobar common parent entry chart',
                                         fixed_source='ZA:Left pulmonary artery proximal-profile protection ball',
                                         target_depth_inside_frozen_ball_mm=2.809660322587563,
                                         exact_target_not_compatible_with_this_frozen_source_ball=True),
                no_new_flow_or_fit_performed=True,geometry_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    path=OUT/'root-protection-compatibility.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),sha256=digest_file(path),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['controls','guards'],default='controls');args=parser.parse_args()
    main()if args.stage=='controls'else guards()
