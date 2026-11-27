"""Restore the recorded root translation to its complete arterial closure.

This is a single source-history-derived rigid initialization, never an
anatomical acceptance or a claim of finite communicating walls.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.contact import self_contact_summary

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT = BASE / 'coherent-pulmonary-root-translation-v28'
FRAME = 'model:male:assembled-metres'


def candidate():
    began = time.monotonic(); OUT.mkdir(exist_ok=True)
    family = json.loads((BASE / 'family-contract-v2.json').read_text())
    assert digest_file(Path(family['source_catalogue'])) == family['source_catalogue_sha256']
    atlas = json.loads(Path(family['source_catalogue']).read_text()); parts = {p['id']:p for p in atlas['parts']}; reader = AtlasReader()
    old = json.loads((BASE / 'takeoff-and-root-controls-v10/evidence.json').read_text())
    ids = old['whole_runtime_qualified_blood_source_ids']
    proximal = ['ZA:Pulmonary trunk', 'ZA:Bifurcation of pulmonary trunk']
    registrations = [parts[i]['cardiacRootRegistration'] for i in proximal]
    transform = np.asarray(registrations[0]['translationMetres'])
    assert np.array_equal(transform, registrations[1]['translationMetres'])
    for pid in proximal:
        assert np.array_equal(parts[pid]['geometryOrigin'], transform)
    moving = [i for i in ids if parts[i]['system'] == 'arterial' and i not in proximal]
    fixed = [i for i in ids if i not in moving]
    assert len(moving) == 32 and len(fixed) == 22
    for pid in moving:
        assert 'cardiacRootRegistration' not in parts[pid]
    a = np.load(BASE / 'takeoff-and-root-controls-v10/literal-correspondence-and-root-profiles.npz', allow_pickle=False)
    profiles = old['native_arterial_root_profiles']
    rp = next(p for p in profiles if p['part_id'] == 'ZA:Right pulmonary artery' and p['boundary'] == 1)
    bp = next(p for p in profiles if p['part_id'] == proximal[1] and p['boundary'] == 0)
    before = a[rp['profile_array']]; recipient = a[bp['profile_array']]
    shifted = before + transform
    distance = np.linalg.norm(shifted[:, None] - recipient[None], axis=2)
    nearest = distance.argmin(1)
    assert len(np.unique(nearest)) == len(before) == len(recipient)
    gap = float(distance[np.arange(len(before)), nearest].max())
    # Source history is supported by matching centres and calibre, but the
    # existing polygon phases are NOT shared annuli. No whole-tree rotation
    # or relaxed interface tolerance is used to conceal that pending repair.
    assert np.linalg.norm(shifted.mean(0)-recipient.mean(0)) < 1e-6
    profile_replay=[]
    for arterial,ab,bifurcation_boundary in [('ZA:Right pulmonary artery',1,0),('ZA:Left pulmonary artery',0,2)]:
        pp=next(p for p in profiles if p['part_id']==arterial and p['boundary']==ab)
        bp2=next(p for p in profiles if p['part_id']==proximal[1]and p['boundary']==bifurcation_boundary)
        aa=a[pp['profile_array']]+transform;bb=a[bp2['profile_array']]
        aa0=aa-aa.mean(0);bb0=bb-bb.mean(0);options=[]
        for reverse in [False,True]:
            for shift in range(len(bb)):
                permutation=np.roll(np.arange(len(bb))[::-1]if reverse else np.arange(len(bb)),shift)
                u,_,vt=np.linalg.svd(aa0.T@bb0[permutation]);rr=u@np.diag([1.,1.,np.linalg.det(u@vt)])@vt
                raw=float(np.linalg.norm(aa-bb[permutation],axis=1).max())
                residual=float(np.linalg.norm(aa0@rr-bb0[permutation],axis=1).max())
                options.append((raw,residual,permutation,rr))
        raw,residual,permutation,rr=min(options,key=lambda x:x[0])
        def area(points):return np.linalg.norm(np.cross(points,np.roll(points,-1,axis=0)).sum(0))*.5
        profile_replay.append(dict(arterial_part=arterial,arterial_profile_array=pp['profile_array'],bifurcation_profile_array=bp2['profile_array'],
                                  maximum_unmodified_phase_point_gap_m=raw,centre_residual_m=float(np.linalg.norm(aa.mean(0)-bb.mean(0))),
                                  mean_radius_before_m=float(np.linalg.norm(aa0,axis=1).mean()),mean_recipient_radius_m=float(np.linalg.norm(bb0,axis=1).mean()),
                                  area_relative_difference=float(area(bb0)/area(aa0)-1),permutation=permutation.tolist(),
                                  diagnostic_phase_rotation_deg=float(np.degrees(np.linalg.norm(Rotation.from_matrix(rr).as_rotvec()))),
                                  diagnostic_phase_rotated_shape_residual_m=residual,phase_rotation_NOT_applied=True,
                                  literal_shared_annulus=False,reciprocal_interface_construction_pending=True))
    contract = dict(schema='whole-pulmonary-arterial-missed-root-translation-contract/v28',
                    source_catalogue_sha256=family['source_catalogue_sha256'],
                    recorded_root_translations={i:parts[i]['cardiacRootRegistration'] for i in proximal},
                    rigid_translation_m=transform.tolist(), moving_source_ids=moving, fixed_family_ids=fixed,
                    actual_RPA_LPA_bifurcation_profile_replay=profile_replay,
                    every_current_venous_heart_skeleton_lung_airway_source_fixed=True,
                    whole_nine_piece_lingular_dependency_closure_moving=True,
                    current_runtime_qualification_sha256=digest_file(BASE / 'actual-runtime-context-v9/runtime-context.json'),
                    whole_1524_physical_context_required=True,
                    RPA_proximal_rim_to_fixed_bifurcation_outlet=dict(source_profile=rp, recipient_profile=bp,
                        maximum_literal_profile_point_gap_after_recorded_translation_m=gap,
                        permutation=nearest.tolist(), finite_join_or_patency_not_inferred=True),
                    finite_wall_and_exact_valve_annulus_still_unestablished=True,
                    no_fit_coefficient_or_radius_or_curvature_change=True,
                    catalogue_changed=False, anatomical_acceptance=False)
    (OUT / 'contract.json').write_text(json.dumps(contract, indent=2) + '\n')
    meshes = []; replay = []
    for pi, pid in enumerate(ids):
        source = reader.read_part(atlas, parts[pid], frame=FRAME)
        vertices = source.vertices + transform if pid in moving else source.vertices
        mesh = Mesh(vertices, source.faces, FRAME, source.attributes,
                    dict(part_id=pid, source_geometry_sha256=source.geometry_sha,
                         recorded_rigid_translation_m=transform.tolist() if pid in moving else [0.,0.,0.],
                         finite_wall_acceptance=False))
        path = OUT / ('whole-family-part-%02d.npz' % pi); mesh.save(path)
        meshes.append(dict(id=pid, mesh=str(path.resolve()), geometry_sha256=mesh.geometry_sha))
        source_edges = source.triangles[:,1:] - source.triangles[:,:1]
        candidate_edges = mesh.triangles[:,1:] - mesh.triangles[:,:1]
        replay.append(dict(id=pid, moving=pid in moving, source_geometry_sha256=source.geometry_sha,
                           candidate_geometry_sha256=mesh.geometry_sha, source_faces_literal=np.array_equal(source.faces,mesh.faces),
                           all_source_attributes_literal=all(np.array_equal(source.attributes[k],mesh.attributes[k])for k in source.attributes),
                           maximum_translation_replay_error_m=float(np.abs(mesh.vertices-source.vertices-(transform if pid in moving else 0.)).max()),
                           maximum_triangle_edge_replay_error_m=float(np.abs(candidate_edges-source_edges).max())))
    evidence = dict(schema='coherent-pulmonary-arterial-root-translation-candidate/v28',
                    source_catalogue_sha256=family['source_catalogue_sha256'],contract_sha256=digest_file(OUT / 'contract.json'),
                    candidate_meshes=meshes,literal_source_replay=replay,
                    moving_arterial_parts=len(moving),fixed_family_parts=len(fixed),
                    original_arterial_incidence_and_relative_anatomy_preserved=True,
                    source_profiles_calibres_and_own_curve_frames_preserved=True,
                    no_new_intrafamily_distortion_beyond_recorded_Float64_roundoff=True,
                    actual_self_context_render_and_finite_join_gates_pending=True,
                    catalogue_changed=False,installed=False,anatomical_acceptance=False,
                    elapsed_seconds=time.monotonic()-began)
    (OUT / 'candidate-evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),sha256=digest_file(OUT/'candidate-evidence.json'),moving=len(moving),fixed=len(fixed),RPA_root_profile_gap_m=gap,elapsed_seconds=evidence['elapsed_seconds'])))


def self_gate(shard):
    began=time.monotonic(); evidence=json.loads((OUT/'candidate-evidence.json').read_text())
    family=json.loads((BASE/'family-contract-v2.json').read_text());assert digest_file(Path(family['source_catalogue']))==family['source_catalogue_sha256']
    atlas=json.loads(Path(family['source_catalogue']).read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    moving=[r for r in evidence['literal_source_replay']if r['moving']];rows=[]
    candidates={r['id']:r for r in evidence['candidate_meshes']}
    for row in moving[shard*4:shard*4+4]:
        source=reader.read_part(atlas,parts[row['id']],frame=FRAME);mesh=Mesh.load(candidates[row['id']]['mesh'])
        before,_=self_contact_summary(source,max_pairs=200000,max_tests=1000000,stop_on_defect=True)
        after,_=self_contact_summary(mesh,max_pairs=200000,max_tests=1000000,stop_on_defect=True)
        rows.append(dict(id=row['id'],source=before,candidate=after, source_defects_not_repaired=True))
    data=dict(candidate_sha256=digest_file(OUT/'candidate-evidence.json'),shard=shard,rows=rows,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    path=OUT/('self-gate-%02d.json'%shard);path.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),rows=[dict(id=r['id'],source_defects=r['source']['crossing_face_pairs']+r['source']['coplanar_overlap_face_pairs'],candidate_defects=r['candidate']['crossing_face_pairs']+r['candidate']['coplanar_overlap_face_pairs'],source_complete=r['source']['complete_within_pair_budget'],candidate_complete=r['candidate']['complete_within_pair_budget'])for r in rows],elapsed_seconds=data['elapsed_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['candidate','self'],required=True);parser.add_argument('--shard',type=int,default=0);args=parser.parse_args()
    candidate()if args.stage=='candidate'else self_gate(args.shard)
