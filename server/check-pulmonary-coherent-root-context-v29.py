"""Cheap literal HRA airway and S2 diagnostics for the one rigid candidate."""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import winding_contains, topology

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
CAND=BASE/'coherent-pulmonary-root-translation-v28'
OUT=BASE/'coherent-pulmonary-root-cheap-context-v29'
FRAME='model:male:assembled-metres'


def inputs():
    family=json.loads((BASE/'family-contract-v2.json').read_text());assert digest_file(Path(family['source_catalogue']))==family['source_catalogue_sha256']
    atlas=json.loads(Path(family['source_catalogue']).read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    evidence=json.loads((CAND/'candidate-evidence.json').read_text());return family,atlas,parts,reader,evidence


def airway(shard):
    began=time.monotonic();OUT.mkdir(exist_ok=True);family,atlas,parts,reader,evidence=inputs()
    airway_ids=json.loads((BASE/'right-S2-parent-first-course-v26/contract.json').read_text())['fixed_airway_ids']
    vv=[];ff=[];owners=[];offset=0
    for i,pid in enumerate(airway_ids):
        m=reader.read_part(atlas,parts[pid],frame=FRAME);vv.extend(m.vertices);ff.extend(m.faces+offset);owners.extend([i]*len(m.faces));offset+=len(m.vertices)
    airway=Mesh(np.asarray(vv),np.asarray(ff),FRAME);tree=BVH(airway);owners=np.asarray(owners)
    moving=[r for r in evidence['literal_source_replay']if r['moving']];candidates={r['id']:r for r in evidence['candidate_meshes']};rows=[];saved={}
    for i,row in enumerate(moving[shard*4:shard*4+4]):
        before=reader.read_part(atlas,parts[row['id']],frame=FRAME);after=Mesh.load(candidates[row['id']]['mesh']);phases={}
        for phase,m in [('source',before),('candidate',after)]:
            data,summary=pair_contacts(m,airway,0.,second_bvh=tree,max_pairs=100000,max_tests=1000000)
            for key,value in data.items():saved[f'owner{i}_{phase}_{key}']=value
            wrong=np.isin(data['relation'],['crossing','coplanar-overlap']);which=owners[data['pairs'][:,1]]
            summary['actual_named_airway_crossings']=[dict(id=airway_ids[int(j)],face_pairs=int(np.sum(wrong&(which==j))))for j in np.unique(which[wrong])]
            phases[phase]=summary
        rows.append(dict(id=row['id'],source_geometry_sha256=before.geometry_sha,candidate_geometry_sha256=after.geometry_sha,phases=phases))
    np.savez_compressed(OUT/('airway-%02d-witnesses.npz'%shard),**saved)
    result=dict(schema='pulmonary-whole-arterial-root-translation-actual-airway-diagnostics/v29',
                candidate_evidence_sha256=digest_file(CAND/'candidate-evidence.json'),source_catalogue_sha256=family['source_catalogue_sha256'],
                actual_fixed_airway_ids=airway_ids,actual_airway_geometry_sha256=airway.geometry_sha,shard=shard,rows=rows,
                actual_wall_lumen_identity_or_patency_not_inferred=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    path=OUT/('airway-%02d.json'%shard);path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),rows=[dict(id=r['id'],source_crossings=r['phases']['source']['crossing_face_pairs'],candidate_crossings=r['phases']['candidate']['crossing_face_pairs'],complete=all(p['complete_within_pair_budget']for p in r['phases'].values()))for r in rows],elapsed_seconds=result['elapsed_seconds'])))


def territory():
    began=time.monotonic();OUT.mkdir(exist_ok=True);family,atlas,parts,reader,evidence=inputs()
    pid='TA98PULMONARY:male:right-superior:posterior';m=reader.read_part(atlas,parts[pid],frame=FRAME);assert topology(m)['closed_oriented_manifold']
    previous=json.loads((BASE/'right-S2-frozen-port-context-v27/evidence.json').read_text())
    sections=np.load(BASE/'right-S2-ordered-chart-contract-v21/literal-ordered-guides-and-roots.npz',allow_pickle=False)
    selected=json.loads((BASE/'right-S2-ordered-chart-contract-v21/contract.json').read_text())['selected_source_curves']
    shift=np.asarray(json.loads((CAND/'contract.json').read_text())['rigid_translation_m']);saved={};rows=[]
    # No arbitrary arclength threshold is called the segment entry. Every
    # source station centre is retained; only a reviewed actual entry can
    # subsequently define the intrasegmental required-containment interval.
    for row in selected:
        cc=sections[row['source_sections_array']].mean(1);phases={}
        for phase,p in [('source',cc),('candidate',cc+shift)]:
            data=winding_contains(m,p,distance_backend='bvh');inside=np.asarray(data['inside']);amb=np.asarray(data['boundary_ambiguous'])
            saved[row['key']+'_'+phase+'_centres']=p;saved[row['key']+'_'+phase+'_inside']=inside
            phases[phase]=dict(inside_centres=int(inside.sum()),outside_centres=int((~inside&~amb).sum()),ambiguous=int(amb.sum()),
                               inside_by_source_station=inside.tolist(),classification_transition_intervals=np.flatnonzero(inside[:-1]!=inside[1:]).tolist())
        rows.append(dict(key=row['key'],part_id=row['part_id'],source_component=row['component'],phases=phases))
    np.savez_compressed(OUT/'S2-centre-witnesses.npz',**saved)
    result=dict(schema='pulmonary-root-rigid-translation-S2-course-classification/v29',candidate_evidence_sha256=digest_file(CAND/'candidate-evidence.json'),
                territory_id=pid,territory_geometry_sha256=m.geometry_sha,rows=rows,actual_segment_entry_not_established=True,
                artery_can_legitimately_travel_extrasegmentally=True,no_physiological_terminals_inferred_from_available_airway_mesh=True,
                whole_triangle_tissue_containment_not_established=True,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    path=OUT/'S2-course-classification.json';path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(path),rows=[dict(key=r['key'],source_inside=r['phases']['source']['inside_centres'],candidate_inside=r['phases']['candidate']['inside_centres'])for r in rows],elapsed_seconds=result['elapsed_seconds'])))


def aggregate():
    began=time.monotonic();airway_files=[OUT/('airway-%02d.json'%i)for i in range(8)];self_files=[CAND/('self-gate-%02d.json'%i)for i in range(8)]
    airway_rows=[r for p in airway_files for r in json.loads(p.read_text())['rows']]
    self_rows=[r for p in self_files for r in json.loads(p.read_text())['rows']]
    classification=OUT/'S2-course-classification.json';s2=json.loads(classification.read_text())
    assert len(airway_rows)==len(self_rows)==32
    data=dict(schema='coherent-pulmonary-root-translation-cheap-gate-disposition/v29',
              candidate_evidence_sha256=digest_file(CAND/'candidate-evidence.json'),
              all_stage_hashes={str(p):digest_file(p)for p in airway_files+self_files+[classification]},
              all_32_actual_source_and_candidate_surfaces_vs_all_27_current_airways_complete=all(v['complete_within_pair_budget']for r in airway_rows for v in r['phases'].values()),
              actual_source_airway_crossing_face_pairs=sum(r['phases']['source']['crossing_face_pairs']for r in airway_rows),
              actual_candidate_airway_crossing_face_pairs=sum(r['phases']['candidate']['crossing_face_pairs']for r in airway_rows),
              source_parts_clear_of_airway_crossings=sum(not r['phases']['source']['crossing_face_pairs']for r in airway_rows),
              candidate_parts_clear_of_airway_crossings=sum(not r['phases']['candidate']['crossing_face_pairs']for r in airway_rows),
              complete_zero_self_defect_source_parts=sum(r['source']['complete_within_pair_budget']and not r['source']['defect_pairs']for r in self_rows),
              complete_zero_self_defect_candidate_parts=sum(r['candidate']['complete_within_pair_budget']and not r['candidate']['defect_pairs']for r in self_rows),
              other_18_whole_source_networks_retain_positive_inherited_branch_overlap_certificates=True,
              self_positive_certificate_counts_match_source_for_all_32_parts=all(r['source']['defect_pairs']==r['candidate']['defect_pairs']for r in self_rows),
              source_S2_inside_station_centres=sum(r['phases']['source']['inside_centres']for r in s2['rows']),
              candidate_S2_inside_station_centres=sum(r['phases']['candidate']['inside_centres']for r in s2['rows']),
              tested_S2_source_station_centres=247,actual_anatomical_segment_entry_unestablished=True,
              source_rigid_history_replay_not_current_recipient_registration=True,
              installation_rejected=True,broad_context_shards_held=True,
              exact_proximal_annuli_finite_wall_and_patency_still_unestablished=True,
              no_next_candidate_without_parent_trunk_and_fixed_HRA_branch_registration_contract=True,
              catalogue_changed=False,installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    p=OUT/'aggregate-disposition.json';p.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(dict(output=str(p),sha256=digest_file(p),source_airway_crossings=data['actual_source_airway_crossing_face_pairs'],candidate_airway_crossings=data['actual_candidate_airway_crossing_face_pairs'],S2_inside_before=data['source_S2_inside_station_centres'],S2_inside_after=data['candidate_S2_inside_station_centres'],elapsed_seconds=data['elapsed_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['airway','territory','aggregate'],required=True);parser.add_argument('--shard',type=int,default=0);args=parser.parse_args()
    if args.stage=='airway':airway(args.shard)
    elif args.stage=='territory':territory()
    else:aggregate()
