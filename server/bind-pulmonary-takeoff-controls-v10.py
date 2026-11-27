"""Bind actual named airway-root charts and existing vascular root profiles."""
import json
import time
from pathlib import Path
from collections import defaultdict
import numpy as np
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.tube_station_topology import recover_station_rings
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')


def main():
    began=time.monotonic();family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue']);assert digest_file(path)==family['source_catalogue_sha256']
    atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader()
    runtime=json.loads((BASE/'actual-runtime-context-v9/runtime-context.json').read_text())
    interfaces=json.loads((BASE/'source-correspondence-contract-v7/contract.json').read_text())['airway_interface_landmarks']
    fixed=np.load(BASE/'source-correspondence-contract-v7/literal-fixed-airway-interface-controls.npz',allow_pickle=False)
    out=BASE/'takeoff-and-root-controls-v10';out.mkdir();arrays={};rows=[]
    cases=[
        ('right','superior','Apical segmental bronchus of right lung (BI)',['right_apical_bronchus']),
        ('right','superior','Posterior segmental bronchus of right lung (BII)',['right_posterior_bronchus']),
        ('right','superior','Anterior segmental bronchus of right lung (BIII)',['right_anterior_bronchus']),
        ('right','middle','Lateral segmental bronchus of right lung (BIV)',['right_lateral_bronchus']),
        ('right','middle','Medial segmental bronchus of right lung (BV)',['right_medial_bronchus']),
        ('right','inferior','Superior segmental bronchus of right lung (BVI)',['right_superior_bronchus']),
        ('right','inferior','Medial basal segmental bronchus of right lung (BVII)',['right_medial_basal_bronchus']),
        ('right','inferior','Anterior basal segmental bronchus of right lung (BVIII)',['right_anterior_basal_bronchus']),
        ('right','inferior','Lateral basal segmental bronchus of right lung (BIX)',['right_lateral_basal_bronchus']),
        ('right','inferior','Posterior basal segmental bronchus of right lung (BX)',['right_posterior_basal_bronchus']),
        ('left','superior','Apicoposterior segmental bronchus of left lung (BI+BII)',['left_apical_bronchus','left_posterior_bronchus']),
        ('left','superior','Anterior segmental bronchus of left lung (BIII)',['left_anterior_bronchus']),
        ('left','superior','Superior lingular segmental bronchus of left lung (BIV)',['left_lingula_superior_bronchus']),
        ('left','superior','Inferior lingular segmental bronchus of left lung (BV)',['left_lingula_inferior_bronchi']),
        ('left','inferior','Superior segmental bronchus of left lung (BVI)',['left_superior_bronchus']),
        ('left','inferior','Medial basal segmental bronchus of left lung (BVII)',['left_medial_basal_bronchus']),
        ('left','inferior','Anterior basal segmental bronchus of left lung (BVIII)',['left_anterior_basal_bronchus']),
        ('left','inferior','Lateral basal segmental bronchus of left lung (BIX)',['left_lateral_basal_bronchus']),
        ('left','inferior','Posterior basal segmental bronchus of left lung (BX)',['left_posterior_basal_bronchus']),
    ]
    for ci,(side,lobe,source_name,targets)in enumerate(cases):
        ident='ZA:'+source_name;m=reader.read_part(atlas,parts[ident],frame='model:male:assembled-metres');v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces]
        tubes=recover_station_rings(v,f);mouths=[];paths=[]
        for ti,tube in enumerate(tubes):
            sections=v[np.array(tube['rings'])];centre=sections.mean(1);radius=np.linalg.norm(sections-centre[:,None],axis=2).mean(1);paths.append(centre)
            for endpoint in [0,-1]:mouths.append(dict(component=ti,end=endpoint,points=sections[endpoint],center=centre[endpoint],radius=float(radius[endpoint])))
        mouths.sort(key=lambda r:-r['radius']);root=mouths[0];unique_largest=len(mouths)==1 or root['radius']>mouths[1]['radius']*1.00001
        target_boundaries=[];target_ids=['LOCAL:male:pulmonary:bronchus:'+x for x in targets]
        parent='LOCAL:male:pulmonary:bronchus:'+side+'_'+lobe+'_lobar_bronchus'
        for pid in target_ids:
            found=[r for r in interfaces if r['part_id']==pid and parent in r['other_branch_ids'] and r['all_edges_have_literal_opposed_other_branch']]
            if len(found)==1:target_boundaries.append(found[0])
        source_root_key=f'case{ci}_source_entry_ring';arrays[source_root_key]=root['points']
        target_key=f'case{ci}_target_named_entry_charts'
        if target_boundaries:arrays[target_key]=np.concatenate([fixed[r['array_key']]for r in target_boundaries])
        qualified=unique_largest and len(target_boundaries)==len(target_ids)
        row=dict(source_id=ident,source_geometry_sha256=m.geometry_sha,side=side,lobe=lobe,
                 source_unique_largest_entry_profile=unique_largest,source_entry_component=root['component'],source_entry_end=root['end'],
                 source_entry_radius_m=root['radius'],source_entry_array=source_root_key,
                 target_ids=target_ids,target_parent_id=parent,target_actual_boundary_charts=target_boundaries,
                 source_point_m=root['center'].tolist(),target_point_m=arrays[target_key].mean(0).tolist()if target_boundaries else None,
                 target_array=target_key if target_boundaries else None,qualified_initialization_landmark=qualified,
                 same_branch_name_not_current_vascular_registration=True,
                 combined_BI_BII_correspondence=len(target_ids)>1,
                 geometry_not_changed=True)
        rows.append(row)
    # Retain full literal arterial root profiles and identify genuine shared
    # edge coordinates. No projected point is promoted into a tissue opening.
    blood_ids=[r['id']for r in runtime['sourceRows']];blood_mesh={pid:reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres')for pid in blood_ids}
    edge_map=defaultdict(list)
    for pid,m in blood_mesh.items():
        for face in m.faces:
            pp=[tuple(x)for x in m.vertices[face]]
            for a,b in zip(pp,pp[1:]+pp[:1]):edge_map[tuple(sorted((a,b)))].append((pid,a,b))
    root_rows=[]
    for pi,pid in enumerate(['ZA:Left pulmonary artery','ZA:Right pulmonary artery','ZA:Bifurcation of pulmonary trunk','ZA:Pulmonary trunk']):
        m=blood_mesh[pid];v,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces]
        try:loops=boundary_loops(f)
        except ValueError as error:root_rows.append(dict(part_id=pid,unsupported_boundary=str(error)));continue
        for li,loop in enumerate(loops):
            pp=v[np.asarray(loop,int)];others=[];opposed=[]
            for a,b in zip(pp,np.roll(pp,-1,axis=0)):
                aa,bb=tuple(a),tuple(b);rr=[x for x in edge_map[tuple(sorted((aa,bb)))]if x[0]!=pid];others.extend(x[0]for x in rr);opposed.append(any(x[1]==bb and x[2]==aa for x in rr))
            key=f'fixed_native_PA_{pi}_profile_{li}';arrays[key]=pp
            root_rows.append(dict(part_id=pid,boundary=li,profile_array=key,profile_points=len(pp),centre_m=pp.mean(0).tolist(),
                                  shared_other_part_ids=sorted(set(others)),opposed_edges=sum(opposed),edges=len(pp),
                                  complete_literal_opposed_profile=bool(all(opposed)),root_or_distal_role_not_assumed_from_centroid=True))
    np.savez_compressed(out/'literal-correspondence-and-root-profiles.npz',**arrays)
    result=dict(schema='named-pulmonary-entry-chart-and-protected-root-control-bindings/v10',
                source_catalogue_sha256=family['source_catalogue_sha256'],ordinary_runtime_qualification_sha256=digest_file(BASE/'actual-runtime-context-v9/runtime-context.json'),
                arrays_sha256=digest_file(out/'literal-correspondence-and-root-profiles.npz'),named_control_pairs=rows,
                native_arterial_root_profiles=root_rows,whole_runtime_qualified_blood_source_ids=blood_ids,
                fixed_heart_skeleton_and_coronaries=True,all_seven_lingular_pieces_explicitly_in_closure=True,
                candidate_not_run=True,geometry_changed=False,catalogue_changed=False,anatomical_acceptance=False,
                constraints=['Qualified named airway-root charts only initialize a new source correspondence field.',
                             'Current vascular roots/ports still require named structural ownership and protected coordinate arrays.',
                             'No source bronchial or vascular opening is inferred from an average boundary point.'],
                elapsed_seconds=time.monotonic()-began)
    (out/'evidence.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(out),named_pairs=len(rows),qualified_pairs=sum(r['qualified_initialization_landmark']for r in rows),
        arterial_profiles=len(root_rows),complete_shared_profiles=sum(r.get('complete_literal_opposed_profile',False)for r in root_rows),
        evidence_sha256=digest_file(out/'evidence.json'),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':main()
