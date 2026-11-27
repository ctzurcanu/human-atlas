"""Bounded actual-current-context diagnostics for the single S2 candidate.

All 1,524 runtime participants retain surface diagnostics. Territory envelopes
have containment roles; attachment references do not silently become occupied
tissue. Same-owner selected faces belong to the separate moving self audit.
"""
import argparse
import json
import math
import time
from pathlib import Path
import numpy as np
from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.spatial import BVH,box_distance_squared
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology,winding_contains

BASE=Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT=BASE/'right-S2-actual-current-context-v24'


def inputs():
    family=json.loads((BASE/'family-contract-v2.json').read_text());path=Path(family['source_catalogue']);assert digest_file(path)==family['source_catalogue_sha256']
    atlas=json.loads(path.read_text());parts={p['id']:p for p in atlas['parts']};participation=json.loads((BASE/'full-context-participation-v11/evidence.json').read_text())
    assert digest_file(Path('.local-models/organ-models/spinal-cord/male-body.json'))==participation['ordinary_body_sha256']
    for path,sha in participation['loader_bindings'].items():assert digest_file(Path(path))==sha
    candidate=json.loads((BASE/'right-S2-terminal-bound-candidate-v22/candidate-evidence.json').read_text());assert candidate['source_catalogue_sha256']==family['source_catalogue_sha256']
    return family,atlas,parts,participation,candidate,AtlasReader()


def moving(family,atlas,parts,candidate,reader):
    before=[];after=[];ff=[];dist=[];face_fragments=[];fragments=[];offset=0
    transport=np.load(BASE/'right-S2-terminal-bound-candidate-v22/literal-proper-ring-transport.npz',allow_pickle=False)
    contract=json.loads((BASE/'right-S2-ordered-chart-contract-v21/contract.json').read_text());rows=contract['selected_source_curves']
    selected={}
    for i,row in enumerate(candidate['selected_fragment_meshes']):
        new=Mesh.load(row['mesh']);full=reader.read_part(atlas,parts[row['part_id']],frame=new.frame);faces=np.asarray(new.source['source_face_ids'],int)
        used,inv=np.unique(full.faces[faces],return_inverse=True);old=full.vertices[used];assert np.array_equal(inv.reshape(-1,3),new.faces)
        source=next(r for r in rows if r['part_id']==row['part_id']and r['component']==row['component']);key=source['key']
        sections=transport[key+'_source_sections'];distance=transport[key+'_root_distance_m'];lookup={tuple(p):d for ring,d in zip(sections,distance)for p in ring}
        dd=np.array([lookup[tuple(p)]for p in old]);before.extend(old);after.extend(new.vertices);dist.extend(dd)
        ff.extend(new.faces+offset);face_fragments.extend([i]*len(new.faces));offset+=len(old)
        selected.setdefault(row['part_id'],[]).extend(faces.tolist());fragments.append(dict(part_id=row['part_id'],component=row['component'],source_face_ids=faces.tolist()))
    attrs=dict(root_distance_m=np.asarray(dist));source=dict(fragments=fragments,face_fragment_indices=face_fragments)
    return Mesh(np.asarray(before),np.asarray(ff,int),'model:male:assembled-metres',attrs,source),Mesh(np.asarray(after),np.asarray(ff,int),'model:male:assembled-metres',attrs,source),selected


def context_mesh(atlas,part,reader,selected):
    m=reader.read_part(atlas,part,frame='model:male:assembled-metres')
    excluded=np.array(selected.get(part['id'],[]),int)
    if not len(excluded):return m,excluded
    keep=np.ones(len(m.faces),bool);keep[excluded]=False
    if not keep.any():return None,excluded
    return Mesh(m.vertices,m.faces[keep],m.frame,m.attributes,dict(m.source,excluded_currently_moving_source_face_ids=excluded.tolist())),excluded


def index():
    began=time.monotonic();family,atlas,parts,participation,candidate,reader=inputs();OUT.mkdir()
    manifest=json.loads((BASE/'full-context-participation-v11/bank-content-hashes-v12.json').read_text())
    for row in manifest['banks']:assert digest_file(Path(row['path']))==row['sha256']
    old,new,selected=moving(family,atlas,parts,candidate,reader);old.save(OUT/'actual-before-moving-support.npz');new.save(OUT/'actual-candidate-moving-support.npz')
    box=np.array([np.minimum(old.bounds[0],new.bounds[0]),np.maximum(old.bounds[1],new.bounds[1])])
    runtime=json.loads((BASE/'actual-runtime-context-v9/runtime-context.json').read_text());declared=np.asarray(runtime['contextBroadphaseBox_m'])
    assert np.all(box[0]>=declared[0])and np.all(box[1]<=declared[1]),'Candidate support escaped the qualified full context envelope'
    records=[];overlapping=[]
    for source in participation['records']:
        if not source['retained']:continue
        pid=source['id'];mesh,excluded=context_mesh(atlas,parts[pid],reader,selected)
        row=dict(id=pid,role=source['role'],source_bank=source['bank_path'],source_record_identity=source['source_record_identity'],
          moving_source_face_ids_removed_from_same_owner_context=excluded.tolist(),source_role_does_not_establish_occupied_tissue=True)
        if mesh is None:row.update(actual_fixed_faces=0,participation='all owned faces are separately in the moving/self audit',overlaps_changed_support=False)
        else:
            bounds=mesh.bounds;overlap=box_distance_squared(box[0],box[1],bounds[0],bounds[1])==0.
            row.update(actual_fixed_faces=len(mesh.faces),actual_indexed_bounds_m=bounds.tolist(),actual_context_geometry_sha256=mesh.geometry_sha,
              overlaps_changed_support=bool(overlap),participation='actual triangle diagnostics retained'if overlap else'exact actual triangle AABB excludes changed support')
            if overlap:overlapping.append(pid)
        records.append(row)
    assert len(records)==1524
    result=dict(schema='right-S2-all-current-physical-context-actual-bounds-index/v24',
      source_catalogue_sha256=family['source_catalogue_sha256'],candidate_evidence_sha256=digest_file(BASE/'right-S2-terminal-bound-candidate-v22/candidate-evidence.json'),
      full_context_evidence_sha256=digest_file(BASE/'full-context-participation-v11/evidence.json'),current_bank_content_hashes_sha256=digest_file(BASE/'full-context-participation-v11/bank-content-hashes-v12.json'),
      all_201_bank_contents_reverified=True,actual_changed_support_union_bounds_m=box.tolist(),candidate_support_inside_original_qualified_context_box=True,
      retained_current_participants=len(records),overlapping_current_participants=len(overlapping),records=records,
      shards=[overlapping[i:i+8]for i in range(0,len(overlapping),8)],
      no_system_provider_or_reference_prefix_exclusion=True,all_territory_envelopes_retained_as_containment_diagnostics=True,
      same_owner_faces_are_in_moving_self_audit=True,
      catalogue_changed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    (OUT/'index.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),participants=len(records),overlapping=len(overlapping),shards=len(result['shards']),
      index_sha256=digest_file(OUT/'index.json'),elapsed_seconds=result['elapsed_seconds'])))


def contacts(shard):
    began=time.monotonic();family,atlas,parts,participation,candidate,reader=inputs();index=json.loads((OUT/'index.json').read_text())
    selected={}
    for row in json.loads((BASE/'right-S2-ordered-chart-contract-v21/contract.json').read_text())['selected_source_curves']:
        arrays=np.load(BASE/'right-S2-ordered-chart-contract-v21/literal-ordered-guides-and-roots.npz',allow_pickle=False)
        selected.setdefault(row['part_id'],[]).extend(arrays[row['source_faces_array']].tolist())
    before=Mesh.load(OUT/'actual-before-moving-support.npz');after=Mesh.load(OUT/'actual-candidate-moving-support.npz');trees={'before':BVH(before),'candidate':BVH(after)}
    current={r['id']:r for r in index['records']};rows=[];saved={}
    for oi,pid in enumerate(index['shards'][shard]):
        context,_=context_mesh(atlas,parts[pid],reader,selected);assert context.geometry_sha==current[pid]['actual_context_geometry_sha256'];tree=BVH(context);phases={}
        for phase,mesh in [('before',before),('candidate',after)]:
            data,summary=pair_contacts(mesh,context,0.,first_bvh=trees[phase],second_bvh=tree,max_tests=500000,max_pairs=100000)
            for key,value in data.items():saved[f'owner{oi}_{phase}_{key}']=value
            phases[phase]=summary
        rows.append(dict(id=pid,role=current[pid]['role'],context_geometry_sha256=context.geometry_sha,phases=phases,
          accidental_versus_anatomical_crossing_review_pending=True,
          declared_source_attachment_candidate=pid in ['ZA:Apical segmental artery of right lung','ZA:Right pulmonary artery'],
          declared_origin_contacts_not_automatically_waived=True))
    np.savez_compressed(OUT/f'shard-{shard:03d}-actual-witnesses.npz',**saved)
    result=dict(schema='right-S2-actual-current-context-triangle-contact-shard/v24',index_sha256=digest_file(OUT/'index.json'),shard=shard,
      rows=rows,witnesses_sha256=digest_file(OUT/f'shard-{shard:03d}-actual-witnesses.npz'),
      all_pair_budgets_complete=all(p['complete_within_pair_budget']for r in rows for p in r['phases'].values()),
      full_context_acceptance=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    (OUT/f'shard-{shard:03d}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(shard=shard,owners=len(rows),all_pair_budgets_complete=result['all_pair_budgets_complete'],
      before_crossings=sum(r['phases']['before']['crossing_face_pairs']for r in rows),candidate_crossings=sum(r['phases']['candidate']['crossing_face_pairs']for r in rows),
      elapsed_seconds=result['elapsed_seconds'])))


def containment():
    began=time.monotonic();family,atlas,parts,participation,candidate,reader=inputs();pid='TA98PULMONARY:male:right-superior:posterior';territory=reader.read_part(atlas,parts[pid],frame='model:male:assembled-metres')
    topo=topology(territory);results=[];saved={}
    for phase,file in [('before','actual-before-moving-support.npz'),('candidate','actual-candidate-moving-support.npz')]:
        mesh=Mesh.load(OUT/file);eligible=np.flatnonzero(np.all(mesh.attributes['root_distance_m'][mesh.faces]>.020,axis=1));used=np.unique(mesh.faces[eligible]);points=mesh.vertices[used]
        row=dict(phase=phase,strictly_distal_source_faces=len(eligible),tested_literal_vertices=len(used),source_root_prefix_threshold_m=.020)
        saved[phase+'_moving_vertex_ids']=used;saved[phase+'_points']=points
        if topo['closed_oriented_manifold']:
            classified=winding_contains(territory,points,distance_backend='bvh');inside=np.asarray(classified['inside']);amb=np.asarray(classified['boundary_ambiguous']);wind=np.array([np.nan if x is None else x for x in classified['winding_number']])
            saved[phase+'_inside']=inside;saved[phase+'_boundary_ambiguous']=amb;saved[phase+'_winding_number']=wind
            row.update(inside_vertices=int(inside.sum()),outside_vertices=int((~inside&~amb).sum()),ambiguous_vertices=int(amb.sum()),actual_all_triangle_solid_angle_winding=True,
              whole_distal_triangle_containment_not_proven_without_crossing_gate=True)
        else:row.update(containment_unestablished='Actual territory is not a closed oriented manifold; no solid volume assumed')
        results.append(row)
    np.savez_compressed(OUT/'S2-territory-containment-witnesses.npz',**saved)
    result=dict(schema='right-S2-current-territory-actual-winding-diagnostics/v24',territory_id=pid,territory_geometry_sha256=territory.geometry_sha,
      territory_topology=topo,phases=results,territory_is_containment_envelope_not_occupied_obstacle=True,
      candidate_evidence_sha256=digest_file(BASE/'right-S2-terminal-bound-candidate-v22/candidate-evidence.json'),
      witness_arrays_sha256=digest_file(OUT/'S2-territory-containment-witnesses.npz'),
      anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    (OUT/'S2-territory-containment.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT/'S2-territory-containment.json'),closed_territory=topo['closed_oriented_manifold'],phases=results,elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['index','contacts','containment'],required=True);parser.add_argument('--shard',type=int);args=parser.parse_args()
    if args.stage=='index':index()
    elif args.stage=='containment':containment()
    else:contacts(args.shard)
