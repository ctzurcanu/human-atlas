"""Exact current thoracic-aorta blood compartment; no wall or atlas mutation."""
from pathlib import Path
import argparse
import json
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file,identity
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import audit,winding_contains
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops,close_with_shared_disks
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

BASE=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-thoracic-aorta-wall-v1'
OUT=BASE/'blood-domain-v1'
CAT=ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
EXPECTED='346bedcd97d42e4514448319b101fd4448cab1257457473754bfd369c36971e6'
HOST='TA98WALL:male:thoracic-aorta:intima'
FRAME='model:local-ta98:assembled-metres'


def prepare():
    began=time.monotonic();assert digest_file(CAT)==EXPECTED and not OUT.exists()
    atlas=json.loads(CAT.read_text());parts={p['id']:p for p in atlas['parts']};record=parts[HOST]
    assert record['sectionSolid'] and record['pairedTissue'] and not record.get('suppressed')
    assert record['tissueDrawGroups']==[dict(count=2592,neighbor='TA98WALL:male:thoracic-aorta:media',start=0),dict(count=2592,start=2592),dict(count=72,start=5184),dict(count=72,start=5256)]
    reader=AtlasReader();host=reader.read_part(atlas,record,frame=FRAME)
    source_faces=np.arange(864,1728);side=host.faces[source_faces][:,[0,2,1]]
    used,inverse=np.unique(side,return_inverse=True);positions=host.vertices[used];faces=inverse.reshape(-1,3)
    fields=np.load(BASE/'fields.npz',allow_pickle=False);native=np.load(BASE/'native-controls.npz',allow_pickle=False)
    # Full literal bank replay: historical Float64 construction fields may not
    # replace the current encoded boundary. Every emitted current facet must
    # equal the actually stored innermost field after its original F32 encoding.
    encoded=fields['fields'][3].astype('<f4').astype(float)
    expected_triangles=encoded[fields['faces']]
    assert np.array_equal(positions[faces],expected_triangles),'Current whole inner chart no longer replays the saved source field.'
    error,field_to_current=cKDTree(positions).query(encoded)
    assert np.all(error==0) and len(set(field_to_current))==444
    loops=boundary_loops(faces);assert len(loops)==2 and all(len(q)==12 for q in loops)
    profiles=positions[field_to_current[native['rings']]];centres=profiles.mean(1)
    assert all(any(set(q)==set(field_to_current[native['rings'][station]])for q in loops)for station in [0,36])
    source_ids=np.asarray(host.source['source_vertex_ids'])[used]
    bank=reader.asset(atlas['chunks'][record['chunk']]['url']);raw=reader._bank(bank)
    words=reader.view(raw,record['normals'],record['vertexCount'],'<i2',3)[source_ids].copy()
    signs=-reader.view(raw,record['normalSigns'],record['vertexCount'],'i1')[source_ids].copy()
    uv=reader.view(raw,record['uvs'],record['vertexCount'],'<f4',2)[source_ids].copy()
    bindings={str(CAT):EXPECTED,str(bank):digest_file(bank),str(BASE/'fields.npz'):digest_file(BASE/'fields.npz'),
              str(BASE/'plan.json'):digest_file(BASE/'plan.json'),str(BASE/'native-controls.npz'):digest_file(BASE/'native-controls.npz'),
              str(BASE/'classified-inventory.json'):digest_file(BASE/'classified-inventory.json')}
    for suffix in ['externa','media','intima']:
        current=reader.read_part(atlas,parts['TA98WALL:male:thoracic-aorta:'+suffix],frame=FRAME)
        assert current.source['bank_sha256']==bindings[str(bank)]
    OUT.mkdir();atomic_npz(OUT/'construction-inputs.npz',positions=positions,faces=faces,source_host_face_ids=source_faces,
        source_host_compact_vertex_ids=used,source_raw_vertex_ids=source_ids,normal_words=words,normal_signs=signs,uvs=uv,
        section_profiles=profiles,section_centres=centres,port0=loops[0],port1=loops[1])
    contract=dict(schema='current-thoracic-aortic-blood-domain-contract/v1',source_catalogue_sha256=EXPECTED,term='A12.0.00.009',
        identity='Blood compartment limited to the current thoracic-aortic lumen; no global blood or microstructure completion.',
        frame=FRAME,host=HOST,host_record=record,host_record_sha256=identity(record),host_geometry_sha256=host.geometry_sha,
        source_faces=source_faces.tolist(),source_facet_count=864,source_controls=444,source_ports=[q.tolist()for q in loops],
        all_current_inner_facets_exactly_match_encoded_source_field=True,source_fields_are_not_used_unencoded=True,
        no_wall_geometry_or_fields_may_move=True,fluid_frontier_policy=dict(expected_ports=2,port_control_count=12,
        caps_are_computational_blood_frontiers=True,render_caps_as_tissue=False,render_caps_by_default=False,adjoining_patency_verified=False),
        input_bindings=bindings,inputs_sha256=digest_file(OUT/'construction-inputs.npz'),
        review_scope_path=str(ROOT/'reports/mesh-audit/REVIEW.md'),review_scope_sha256=digest_file(ROOT/'reports/mesh-audit/REVIEW.md'),
        installed=False,anatomical_acceptance=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'contract-before-construction.json',contract)
    print(json.dumps(dict(prepared=str(OUT),source_controls=444,source_facets=864,elapsed_seconds=contract['elapsed_seconds'])))


def construct():
    began=time.monotonic();c=json.loads((OUT/'contract-before-construction.json').read_text())
    assert not (OUT/'quality.json').exists()
    for path,expected in c['input_bindings'].items():assert digest_file(path)==expected
    assert digest_file(OUT/'construction-inputs.npz')==c['inputs_sha256']
    x=np.load(OUT/'construction-inputs.npz',allow_pickle=False)
    v,f,caps=close_with_shared_disks(x['positions'],x['faces'],expected_mouths=2)
    assert np.array_equal(v[:444],x['positions']) and np.array_equal(f[:864],x['faces'])
    mesh=Mesh(v,f,FRAME,source=dict(kind='actual-current-thoracic-aortic-fluid',catalogue_sha256=EXPECTED,
        host=HOST,host_bank_sha256=c['input_bindings'][str(BASE/'packed/cells.bin')],tissue=False,cap_faces_are_computational_only=True))
    mesh.save(OUT/'blood-domain.npz');check=audit(mesh,require_closed=True,expected_components=1,expected_genus=[0])
    positive=float(check['topology']['oriented_triangle_integral_m3'])>0
    contains=winding_contains(mesh,x['section_centres'][1:-1])
    inside=all(contains['inside']) and not any(contains['boundary_ambiguous'])
    fluid,_=make_solid(v,f);atlas=json.loads(CAT.read_text());reader=AtlasReader();parts={p['id']:p for p in atlas['parts']}
    overlaps=[]
    for suffix in ['externa','media','intima']:
        pid='TA98WALL:male:thoracic-aorta:'+suffix;wall=reader.read_part(atlas,parts[pid],frame=FRAME)
        solid,_=make_solid(wall.vertices,wall.faces);overlap=float((fluid^solid).volume())
        overlaps.append(dict(part_id=pid,wall_geometry_sha256=wall.geometry_sha,wall_attributes_sha256=wall.attributes_sha,
                             intersection_volume_m3=overlap))
    # Per-corner fields preserve the actual source words and UVs. A separate
    # flat normal chart belongs only to computational fluid-port caps.
    side_words=x['normal_words'][x['faces']];side_signs=x['normal_signs'][x['faces']];side_uv=x['uvs'][x['faces']]
    ct=v[caps];cn=np.cross(ct[:,1]-ct[:,0],ct[:,2]-ct[:,0]);cn/=np.linalg.norm(cn,axis=1)[:,None]
    atomic_npz(OUT/'boundary-corner-fields.npz',normal_words=side_words,normal_signs=side_signs,uvs=side_uv,
               world_positions=v[f[:864]],source_host_face_ids=x['source_host_face_ids'],
               computational_cap_face_ids=np.arange(864,len(f)),computational_cap_normal_words=np.repeat(np.rint(cn*32767).astype('<i2')[:,None],3,axis=1))
    atomic_json(OUT/'fluid-frontiers.json',dict(schema='fluid-frontiers/v1',faces=list(range(864,len(f))),
        positions=v[caps],source_port_controls=[x['positions'][x['port0']],x['positions'][x['port1']]],
        tissue=False,render_as_tissue=False,render_by_default=False,adjoining_vessel_continuity_verified=False))
    overlap_pass=all(abs(r['intersection_volume_m3'])<=1e-18 for r in overlaps)
    quality=dict(schema='current-thoracic-aortic-blood-domain-quality/v1',passed=check['passed'] and positive and inside and overlap_pass,
        source_catalogue_sha256=EXPECTED,contract_sha256=digest_file(OUT/'contract-before-construction.json'),candidate_sha256=digest_file(OUT/'blood-domain.npz'),
        term='A12.0.00.009',representation='fluid-compartment',audit=check,positive_volume=positive,
        exact_reciprocal_intimal_facets=864,actual_source_normal_words_UV_literal_with_opposite_signs=True,
        all444wall_boundary_controls_literal=True,fluid_frontiers=2,computational_cap_triangles=len(caps),
        no_tissue_caps_emitted=True,inner_section_centres=contains,wall_overlap=overlaps,wall_overlap_passed=overlap_pass,
        wall_geometry_attributes_and_all_other_anatomy_unchanged=True,source_fields_full_binding_verified=True,
        root_and_branch_patency_unverified=True,histological_thickness_illustrative=True,anatomical_acceptance=False,installed=False,
        elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT/'quality.json',quality)
    print(json.dumps(dict(output=str(OUT),passed=quality['passed'],faces=len(f),volume_m3=check['topology']['oriented_triangle_integral_m3'],wall_overlap=overlaps,elapsed_seconds=quality['elapsed_seconds'])))
    assert quality['passed'],'Stop after actual geometry/containment failure; no corrective loop.'


def pack_volume():
    from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
    from geometry_tools.tissue import decode_packed,validate_cells
    from geometry_tools.common import atomic_bytes
    began=time.monotonic();replay=json.loads((OUT/'existing-cell-validator-replay-v1.json').read_text())
    assert replay['passed'] and replay['candidate_sha256']==digest_file(OUT/'blood-domain.npz')
    contract=json.loads((OUT/'contract-before-construction.json').read_text())
    for path,expected in contract['input_bindings'].items():assert digest_file(path)==expected
    policy=dict(schema='aortic-blood-display-policy/v2',supersedes_only_display_requirement_in_contract_sha256=digest_file(OUT/'contract-before-construction.json'),
        closed_fluid_compartment=True,render_all888_fluid_boundary_faces=True,displayed_truncation_frontiers=2,
        computational_frontiers_are_anatomical_walls=False,frontiers_are_vessel_obstructions=False,
        frontier_qualification='Displayed truncation boundaries of the modeled blood region; adjoining vessel patency is unverified.',
        material=dict(color=[105,15,27],opacity=.6,roughness=.5,metalness=0),pairedTissue=False,tissueWallAssembly=None,
        source_wall_or_pose_changes=False,installed=False)
    atomic_json(OUT/'display-policy-v2.json',policy)
    target=OUT/'packed-blood-v1';assert not target.exists()
    mesh=Mesh.load(OUT/'blood-domain.npz');x=np.load(OUT/'construction-inputs.npz',allow_pickle=False);packet=np.load(OUT/'boundary-corner-fields.npz',allow_pickle=False)
    assert np.all(packet['normal_signs']==1)
    origin=np.array([0.,1.375,0.]);actual=(mesh.vertices-origin).astype('<f4').astype(float)+origin
    surfaces={'blood-intimal-boundary':mesh.faces[:864],'fluid-truncation-frontiers':mesh.faces[864:]}
    ct=actual[surfaces['fluid-truncation-frontiers']];normal=np.cross(ct[:,1]-ct[:,0],ct[:,2]-ct[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    words={'blood-intimal-boundary':packet['normal_words'],'fluid-truncation-frontiers':np.repeat(np.rint(normal*32767).astype('<i2')[:,None],3,axis=1)}
    uvs={'blood-intimal-boundary':packet['uvs'],'fluid-truncation-frontiers':np.zeros((24,3,2),'<f4')}
    refs={'blood':[dict(surface=s,side=1,neighbor=None)for s in surfaces]};pid='TA98FLUID:male:thoracic-aorta:blood'
    packed=pack_finite_cells(mesh.vertices,surfaces,words,refs,{'blood':pid},origin,corner_uvs=uvs)
    decoded=decode_packed(packed,FRAME,origin)['blood']
    assert np.array_equal(decoded.triangles[:864],mesh.triangles[:864])
    assert np.array_equal(decoded.attributes['normals'][decoded.faces[:864]],packet['normal_words'])
    assert np.array_equal(decoded.attributes['uvs'][decoded.faces[:864]],packet['uvs'])
    assert np.all(decoded.attributes['normalSigns'][decoded.faces[:864]]==packet['normal_signs'])
    a=json.loads(CAT.read_text());byid={p['id']:p for p in a['parts']};reader=AtlasReader();cells={'blood':decoded}
    for suffix in ['externa','media','intima']:cells[suffix]=reader.read_part(a,byid['TA98WALL:male:thoracic-aorta:'+suffix],frame=FRAME)
    check=validate_cells(cells);assert check['passed'],check['failures']
    target.mkdir();atomic_bytes(target/'cells.bin',packed.data);decoded.save(target/'blood-domain-decoded.npz')
    record=dict(packed.cells['blood']);groups=record.pop('tissueDrawGroups')
    record.update(conceptId=pid,name='Blood · thoracic-aortic compartment',ta98Term='A12.0.00.009',sourceId='Blood compartment inside current native thoracic aorta',
        system='arterial',material='FLUID:male:thoracic-aorta:blood',sectionSolid=True,anatomicalEntityType='fluid-compartment',
        fluidBoundaryGroups=groups,pairedTissue=False,
        coverageLimitation='Blood compartment within the modeled thoracic aorta only. Displayed end faces are fluid-region truncation frontiers, not tissue walls or vessel obstructions; adjoining patency is unverified.')
    assert 'tissueWallAssembly'not in record
    atomic_json(target/'record.json',record)
    proof=dict(schema='packed-current-aortic-blood-domain/v1',passed=True,source_catalogue_sha256=EXPECTED,
        source_candidate_sha256=digest_file(OUT/'blood-domain.npz'),source_standard_replay_sha256=digest_file(OUT/'existing-cell-validator-replay-v1.json'),
        display_policy_sha256=digest_file(OUT/'display-policy-v2.json'),bank_sha256=digest_file(target/'cells.bin'),bank_bytes=len(packed.data),
        frame=FRAME,origin=origin,record=record,material=policy['material'],actual_encoded_validation=check,
        exact_opposed_intimal_facets=864,source_intimal_coordinates_normal_words_UV_opposed_signs_literal=True,
        source_wall_geometry_attributes_unchanged=True,displayed_non_tissue_fluid_truncation_faces=24,
        candidate_cap_encoding_error_max_m=float(np.linalg.norm(actual[444:]-mesh.vertices[444:],axis=1).max()),
        anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-began)
    atomic_json(target/'packing-proof.json',proof)
    print(json.dumps(dict(output=str(target),passed=True,bank_sha256=proof['bank_sha256'],bytes=len(packed.data),elapsed_seconds=proof['elapsed_seconds'])))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['prepare','construct','pack']);a=p.parse_args()
    {'prepare':prepare,'construct':construct,'pack':pack_volume}[a.stage]()
