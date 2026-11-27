"""Guarded three-course ACA/ACom exterior and five-piece TA98 catalogue.

No work runs on import. Root must review the exact emitted candidate, actual
fixed context, surfaces and sections before this constructor can run. This
script writes a candidate only; it never installs or alters the source atlas.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.native_tree_exterior import construct
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.audit import audit, topology, winding_contains
from geometry_tools.model import AtlasReader, Mesh

CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1'
OLD = ['ZA:Anterior cerebral artery.l', 'ZA:Anterior cerebral artery.r',
       'ZA:Anterior communicating artery']
NEW = ['TA98PART:ZA:Anterior cerebral artery.l:A12.2.07.023',
       'TA98PART:ZA:Anterior cerebral artery.l:A12.2.07.033',
       'TA98PART:ZA:Anterior cerebral artery.r:A12.2.07.023',
       'TA98PART:ZA:Anterior cerebral artery.r:A12.2.07.033',
       'TA98PART:ZA:Anterior communicating artery:communicating-course']
TERMS = ['A12.2.07.023', 'A12.2.07.033', 'A12.2.07.023',
         'A12.2.07.033', 'A12.2.07.029']
LIMIT = ('Macroscopic arterial exterior with shared open communicating ostia. '
         'The ACA course is fitted with source-calibre-preserving sections; '
         'the ACom bridge is independently constructed. ICA proximal and '
         'distal branch interfaces remain fixed source continuations. '
         'No finite wall, complete perfusion or microscopic anatomy claim.')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=lambda x:
        x.tolist() if hasattr(x, 'tolist') else str(x)) + '\n')


def classifier(mesh):
    """Temporary closing fans for numerical interior tests, never emitted."""
    vertices, faces = mesh.vertices.tolist(), mesh.faces.tolist()
    for loop in boundary_loops(mesh.faces):
        index = len(vertices)
        vertices.append(mesh.vertices[loop].mean(0).tolist())
        faces.extend([[index, int(loop[(j + 1) % len(loop)]), int(loop[j])]
                      for j in range(len(loop))])
    return Mesh(np.asarray(vertices), np.asarray(faces), mesh.frame)


def require_inside(mesh, points, message):
    result = winding_contains(mesh, np.asarray(points))
    assert all(result['inside']) and not any(result['boundary_ambiguous']), message
    return result


def orient_bridge(mesh, fields, centres):
    """Only the independently generated bridge chart may change winding."""
    assert np.array_equal(mesh.vertices, fields.reshape(-1, 3))
    radial = fields - centres[:, None]
    normal = np.cross(mesh.triangles[:, 1] - mesh.triangles[:, 0],
                      mesh.triangles[:, 2] - mesh.triangles[:, 0])
    alignment = np.einsum('ij,ij->i', normal,
                          radial.reshape(-1, 3)[mesh.faces].mean(1))
    reverse = bool(alignment.sum() < 0)
    if reverse:
        alignment = -alignment
    assert np.all(alignment > 0), 'Locally inverted independent bridge chart'
    faces = mesh.faces[:, ::-1].copy() if reverse else mesh.faces.copy()
    return Mesh(mesh.vertices.copy(), faces, mesh.frame,
                {k: v.copy() for k, v in mesh.attributes.items()},
                dict(mesh.source, independent_bridge_winding_reversed=reverse)), reverse


def run(candidate, approval_file, output):
    started = time.monotonic()
    candidate, output = candidate.resolve(), output.resolve()
    approval = json.loads(approval_file.read_text())
    manifest_path = candidate / 'candidate-manifest.json'
    gates_path = candidate / 'actual-geometry-gates.json'
    sections_path = candidate / 'section-isometry-proof.json'
    manifest = json.loads(manifest_path.read_text())
    gates = json.loads(gates_path.read_text())
    assert gates.get('complete') and gates.get('passed'), 'Candidate geometry gates failed'
    assert approval['candidate_manifest_sha256'] == sha(manifest_path)
    assert approval['candidate_geometry_gates_sha256'] == sha(gates_path)
    assert approval['candidate_section_proof_sha256'] == sha(sections_path)
    for key in ['course_passed', 'root_surface_course_review_passed',
                'root_true_section_review_passed', 'fixed_context_review_passed']:
        assert approval.get(key) is True, 'Missing exact root review: ' + key
    assert approval.get('context_evidence'), 'Actual fixed context evidence required'
    for item in approval['context_evidence']:
        assert sha(item['path']) == item['sha256'], 'Context evidence changed'
    assert sha(CAT) == manifest['source_catalogue_sha256'] == approval['source_catalogue_sha256']
    inventory_path = WORK / 'complete-authority-effective-inventory.json'
    inventory = json.loads(inventory_path.read_text())
    assert inventory['term_count'] == len(inventory['terms']) == 24
    assert inventory['source_catalogue_sha256'] == sha(CAT)
    assert not output.exists(), 'Use a fresh output directory'
    output.mkdir(parents=True)
    write(output / 'construction-state.json', dict(complete=False, passed=False,
          installed=False, anatomical_acceptance=False, source_catalogue_sha256=sha(CAT)))
    (output / 'authority-inventory.json').write_bytes(inventory_path.read_bytes())
    atlas = json.loads(CAT.read_text())
    records = {p['id']: p for p in atlas['parts']}
    assert all(x in records and not records[x].get('suppressed') for x in OLD)
    assert not set(NEW) & set(records), 'Candidate part ID already exists'
    fields = np.load(candidate / 'source-and-candidate-fields.npz')
    meshes = [Mesh.load(candidate / 'aca-left.npz'),
              Mesh.load(candidate / 'aca-right.npz'),
              Mesh.load(candidate / 'acom-open-bridge.npz')]
    frame = meshes[0].frame
    assert all(m.frame == frame for m in meshes)
    meshes[2], reversed_bridge = orient_bridge(meshes[2], fields['bridge_fields'],
                                                fields['bridge_centres'])
    for i, label in enumerate(['left', 'right']):
        source_rings = fields[label + '_rings']
        source = Mesh.load(WORK / ('source-01.npz' if i == 0 else 'source-02.npz'))
        assert np.array_equal(meshes[i].faces, source.faces)
        for ring in [source_rings[0], source_rings[-1]]:
            assert np.array_equal(meshes[i].vertices[ring], source.vertices[ring])
            assert np.array_equal(meshes[i].attributes['normals'][ring], source.attributes['normals'][ring])
            assert np.array_equal(meshes[i].attributes['uvs'][ring], source.attributes['uvs'][ring])
    nv = np.cumsum([0] + [len(m.vertices) for m in meshes])
    nf = np.cumsum([0] + [len(m.faces) for m in meshes])
    composite = Mesh(np.concatenate([m.vertices for m in meshes]),
                     np.concatenate([m.faces + nv[i] for i, m in enumerate(meshes)]),
                     frame, {k: np.concatenate([m.attributes[k] for m in meshes])
                             for k in ['normals', 'uvs']})
    mouths, anchors = [], {}
    bridge_centres = fields['bridge_centres']
    bridge_direction = bridge_centres[-1] - bridge_centres[0]
    bridge_direction /= np.linalg.norm(bridge_direction)
    bridge_radius = float(manifest['bridge_radius_m'])
    explicit_anchors = [bridge_centres[0] - bridge_direction * bridge_radius / 2,
                        bridge_centres[-1] + bridge_direction * bridge_radius / 2]
    # Verify complete continuation cones, not merely their apex coordinates.
    for end, side in [(0, 0), (-1, 1)]:
        ring = fields['bridge_fields'][end]
        point = explicit_anchors[side]
        samples = np.concatenate([ring, point[None], (ring + point) / 2])
        require_inside(classifier(meshes[side]), samples, 'Bridge root/cone not fully buried in ACA')
    for role, mesh in enumerate(meshes):
        loops = boundary_loops(mesh.faces)
        assert len(loops) == 2
        for loop in loops:
            controls = loop + nv[role]
            mouths.append(dict(positions_m=composite.vertices[controls],
                               source_control_ids=controls,
                               kind='original' if role < 2 else 'internal communication'))
            if role == 2:
                centre = mesh.vertices[loop].mean(0)
                side = int(np.argmin(np.linalg.norm(bridge_centres[[0, -1]] - centre, axis=1)))
                anchors[len(mouths) - 1] = explicit_anchors[side]
    assert len(anchors) == 2
    plan = dict(schema='ACA-ACom-shared-exterior-five-part/v1',
                whole_organ='Bilateral source-calibre ACA courses with precommunicating and postcommunicating regions joined through one independently constructed anterior communicating exterior.',
                source_catalogue_sha256=sha(CAT), candidate_manifest_sha256=sha(manifest_path),
                approval_sha256=sha(approval_file), changed_sources=OLD, new_parts=NEW,
                component_roles=['left ACA', 'right ACA', 'independent ACom bridge'],
                original_outer_mouths=4, internal_bridge_roots=2, mouths=mouths,
                anchors=anchors, every_output_face_owned_once=True,
                adjoining_interfaces=[
                    'Left and right ACom side ostia meet the ACA pre/post regional boundary through shared exterior coordinates and normal words; no lumen disk.',
                    'Four source ACA end profiles retain literal positions, normal words and UVs.',
                    'Proximal cerebral ICA interface remains an unjoined source continuation; the measured original 2.376 mm gap is not repaired by this cohort.',
                    'Distal pericallosal and callosomarginal root centres retain their original matching source positions; those separate source mouths are not newly welded.',
                    'Optic and olfactory nerves, brain surfaces, foreign vessels and the skeleton remain fixed context and require independent crossing review.',
                ],
                regional_boundary='Source station 13, at the explicitly fitted ACom interface',
                fixed_anatomy='All outside anatomy, skeleton, optic nerves, basilar artery and source ACA end profiles',
                finite_wall=False, physiological_patency_claim=False,
                storage_origin_m=[0., 1.71875, .03125],
                storage_rationale='Dyadic origin near the candidate. Literal fixed Float32 source mouths subtract exactly; one bank preserves shared controls.',
                installed=False, anatomical_acceptance=False)
    write(output / 'assembly-plan-before-union.json', plan)
    composite.save(output / 'composite-input.npz')
    arrays, proof = construct(composite, [np.arange(nf[i], nf[i + 1]) for i in range(3)],
                              mouths, anchors, wall_seconds=35, max_pairs=600000,
                              max_tests=2000000)
    np.savez_compressed(output / 'union-source-provenance.npz', **arrays)
    write(output / 'native-tree-proof.json', proof)
    assert proof['passed'], proof['failures']
    assert proof['derived_fragments'] == 0, 'Temporary conical closure remains visible'
    vertices, faces = arrays['vertices'], arrays['faces']
    roles, parents = arrays['source_component_ids'], arrays['source_parent_face_ids']
    owners = np.full(len(faces), -1, dtype=np.int64)
    for role, label in [(0, 'left'), (1, 'right')]:
        station = np.full(len(meshes[role].vertices), -1, dtype=np.int64)
        rings = fields[label + '_rings']
        station[rings] = np.arange(len(rings))[:, None]
        assert np.all(station >= 0)
        source_face_station = station[meshes[role].faces].min(1)
        selected = np.flatnonzero(roles == role)
        local_parent = parents[selected] - nf[role]
        assert local_parent.min() >= 0 and local_parent.max() < len(meshes[role].faces)
        owners[selected] = role * 2 + (source_face_station[local_parent] >= 13)
    owners[roles == 2] = 4
    assert set(owners) == set(range(5)) and np.all(parents >= 0)
    normals, uvs = arrays['corner_normal_words'].copy(), arrays['corner_uvs'].copy()
    shared = []
    for vertex in np.unique(faces):
        where = np.argwhere(faces == vertex)
        if len(set(owners[where[:, 0]])) < 2:
            continue
        word = normals[where[:, 0], where[:, 1]].astype(float).sum(0)
        assert np.linalg.norm(word) > 0
        word = np.rint(word / np.linalg.norm(word) * 32767).astype(np.int16)
        normals[where[:, 0], where[:, 1]] = word
        shared.append(vertex)
    origin = np.asarray(plan['storage_origin_m'])
    stored = (vertices - origin).astype('<f4')
    actual = stored.astype(float) + origin
    merged = Mesh(actual, faces, frame)
    compiled = audit(merged, require_closed=False, expected_components=1,
                     max_pairs=600000, max_tests=2000000)
    write(output / 'compiled-merged-audit.json', compiled)
    assert compiled['passed'], compiled['failures']
    assert len(boundary_loops(faces)) == 4 and topology(merged)['components'][0]['euler'] == -2
    boundary_keys = [k for k in arrays if k.startswith('external_boundary_') and k.endswith('_vertex_ids')]
    assert len(boundary_keys) == 4, 'Constructor did not emit the four declared exterior profiles'
    source_mouths = {frozenset(map(tuple, m.vertices[loop]))
                     for m in meshes[:2] for loop in boundary_loops(m.faces)}
    actual_mouths = {frozenset(map(tuple, actual[loop])) for loop in boundary_loops(faces)}
    assert len(source_mouths) == 4 and source_mouths == actual_mouths, 'Actual decoded exterior mouths differ from all four literal source mouths'
    for key in boundary_keys:
        port = arrays[key]
        assert np.array_equal(actual[port], vertices[port]), 'Fixed external mouth lost Float32 identity'
    blood = classifier(merged)
    blood_audit = audit(blood, require_closed=True, expected_components=1, expected_genus=[0],
                        max_pairs=600000, max_tests=2000000)
    write(output / 'diagnostic-blood-classifier-audit.json', blood_audit)
    assert blood_audit['passed'], blood_audit['failures']
    probes = [arrays['communication_probe_points_m']]
    probes.extend(fields[label + '_centres'][1:-1] for label in ['left', 'right'])
    probes.append(fields['bridge_centres'])
    require_inside(blood, np.concatenate(probes), 'A declared centre-course probe is outside the merged cast')
    blood.save(output / 'diagnostic-only-blood-classifier.npz')
    after = copy.deepcopy(atlas)
    data = bytearray()
    chunk = len(atlas['chunks'])
    def append(array):
        array = np.ascontiguousarray(array)
        while len(data) % array.dtype.itemsize:
            data.append(0)
        offset = len(data)
        data.extend(array.tobytes())
        return offset
    for owner, pid in enumerate(NEW):
        group = np.flatnonzero(owners == owner)
        cp, cn, cu = stored[faces[group]].reshape(-1, 3), normals[group].reshape(-1, 3), uvs[group].reshape(-1, 2)
        dtype = np.dtype([('p', '<f4', (3,)), ('n', '<i2', (3,)), ('u', '<f4', (2,))])
        keys = np.empty(len(cp), dtype)
        keys['p'], keys['n'], keys['u'] = cp, cn, cu
        _, first, indices = np.unique(keys, return_index=True, return_inverse=True)
        q, n, uv, ff = cp[first], cn[first], cu[first], indices.reshape(-1, 3).astype('<u4')
        oldid = OLD[owner // 2] if owner < 4 else OLD[2]
        record = copy.deepcopy(records[oldid])
        side = 'left' if owner < 2 else 'right'
        region = 'Precommunicating' if owner % 2 == 0 else 'Postcommunicating'
        name = f'{region} part of anterior cerebral artery ({side})' if owner < 4 else 'Anterior communicating artery'
        for key in ['sourceOffset', 'geometryScale', 'normalSigns', 'geometryInstance', 'suppressed', 'defaultHidden']:
            record.pop(key, None)
        record.update(id=pid, conceptId=pid, name=name, ta98Term=TERMS[owner],
                      positions=append(q), normals=append(n), uvs=append(uv), indices=append(ff),
                      chunk=chunk, geometryOrigin=origin.tolist(), vertexCount=len(q), indexCount=int(ff.size),
                      bounds=[actual[faces[group]].min((0, 1)).tolist(), actual[faces[group]].max((0, 1)).tolist()],
                      sourceOuterSurfaceOnly=True, coverageLimitation=LIMIT,
                      reviewStatus='local macroscopic exterior; source regional continuation interfaces qualified',
                      nativeACACommunicationAssembly=dict(schema='ACA-ACom-five-piece/v1',
                          sourceCatalogueSha256=sha(CAT), role=owner, finiteWall=False,
                          sourceCalibreHeld=True, sourceEndProfilesLiteral=True))
        if owner == 4:
            record['provenance']['detail'] += ' Communicating bridge independently constructed from source course and calibre evidence; original closed ellipsoid retained only as hidden reference.'
            record['uvMethod'] = 'Procedural longitudinal/circumferential coordinates'
        after['parts'].append(record)
        after['concepts'].append(dict(id=pid, name=name, elements=[pid], partOf=oldid,
             parentId=oldid, anatomicalParentId=oldid, ta98Term=TERMS[owner], ta98Kind='exact',
             coverageLimitation=LIMIT))
        after['ta98TermBindings'].append(dict(conceptId=pid, partIds=[pid], term=TERMS[owner], kind='exact'))
        np.savez_compressed(output / f'owned-{owner}.npz', vertices=actual, faces=faces[group],
                            merged_face_ids=group, corner_normal_words=normals[group], corner_uvs=uvs[group])
    for index, oldid in enumerate(OLD):
        children = NEW[index * 2:index * 2 + 2] if index < 2 else NEW[4:]
        old = next(p for p in after['parts'] if p['id'] == oldid)
        old.update(suppressed=True, defaultHidden=True, replacementParts=children)
        if oldid not in after['defaultHidden']:
            after['defaultHidden'].append(oldid)
        concept = next(c for c in after['concepts'] if c['id'] == oldid)
        concept.update(elements=children, componentElements=children,
                       ta98Term='A12.2.07.022' if index < 2 else 'A12.2.07.029',
                       ta98Kind='exact', coverageLimitation=LIMIT)
    unrelated = [p for p in after['parts'] if p['id'] not in set(OLD + NEW)]
    assert unrelated == [p for p in atlas['parts'] if p['id'] not in OLD]
    bank = output / 'aca-communicating-exterior.bin'
    bank.write_bytes(data)
    after['chunks'].append(dict(url='/local-models/' + str(bank.relative_to(ROOT / '.local-models')), bytes=len(data)))
    after['partsCount'] = len(after['parts'])
    after['triangles'] = sum(p['indexCount'] // 3 for p in after['parts'])
    write(output / 'proposed-atlas.json', after)
    for owner, pid in enumerate(NEW):
        record = next(p for p in after['parts'] if p['id'] == pid)
        AtlasReader().read_part(after, record, frame=frame).save(output / f'packed-{owner}.npz')
    np.savez_compressed(output / 'compiled-shared-geometry.npz', vertices=actual, faces=faces,
                        face_owner=owners, corner_normal_words=normals, corner_uvs=uvs,
                        shared_interface_control_ids=np.asarray(shared))
    proof = dict(schema='ACA-ACom-five-piece-packing-proof/v1', complete=True, passed=True,
                 installed=False, anatomical_acceptance=False, source_catalogue_sha256=sha(CAT),
                 proposed_catalogue_sha256=sha(output / 'proposed-atlas.json'), bank_sha256=sha(bank),
                 bank_bytes=len(data), changed_parts=OLD + NEW, new_TA98_terms=['A12.2.07.023', 'A12.2.07.033'],
                 unrelated_part_records_literal=len(unrelated), five_active_owned_pieces=True,
                 shared_interface_controls=len(shared), original_outer_mouths=4,
                 every_face_owned_once=True, temporary_classifier_faces_visible=0,
                 bridge_winding_reversed=reversed_bridge, finite_wall=False,
                 elapsed_seconds=time.monotonic() - started)
    write(output / 'packing-proof.json', proof)
    write(output / 'construction-state.json', proof)
    print(json.dumps(proof))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', required=True, type=Path)
    parser.add_argument('--approval', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    run(args.candidate, args.approval, args.output)
