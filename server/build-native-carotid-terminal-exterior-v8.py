"""Guarded source-bound ICA/MCA/ACA/ACom joint exterior candidate.

No work runs on import. This script writes a candidate only. It cannot install,
modify the current catalogue, move a fixed organ or perform a course solve.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
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

SPEC = importlib.util.spec_from_file_location('aca_builder_helpers',
            ROOT / 'server/build-aca-communicating-exterior-v1.py')
HELPERS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(HELPERS)
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1'
COHORT = WORK / 'native-carotid-terminal-cohort-v7b'
OLD = HELPERS.OLD
NEW = HELPERS.NEW
HOSTS = ['TA98PART:ZA:Internal carotid artery.l:A12.2.06.015',
         'ZA:Middle cerebral artery (M1-segment).l',
         'TA98PART:ZA:Internal carotid artery.r:A12.2.06.015',
         'ZA:Middle cerebral artery (M1-segment).r']
ACTIVE = HOSTS + NEW
TERMS = ['A12.2.06.015', 'A12.2.07.047', 'A12.2.06.015', 'A12.2.07.047'] + HELPERS.TERMS
LIMIT = ('Source-bound macroscopic joint carotid-terminal exterior. The source '
         'MCA-labelled proximal course is retained as terminal ICA tissue; '
         'station 17 is a provisional geometric ownership boundary, not a '
         'universal clinical M1 level. Shared openings and six exterior source '
         'mouths are preserved. No finite wall, microscopic or complete '
         'perfusion claim. Fixed context and source-interface anatomy require '
         'separate review.')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=lambda x:
        x.tolist() if hasattr(x, 'tolist') else str(x)) + '\n')


def check_source(reader, atlas, record, frozen):
    actual = reader.read_part(atlas, record, frame=frozen.frame)
    assert actual.geometry_sha == frozen.geometry_sha, 'Frozen source geometry differs from live source'
    assert actual.attributes_sha == frozen.attributes_sha, 'Frozen source fields differ from live source'


def sew_carrier(ica, mca, side_plan):
    """Sew one canonical geometric ring while keeping both source UV charts."""
    correspondence = side_plan['correspondence']
    ii = np.asarray(correspondence['ica_distal_control_ids'], int)
    mi = np.asarray(correspondence['mca_control_ids_in_ica_order'], int)
    assert len(ii) == len(mi) == 12 and len(set(ii)) == len(set(mi)) == 12
    iv, inn = ica.vertices.copy(), ica.attributes['normals'].copy()
    iv[ii] = mca.vertices[mi]
    inn[ii] = mca.attributes['normals'][mi]
    outside = np.setdiff1d(np.arange(len(iv)), ii)
    assert np.array_equal(iv[outside], ica.vertices[outside])
    assert np.array_equal(inn[outside], ica.attributes['normals'][outside])
    original_corner_uvs = np.concatenate([ica.attributes['uvs'][ica.faces],
                                          mca.attributes['uvs'][mca.faces]])
    all_vertices = np.concatenate([iv, mca.vertices])
    all_faces = np.concatenate([ica.faces, mca.faces + len(iv)])
    all_normals = np.concatenate([inn, mca.attributes['normals']])
    all_uvs = np.concatenate([ica.attributes['uvs'], mca.attributes['uvs']])
    vertices, first, inverse = np.unique(all_vertices, axis=0, return_index=True,
                                        return_inverse=True)
    faces = inverse[all_faces]
    # Only the twelve explicitly corresponding ring controls may be welded.
    assert len(all_vertices) - len(vertices) == 12
    assert np.array_equal(inverse[ii], inverse[len(iv) + mi])
    for a, b in zip(ii, mi):
        assert np.array_equal(all_normals[a], all_normals[len(iv) + b])
    mesh = Mesh(vertices, faces, ica.frame,
                dict(normals=all_normals[first], uvs=all_uvs[first]),
                dict(only_ICA_distal_12_positions_and_normal_words_changed=True,
                     MCA_all_positions_normal_words_UVs_literal=True,
                     source_corner_UV_charts_preserved_separately=True))
    assert len(boundary_loops(mesh.faces)) == 2
    return mesh, original_corner_uvs, dict(
        ICA_changed_control_ids=ii, MCA_canonical_control_ids=mi,
        canonical_positions_m=mca.vertices[mi], canonical_normal_words=mca.attributes['normals'][mi],
        ICA_other_controls_literal=True, all_MCA_source_controls_literal=True,
        shared_geometric_control_count=12, source_UV_charts_separate=True,
        ICA_face_count=len(ica.faces), MCA_face_count=len(mca.faces))


def extend_aca_root(aca, rings, root_profile, side_plan):
    """Add only the explicitly planned short frustum; original ACA is literal."""
    profile = np.asarray(root_profile['profile'], float)
    root = rings[0]
    assert profile.shape == (12, 3) and len(root) == 12
    assert np.array_equal(root_profile['anchor'], np.asarray(side_plan['internal_profile']['anchor_m']))
    centre = profile.mean(0)
    pdir = profile - centre
    rdir = aca.vertices[root] - aca.vertices[root].mean(0)
    pdir /= np.linalg.norm(pdir, axis=1)[:, None]
    rdir /= np.linalg.norm(rdir, axis=1)[:, None]
    assert np.einsum('ij,ij->i', pdir, rdir).min() > .95, 'Root profile phase does not match source ACA ring'
    vertices = np.concatenate([aca.vertices, profile])
    small = np.arange(len(aca.vertices), len(vertices))
    frustum = []
    for j in range(12):
        k = (j + 1) % 12
        frustum.extend([[small[j], small[k], root[j]],
                        [small[k], root[k], root[j]]])
    frustum = np.asarray(frustum, int)
    # The new root edge must run opposite the original ACA boundary edge.
    directed = {(int(loop[j]), int(loop[(j + 1) % len(loop)]))
                for loop in boundary_loops(aca.faces) for j in range(len(loop))}
    if (int(root[0]), int(root[1])) not in directed:
        assert (int(root[1]), int(root[0])) in directed
        frustum = frustum[:, ::-1]
    triangles = vertices[frustum]
    face_normals = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    small_normals = np.zeros((12, 3))
    for face, normal in zip(frustum, face_normals):
        for index in face:
            if index >= len(aca.vertices):
                small_normals[index - len(aca.vertices)] += normal
    assert np.all(np.linalg.norm(small_normals, axis=1) > 0)
    small_normals /= np.linalg.norm(small_normals, axis=1)[:, None]
    normal_words = np.concatenate([aca.attributes['normals'],
                                   np.rint(small_normals * 32767).astype(np.int16)])
    # Only new controls need a new chart. The original ACA chart is literal.
    small_uvs = np.column_stack([np.arange(12) / 12, np.zeros(12)]).astype(np.float32)
    uvs = np.concatenate([aca.attributes['uvs'], small_uvs])
    faces = np.concatenate([aca.faces, frustum])
    mesh = Mesh(vertices, faces, aca.frame, dict(normals=normal_words, uvs=uvs),
                dict(original_ACA_positions_faces_normal_words_UVs_literal=True,
                     root_frustum_source_faces=list(range(len(aca.faces), len(faces)))))
    assert len(boundary_loops(mesh.faces)) == 2
    return mesh, small, dict(source_ACA_faces=len(aca.faces), root_frustum_faces=len(frustum),
                              explicit_internal_profile_controls=small,
                              source_original_controls_literal=True)


def run(candidate, approval_file, interface_file, output):
    started = time.monotonic()
    candidate, interface_file, output = candidate.resolve(), interface_file.resolve(), output.resolve()
    manifest_path = candidate / 'candidate-manifest.json'
    gates_path = candidate / 'actual-geometry-gates.json'
    sections_path = candidate / 'section-isometry-proof.json'
    fields_path = candidate / 'source-and-candidate-fields.npz'
    manifest, gates = json.loads(manifest_path.read_text()), json.loads(gates_path.read_text())
    approval, interfaces = json.loads(approval_file.read_text()), json.loads(interface_file.read_text())
    assert gates.get('complete') and gates.get('passed'), 'Complete course geometry gates required'
    for key, path in [('candidate_manifest_sha256', manifest_path),
                      ('candidate_geometry_gates_sha256', gates_path),
                      ('candidate_section_proof_sha256', sections_path),
                      ('candidate_fields_sha256', fields_path),
                      ('joint_interface_plan_sha256', interface_file)]:
        assert approval[key] == sha(path), 'Changed approved input: ' + key
    for key in ['course_passed', 'root_surface_course_review_passed',
                'root_true_section_review_passed', 'fixed_context_review_passed',
                'joint_carotid_terminal_interface_approved']:
        assert approval.get(key) is True, 'Missing root approval: ' + key
    assert approval.get('context_evidence')
    for evidence in approval['context_evidence']:
        assert sha(evidence['path']) == evidence['sha256']
    source_hash = sha(CAT)
    assert source_hash == manifest['source_catalogue_sha256'] == interfaces['source_catalogue_sha256'] == approval['source_catalogue_sha256']
    assert interfaces.get('complete')
    inventory_path = Path(interfaces['inventory_path'])
    inventory = json.loads(inventory_path.read_text())
    assert inventory['term_count'] == len(inventory['terms']) == 129
    assert approval['joint_inventory_sha256'] == sha(inventory_path)
    assert inventory.get('source_catalogue_sha256', inventory.get('current_catalogue_sha256')) == source_hash
    assert not output.exists(), 'Use a fresh candidate directory'
    output.mkdir(parents=True)
    write(output / 'construction-state.json', dict(complete=False, passed=False, installed=False))
    (output / 'authority-inventory.json').write_bytes(inventory_path.read_bytes())
    atlas, reader = json.loads(CAT.read_text()), AtlasReader()
    records = {p['id']: p for p in atlas['parts']}
    assert all(pid in records and not records[pid].get('suppressed') for pid in OLD + HOSTS)
    assert not set(NEW) & set(records)
    fields = np.load(fields_path, allow_pickle=False)
    raw_aca = [Mesh.load(candidate / 'aca-left.npz'), Mesh.load(candidate / 'aca-right.npz')]
    bridge = Mesh.load(candidate / 'acom-open-bridge.npz')
    frame = raw_aca[0].frame
    bridge, reverse_bridge = HELPERS.orient_bridge(bridge, fields['bridge_fields'], fields['bridge_centres'])
    assert bridge.frame == frame
    carriers, acas, corner_uvs, seam_evidence, root_evidence = [], [], [], [], []
    local_owner, original_mouths, anchor_by_controls = [], [], []
    side_plans = {p['side']: p for p in interfaces['sides']}
    for role, side, label in [(0, 'l', 'left'), (1, 'r', 'right')]:
        ica, mca = Mesh.load(COHORT / f'ica-{side}.npz'), Mesh.load(COHORT / f'mca-{side}.npz')
        assert ica.frame == mca.frame == frame
        check_source(reader, atlas, records[HOSTS[role * 2]], ica)
        check_source(reader, atlas, records[HOSTS[role * 2 + 1]], mca)
        carrier, cu, seam = sew_carrier(ica, mca, side_plans[side])
        carrier_audit = audit(carrier, require_closed=False, expected_components=1, max_pairs=600000, max_tests=2000000)
        write(output / f'{side}-sewn-carrier-audit.json', carrier_audit)
        assert carrier_audit['passed'], carrier_audit['failures']
        carriers.append(carrier); corner_uvs.append(cu); seam_evidence.append(seam)
        stations = np.load(COHORT / f'mca-{side}-stations.npz', allow_pickle=False)
        mca_station = np.full(len(mca.vertices), -1, dtype=np.int64)
        mca_station[stations['rings']] = np.arange(len(stations['rings']))[:, None]
        assert np.all(mca_station >= 0)
        station_owner = mca_station[mca.faces].min(1)
        carrier_owner = np.r_[np.full(len(ica.faces), role * 2),
                               role * 2 + (station_owner >= 17)]
        local_owner.append(carrier_owner)
        source_aca = Mesh.load(COHORT / f'aca-{side}.npz')
        rings = fields[label + '_rings']
        retained = np.r_[np.arange(3), np.arange(23, len(rings))]
        assert np.array_equal(raw_aca[role].faces, source_aca.faces)
        assert np.array_equal(raw_aca[role].vertices[rings[retained]], source_aca.vertices[rings[retained]])
        assert np.array_equal(raw_aca[role].attributes['normals'][rings[retained]], source_aca.attributes['normals'][rings[retained]])
        assert np.array_equal(raw_aca[role].attributes['uvs'], source_aca.attributes['uvs'])
        planned_root = np.load(COHORT / f'{side}-planned-internal-root-profile.npz', allow_pickle=False)
        aca, small, root = extend_aca_root(raw_aca[role], rings, planned_root, side_plans[side])
        acas.append(aca); root_evidence.append(root)
        anchor_by_controls.append((2 + role, small, np.asarray(side_plans[side]['internal_profile']['anchor_m'])))
        root_cone = np.r_[aca.vertices[small], np.asarray(side_plans[side]['internal_profile']['anchor_m'])[None]]
        HELPERS.require_inside(HELPERS.classifier(carrier), root_cone, 'ACA root closure is not buried in native terminal carrier')
        original_mouths.extend([ica.vertices[planned_root['ica_proximal_ring']],
                                mca.vertices[stations['rings'][-1]],
                                source_aca.vertices[rings[-1]]])
    meshes = carriers + acas + [bridge]
    # Preserve original source ACA faces and assign the added frustum to A1.
    for role, label in [(0, 'left'), (1, 'right')]:
        rings = fields[label + '_rings']
        station = np.full(len(raw_aca[role].vertices), -1, int)
        station[rings] = np.arange(len(rings))[:, None]
        source_station = station[raw_aca[role].faces].min(1)
        local_owner.append(np.r_[4 + role * 2 + (source_station >= 13),
                                 np.full(len(acas[role].faces) - len(raw_aca[role].faces), 4 + role * 2)])
        corner_uvs.append(acas[role].attributes['uvs'][acas[role].faces])
    local_owner.append(np.full(len(bridge.faces), 8))
    corner_uvs.append(bridge.attributes['uvs'][bridge.faces])
    nv, nf = np.cumsum([0] + [len(m.vertices) for m in meshes]), np.cumsum([0] + [len(m.faces) for m in meshes])
    composite = Mesh(np.concatenate([m.vertices for m in meshes]),
        np.concatenate([m.faces + nv[i] for i, m in enumerate(meshes)]), frame,
        {k: np.concatenate([m.attributes[k] for m in meshes]) for k in ['normals', 'uvs']})
    source_corner_uvs, source_owners = np.concatenate(corner_uvs), np.concatenate(local_owner)
    mouths, anchors = [], {}
    bcentres = fields['bridge_centres']
    bdir = bcentres[-1] - bcentres[0]; bdir /= np.linalg.norm(bdir)
    radius = float(manifest['bridge_radius_m'])
    for endpoint, side in [(0, 0), (-1, 1)]:
        profile = fields['bridge_fields'][endpoint]
        anchor = bcentres[endpoint] + (-1 if endpoint == 0 else 1) * bdir * radius / 2
        HELPERS.require_inside(HELPERS.classifier(acas[side]), np.r_[profile, anchor[None], (profile + anchor) / 2], 'ACom bridge closure is not buried in the corrected ACA')
        ids = np.arange(16) if endpoint == 0 else np.arange(len(bridge.vertices) - 16, len(bridge.vertices))
        anchor_by_controls.append((4, ids, anchor))
    for role, mesh in enumerate(meshes):
        for loop in boundary_loops(mesh.faces):
            internal = [(ids, point) for r, ids, point in anchor_by_controls
                        if r == role and set(map(int, ids)) == set(map(int, loop))]
            assert len(internal) <= 1
            ids = loop + nv[role]
            mouths.append(dict(positions_m=composite.vertices[ids], source_control_ids=ids,
                               kind='internal communicating root' if internal else 'original'))
            if internal:
                anchors[len(mouths) - 1] = internal[0][1]
    assert len(mouths) == 10 and len(anchors) == 4
    plan = dict(schema='native-carotid-terminal-joint-exterior/v8',
        whole_organ='Both cerebral ICA terminal courses, native M1 courses and corrected ACA pre/post regions connected by one ACom bridge.',
        source_catalogue_sha256=source_hash, candidate_manifest_sha256=sha(manifest_path),
        joint_interface_plan=str(interface_file), joint_interface_plan_sha256=sha(interface_file),
        approval_sha256=sha(approval_file), changed_existing_parts=HOSTS + OLD, new_parts=NEW,
        active_parts=ACTIVE, component_roles=['left ICA/MCA carrier', 'right ICA/MCA carrier', 'left ACA/root frustum', 'right ACA/root frustum', 'ACom bridge'],
        original_exterior_mouths=6, internal_roots=4, mouths=mouths, anchors=anchors,
        original_exterior_profiles=original_mouths, ICA_MCA_seams=seam_evidence,
        ACA_root_frusta=root_evidence, every_output_face_owned_once=True,
        adjoining_interfaces=[
            'ICA distal 12-control profile adopts the literal canonical MCA proximal positions/normal words. All other ICA controls and all MCA source controls remain literal.',
            'Native ICA/MCA source corner UV charts remain separate across the welded geometric ring.',
            'The proximal source MCA course through strip 16 is terminal ICA ownership; station 17 is a provisional geometric boundary, not a universal clinical M1 landmark.',
            'Each retained ACA root is extended to the explicitly verified buried 0.4 mm profile through one short frustum. The closure cone is diagnostic only; visible frustum fragments need actual surface and cut review.',
            'The bilateral ACA and ACom share actual open side ostia; no internal disk is emitted.',
            'Six source continuations remain open and literal: ICA proximal, M1 distal, ACA distal on each side. Fixed brain, nerve, skull and foreign-vessel interfaces need separate actual context review.'
        ], fixed_anatomy='All sources outside the seven declared source records, including basilar artery, skeleton, brain and nerves.',
        storage_origin_m=[0., 1.71875, .03125], source_corner_UV_replay=True,
        finite_wall=False, physiological_patency_claim=False, installed=False,
        anatomical_acceptance=False)
    write(output / 'assembly-plan-before-union.json', plan)
    composite.save(output / 'composite-input.npz')
    np.savez_compressed(output / 'source-face-corner-fields.npz', corner_uvs=source_corner_uvs,
                        source_face_owner=source_owners, source_component_face_ranges=nf)
    arrays, native = construct(composite, [np.arange(nf[i], nf[i + 1]) for i in range(5)],
        mouths, anchors, wall_seconds=35, max_pairs=600000, max_tests=2000000)
    np.savez_compressed(output / 'union-source-provenance.npz', **arrays)
    write(output / 'native-tree-proof.json', native)
    assert native['passed'], native['failures']
    assert native['derived_fragments'] == 0, 'A temporary closure cone remains visible'
    v, f, parents = arrays['vertices'], arrays['faces'], arrays['source_parent_face_ids']
    assert np.all(parents >= 0)
    owners = source_owners[parents]
    assert set(owners) == set(range(9))
    normals = arrays['corner_normal_words'].copy()
    uv = np.empty((len(f), 3, 2), np.float32)
    for i, face in enumerate(f):
        parent = int(parents[i]); chart = composite.triangles[parent]
        uv[i] = (arrays['corner_source_barycentric'][i] @ source_corner_uvs[parent]).astype(np.float32)
        for corner, p in enumerate(v[face]):
            exact = np.flatnonzero(np.all(chart == p, axis=1))
            if len(exact):
                uv[i, corner] = source_corner_uvs[parent, exact[0]]
    shared = []
    for vi in np.unique(f):
        where = np.argwhere(f == vi)
        if len(set(owners[where[:, 0]])) < 2:
            continue
        words = normals[where[:, 0], where[:, 1]]
        if len(np.unique(words, axis=0)) == 1:
            word = words[0]
        else:
            word = words.astype(float).sum(0)
            assert np.linalg.norm(word) > 0
            word = np.rint(word / np.linalg.norm(word) * 32767).astype(np.int16)
        normals[where[:, 0], where[:, 1]] = word
        shared.append(int(vi))
    origin = np.asarray(plan['storage_origin_m'])
    stored, actual = (v - origin).astype('<f4'), None
    actual = stored.astype(float) + origin
    compiled = Mesh(actual, f, frame)
    compiled_audit = audit(compiled, require_closed=False, expected_components=1,
                           max_pairs=600000, max_tests=2000000)
    write(output / 'compiled-merged-audit.json', compiled_audit)
    assert compiled_audit['passed'], compiled_audit['failures']
    assert topology(compiled)['components'][0]['euler'] == -4
    loops = boundary_loops(f)
    expected_ports = {frozenset(map(tuple, p)) for p in original_mouths}
    actual_ports = {frozenset(map(tuple, actual[loop])) for loop in loops}
    assert len(expected_ports) == len(loops) == 6 and actual_ports == expected_ports
    boundary_keys = [key for key in arrays if key.startswith('external_boundary_') and key.endswith('_vertex_ids')]
    assert len(boundary_keys) == 6
    for key in boundary_keys:
        assert np.array_equal(actual[arrays[key]], v[arrays[key]])
    blood = HELPERS.classifier(compiled)
    blood_audit = audit(blood, require_closed=True, expected_components=1, expected_genus=[0], max_pairs=600000, max_tests=2000000)
    write(output / 'diagnostic-blood-classifier-audit.json', blood_audit)
    assert blood_audit['passed'], blood_audit['failures']
    HELPERS.require_inside(blood, arrays['communication_probe_points_m'], 'Joint cast communication probe is outside')
    blood.save(output / 'diagnostic-only-blood-classifier.npz')
    after, data = copy.deepcopy(atlas), bytearray()
    chunk = len(atlas['chunks'])
    def append(array):
        array = np.ascontiguousarray(array)
        while len(data) % array.dtype.itemsize:
            data.append(0)
        offset = len(data); data.extend(array.tobytes()); return offset
    names = [records[p]['name'] for p in HOSTS] + [
        'Precommunicating part of anterior cerebral artery (left)',
        'Postcommunicating part of anterior cerebral artery (left)',
        'Precommunicating part of anterior cerebral artery (right)',
        'Postcommunicating part of anterior cerebral artery (right)',
        'Anterior communicating artery']
    for owner, pid in enumerate(ACTIVE):
        fi = np.flatnonzero(owners == owner)
        cp, cn, cu = stored[f[fi]].reshape(-1, 3), normals[fi].reshape(-1, 3), uv[fi].reshape(-1, 2)
        dtype = np.dtype([('p', '<f4', (3,)), ('n', '<i2', (3,)), ('u', '<f4', (2,))])
        keys = np.empty(len(cp), dtype); keys['p'], keys['n'], keys['u'] = cp, cn, cu
        _, first, indices = np.unique(keys, return_index=True, return_inverse=True)
        q, n, u, ff = cp[first], cn[first], cu[first], indices.reshape(-1, 3).astype('<u4')
        parent = None if owner < 4 else OLD[(owner - 4) // 2] if owner < 8 else OLD[2]
        record = next(p for p in after['parts'] if p['id'] == pid) if owner < 4 else copy.deepcopy(records[parent])
        for key in ['sourceOffset', 'geometryScale', 'normalSigns', 'geometryInstance', 'suppressed', 'defaultHidden']:
            record.pop(key, None)
        record.update(id=pid, conceptId=pid, name=names[owner], ta98Term=TERMS[owner],
            positions=append(q), normals=append(n), uvs=append(u), indices=append(ff),
            chunk=chunk, geometryOrigin=origin.tolist(), vertexCount=len(q), indexCount=int(ff.size),
            bounds=[actual[f[fi]].min((0, 1)).tolist(), actual[f[fi]].max((0, 1)).tolist()],
            sourceOuterSurfaceOnly=True, coverageLimitation=LIMIT,
            nativeCarotidTerminalAssembly=dict(schema='joint-terminal-carotid-v8', role=owner,
                sourceCatalogueSha256=source_hash, MCAownershipBoundaryStation=17,
                boundaryProvisional=True, finiteWall=False))
        if owner >= 4:
            after['parts'].append(record)
            after['concepts'].append(dict(id=pid, name=names[owner], elements=[pid], partOf=parent,
                parentId=parent, anatomicalParentId=parent, ta98Term=TERMS[owner],
                ta98Kind='exact', coverageLimitation=LIMIT))
        elif owner in [1, 3]:
            concept = next(c for c in after['concepts'] if c['id'] == pid)
            concept.update(ta98Term=TERMS[owner], ta98Kind='exact', coverageLimitation=LIMIT)
        binding = dict(conceptId=pid, partIds=[pid], term=TERMS[owner], kind='exact')
        previous = [b for b in after['ta98TermBindings'] if b.get('conceptId') == pid]
        assert not previous or all(b['term'] == TERMS[owner] for b in previous)
        if not previous:
            after['ta98TermBindings'].append(binding)
        np.savez_compressed(output / f'owned-{owner}.npz', vertices=actual, faces=f[fi],
                            merged_face_ids=fi, corner_normal_words=normals[fi], corner_uvs=uv[fi])
    for i, pid in enumerate(OLD):
        children = NEW[i * 2:i * 2 + 2] if i < 2 else NEW[4:]
        record = next(p for p in after['parts'] if p['id'] == pid)
        record.update(suppressed=True, defaultHidden=True, replacementParts=children)
        if pid not in after['defaultHidden']:
            after['defaultHidden'].append(pid)
        concept = next(c for c in after['concepts'] if c['id'] == pid)
        concept.update(elements=children, componentElements=children,
            ta98Term='A12.2.07.022' if i < 2 else 'A12.2.07.029', ta98Kind='exact', coverageLimitation=LIMIT)
    outside = [p for p in after['parts'] if p['id'] not in set(OLD + ACTIVE)]
    assert outside == [p for p in atlas['parts'] if p['id'] not in set(OLD + HOSTS)]
    bank = output / 'native-carotid-terminal-exterior.bin'; bank.write_bytes(data)
    after['chunks'].append(dict(url='/local-models/' + str(bank.relative_to(ROOT / '.local-models')), bytes=len(data)))
    after['partsCount'] = len(after['parts'])
    after['triangles'] = sum(p['indexCount'] // 3 for p in after['parts'])
    write(output / 'proposed-atlas.json', after)
    for owner, pid in enumerate(ACTIVE):
        record = next(p for p in after['parts'] if p['id'] == pid)
        reader.read_part(after, record, frame=frame).save(output / f'packed-{owner}.npz')
    np.savez_compressed(output / 'compiled-shared-geometry.npz', vertices=actual, faces=f,
        face_owner=owners, corner_normal_words=normals, corner_uvs=uv,
        shared_interface_control_ids=np.asarray(shared), original_outer_profile_count=np.asarray(6))
    report = dict(schema='joint-terminal-carotid-packing-proof/v8', complete=True, passed=True,
        source_catalogue_sha256=source_hash, proposed_catalogue_sha256=sha(output / 'proposed-atlas.json'),
        bank_sha256=sha(bank), bank_bytes=len(data), changed_parts=HOSTS + OLD + NEW,
        changed_existing_parts=HOSTS + OLD, new_parts=NEW, active_parts=ACTIVE,
        explicitly_bound_terms=sorted(set(TERMS)), unchanged_existing_parts=len(outside),
        nine_active_owned_parts=True, every_face_owned_once=True, shared_controls=len(shared),
        source_corner_UVs_replayed=True, literal_original_exterior_mouths=6,
        temporary_classifier_faces_visible=0, derived_closure_faces_visible=0,
        MCA_station17_boundary_provisional=True, source_outer_surface_only=True,
        finite_wall=False, full_perfusion_claim=False, bridge_winding_reversed=reverse_bridge,
        installed=False, anatomical_acceptance=False, elapsed_seconds=time.monotonic() - started)
    write(output / 'packing-proof.json', report)
    write(output / 'construction-state.json', report)
    print(json.dumps(report))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--approval', type=Path, required=True)
    parser.add_argument('--interface-plan', type=Path, default=COHORT / 'joint-interface-plan-v8.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(args.candidate, args.approval, args.interface_plan, args.output)
