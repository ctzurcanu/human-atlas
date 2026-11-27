"""Literal mandibular defect localization only; not an occupied classifier."""
from pathlib import Path
from collections import defaultdict, Counter
import importlib.util
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts, patch_components


def main():
    began = time.monotonic()
    work = ROOT / '.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17'
    out = work / 'mandible-source-defect-localization-v32'
    assert not out.exists()
    out.mkdir()
    binding_path = work / 'source-binding-v19/evidence.json'
    binding = json.loads(binding_path.read_bytes())
    source = ROOT / binding['snapshot']
    expected = 'd73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
    assert digest_file(source) == expected
    assert digest_file(ROOT / binding['current_catalogue']) == expected
    atlas = json.loads(source.read_bytes())
    reader = AtlasReader()
    wall = Mesh.load(work / 'finite-wall.npz')
    bank_path = work / 'surfaces-and-ancestry.npz'
    bank = np.load(bank_path)
    names = json.loads((work / 'ports.json').read_bytes())['component_names']
    spec = importlib.util.spec_from_file_location('mandible_exact_surface_normalization_v32', ROOT / 'server/normalize-facial-source-classifiers-v26.py')
    normal = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(normal)
    selected = [p for p in atlas['parts'] if p.get('sourceId') == 'Mandible' and not p.get('suppressed')]
    assert [p['id'] for p in selected] == ['ZA:Mandible:Bone', 'ZA:Mandible:Cartilage', 'ZA:Mandible:Suture-1']
    vv, ff, face_part, face_local, members = [], [], [], [], []
    offset = 0
    for pi, p in enumerate(selected):
        mesh = reader.read_part(atlas, p, frame=wall.frame, catalogue_path=source, catalogue_sha=expected)
        v0, inv0 = np.unique(mesh.vertices, axis=0, return_inverse=True)
        c0, n0 = normal.exact_components(v0, inv0[mesh.faces])
        members.append(dict(part_id=p['id'], source_object_id=p['sourceId'], material=p['material'], system=p['system'], geometry_sha256=mesh.geometry_sha, source_bank_sha256=mesh.source['bank_sha256'], original_topology=topology(mesh), exact_surface_edge_components=[dict(face_count=len(f), vertex_count=len(u), topology=topology(Mesh(v0[u],f,wall.frame))) for u,f,_ in c0], normalization=n0))
        vv.append(mesh.vertices)
        ff.append(mesh.faces + offset)
        face_part.extend([pi] * len(mesh.faces))
        face_local.extend(range(len(mesh.faces)))
        offset += len(mesh.vertices)
    v, inv = np.unique(np.concatenate(vv), axis=0, return_inverse=True)
    rawf = inv[np.concatenate(ff)]
    components, normalization = normal.exact_components(v, rawf)
    assert len(components) == 1
    used, f, original_faces = components[0]
    vertices = v[used]
    joined = Mesh(vertices, f, wall.frame)
    topo = topology(joined)
    assert not topo['closed_oriented_manifold']
    parent_part = np.asarray(face_part)[original_faces]
    parent_face = np.asarray(face_local)[original_faces]
    edges = defaultdict(list)
    for fi, face in enumerate(f):
        for a,b in zip(face, np.roll(face,-1)):
            edges[tuple(sorted((int(a),int(b))))].append(fi)
    open_edges = np.array([e for e,o in edges.items() if len(o) == 1], int)
    nonmanifold_edges = np.array([e for e,o in edges.items() if len(o) > 2], int)
    badlinks = np.array(topo['invalid_vertex_links'], int)
    # topology() reindexes every exact-coordinate vertex in sorted order;
    # vertices here are already the identical sorted unique source controls.
    assert np.array_equal(vertices, np.unique(vertices,axis=0))
    types = {}
    for key, ed in [('open-edge', open_edges), ('nonmanifold-edge', nonmanifold_edges)]:
        ids = sorted(set(fi for e in ed for fi in edges[tuple(e)]))
        types[key] = np.array(ids, int)
    types['invalid-vertex-link'] = np.flatnonzero(np.isin(f,badlinks).any(1))
    all_bad_faces = np.unique(np.concatenate(list(types.values())))
    patches = patch_components(joined, all_bad_faces)
    patch_of = {face: pi for pi, patch in enumerate(patches) for face in patch}
    rows = []
    for pi, patch in enumerate(patches):
        ids = np.asarray(patch, int)
        pts = vertices[np.unique(f[ids])]
        counts = Counter(parent_part[ids])
        rows.append(dict(patch=pi, normalized_source_faces=ids.tolist(), bounds_m=[pts.min(0).tolist(),pts.max(0).tolist()], axis_extent_mm=((pts.max(0)-pts.min(0))*1000).tolist(), bounds_diagonal_mm=float(np.linalg.norm(pts.max(0)-pts.min(0))*1000), defect_face_counts={k:int(np.isin(ids,x).sum()) for k,x in types.items()}, source_material_face_counts={selected[int(k)]['id']:int(c) for k,c in counts.items()}, branch_outer_surface_distances=[]))
    defect_mesh = Mesh(vertices,f[all_bad_faces],wall.frame)
    branch_checks = []
    complete = True
    threshold = 0.005
    for ci, name in enumerate(names[1:],1):
        assert time.monotonic()-began < 45, '45-second saved-stage budget'
        faces = bank['outer_faces'][bank['outer_component']==ci]
        branch = Mesh(bank['vertices'],faces,wall.frame)
        pairs, summary = pair_contacts(defect_mesh, branch, threshold, max_pairs=500000, max_tests=4000000)
        complete &= bool(summary['complete_within_pair_budget'])
        artifact = out / (name+'-defect-distance-pairs.npz')
        atomic_npz(artifact, **pairs)
        for row in rows:
            if len(pairs['pairs']):
                local_faces = pairs['pairs'][:,0]
                source_faces = all_bad_faces[local_faces]
                hit = np.isin(source_faces,row['normalized_source_faces'])
            else:
                hit = np.zeros(0,bool)
            d = dict(branch=name, search_radius_mm=threshold*1000, complete_within_pair_budget=summary['complete_within_pair_budget'])
            if np.any(hit):
                positions = np.flatnonzero(hit)
                nearest = int(positions[np.argmin(pairs['distance_m'][hit])])
                d.update(minimum_triangle_distance_mm=float(pairs['distance_m'][nearest]*1000), normalized_source_face=int(all_bad_faces[pairs['pairs'][nearest,0]]), branch_outer_face=int(pairs['pairs'][nearest,1]), closest_source_point_m=pairs['first_point'][nearest].tolist() if 'first_point' in pairs else None, closest_branch_point_m=pairs['second_point'][nearest].tolist() if 'second_point' in pairs else None)
            else:
                d.update(minimum_triangle_distance_lower_bound_mm=threshold*1000 if summary['complete_within_pair_budget'] else None)
            row['branch_outer_surface_distances'].append(d)
        branch_checks.append(dict(branch=name, actual_emitted_outer_geometry_sha256=branch.geometry_sha, actual_emitted_outer_face_count=len(faces), summary=summary, artifact=str(artifact.relative_to(ROOT)),artifact_sha256=digest_file(artifact), minimum_triangle_distance_mm=float(pairs['distance_m'].min()*1000) if len(pairs['pairs']) else None, minimum_triangle_distance_lower_bound_mm=threshold*1000 if not len(pairs['pairs']) and summary['complete_within_pair_budget'] else None))
        atomic_json(out/'progress.json',dict(defect_patches=rows,branch_checks=branch_checks,elapsed_seconds=time.monotonic()-began),overwrite=True)
    arrays = out/'literal-source-defects.npz'
    atomic_npz(arrays,vertices=vertices,faces=f,original_joined_face_indices=original_faces,source_part_indices=parent_part,source_part_local_face_indices=parent_face,open_edges=open_edges,nonmanifold_edges=nonmanifold_edges,invalid_vertex_link_vertices=badlinks,defect_faces=all_bad_faces)
    for member in members:
        # Keep identifiers and exact topology, while reducing repeated lists.
        for t in [member['original_topology']] + [c['topology'] for c in member['exact_surface_edge_components']]:
            t['invalid_vertex_link_count']=len(t.pop('invalid_vertex_links'))
    t=topo.copy();t['invalid_vertex_link_count']=len(t.pop('invalid_vertex_links'))
    result=dict(schema='literal-mandible-source-defect-localization/v32',source_catalogue_sha256=expected,registry_revision='31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2',candidate_geometry_sha256=wall.geometry_sha,source_ancestry_bank_sha256=digest_file(bank_path),frame=wall.frame,source_members=members,whole_joined_source_geometry_sha256=joined.geometry_sha,whole_source_topology=t,whole_edge_connected_component_count=len(components),whole_component_source_material_face_counts={selected[int(k)]['id']:int(c)for k,c in Counter(parent_part).items()},exact_surface_normalization=normalization,defect_patch_count=len(rows),defect_patches=rows,branch_checks=branch_checks,artifact=str(arrays.relative_to(ROOT)),artifact_sha256=digest_file(arrays),complete=complete,elapsed_seconds=time.monotonic()-began,occupied_classification='unknown-no-valid-literal-enclosure',original_source_grouping_does_not_establish_closed_bone=True,no_caps_or_new_classifier=True,no_native_shape_or_visibility_changes=True,no_installation=True)
    assert digest_file(ROOT/binding['current_catalogue'])==expected
    assert Mesh.load(work/'finite-wall.npz').geometry_sha==wall.geometry_sha
    atomic_json(out/'evidence.json',result)
    print(json.dumps(dict(output=str((out/'evidence.json').relative_to(ROOT)),complete=complete,elapsed_seconds=result['elapsed_seconds'],defect_patch_count=len(rows),patch_extents_mm=[r['axis_extent_mm']for r in rows],branch_minimums=[dict(branch=r['branch'],minimum_mm=r['minimum_triangle_distance_mm'],lower_bound_mm=r['minimum_triangle_distance_lower_bound_mm'])for r in branch_checks],occupied_classification=result['occupied_classification'])))


if __name__=='__main__':
    main()
