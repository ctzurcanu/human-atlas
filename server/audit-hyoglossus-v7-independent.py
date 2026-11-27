"""Read-only, actual-array audit of the selected nonplanar HG origin repair."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.audit import audit, topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature

WORK = ROOT / '.local-models/geometry-workbench/tongue-completion'
DOMAINS = {'left': [[60, 61, 103, 104], [70, 71, 113, 114]],
           'right': [[62, 63, 106, 107], [89, 90, 133, 134]]}


def canonical_triangle(triangle):
    points = [tuple(map(float, p)) for p in triangle]
    return min(tuple(points[k:] + points[:k]) for k in range(3))


def clear(summary):
    return (summary['complete_within_pair_budget'] and not summary['truncated']
            and summary['crossing_face_pairs'] == 0
            and summary['coplanar_overlap_face_pairs'] == 0)


def main(stage, output):
    start = time.monotonic()
    if output.exists():
        raise ValueError('Independent evidence output must be new')
    candidate_path = stage / 'candidate.json'
    report = json.loads(candidate_path.read_text())
    catalogue = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
    source_audit = WORK / 'source-v1/source-audit.json'
    catalogue_sha = digest_file(catalogue)
    if catalogue_sha != 'cbed9740aca9cd20ca12bf43a051eb7b36f60f90dc1919b64818dde8848f963f':
        raise ValueError('Parent source catalogue is stale')
    bones = {name: (WORK / f'source-v1/meshes/native-{idx:02}.npz') for name, idx in [
        ('mandible', 5), ('temporal-left', 6), ('temporal-right', 7),
        ('hyoid-body', 13), ('lesser-left', 14), ('lesser-right', 15),
        ('greater-left', 16), ('greater-right', 17)]}
    hyoid_v = []
    hyoid_f = []
    offset = 0
    for name, path in bones.items():
        if name in ('mandible', 'temporal-left', 'temporal-right'):
            continue
        bone = Mesh.load(path)
        hyoid_v.append(bone.vertices)
        hyoid_f.append(bone.faces + offset)
        offset += len(bone.vertices)
    hyoid = Mesh(np.vstack(hyoid_v), np.vstack(hyoid_f), 'model:local-ta98:assembled-metres')
    mucosa_path = WORK / 'exposed-mucosa-v1/packed/TA98TONGUE_male_mucosa.npz'
    contexts = {**bones, 'mucosa': mucosa_path}
    context_bindings = {name: dict(path=str(p), sha256=digest_file(p),
                                  geometry_sha256=Mesh.load(p).geometry_sha)
                        for name, p in contexts.items()}
    output.mkdir(parents=True)
    rows = []
    for row in report['rows']:
        if time.monotonic() - start > 45:
            raise TimeoutError('Independent HG audit reached its finite stage budget')
        side = row['side']
        path = Path(row['candidate'])
        source_path = Path(row['source'])
        footprint_path = Path(row['nativeFootprint'])
        if digest_file(path) != row['candidateSha256'] or digest_file(source_path) != row['sourceSha256']:
            raise ValueError('Actual HG candidate/source hash mismatch')
        if digest_file(footprint_path) != row['nativeFootprintSha256']:
            raise ValueError('Actual cortical origin footprint hash mismatch')
        mesh = Mesh.load(path)
        source = Mesh.load(source_path)
        source_words = source.attributes['normals'][source.faces]
        words = mesh.attributes['normals'][mesh.faces]
        source_map = {canonical_triangle(t): i for i, t in enumerate(source.triangles)}
        target_map = {canonical_triangle(t): i for i, t in enumerate(mesh.triangles)}
        protected = np.flatnonzero(np.all(source.triangles[:, :, 1] >= 1.621, axis=1))
        preservation = []
        for face in protected:
            target = target_map.get(canonical_triangle(source.triangles[face]))
            equal = target is not None
            if equal:
                for corner, point in enumerate(source.triangles[face]):
                    target_corner = np.flatnonzero(np.all(mesh.triangles[target] == point, axis=1))
                    equal = equal and len(target_corner) == 1 and np.array_equal(
                        source_words[face, corner], words[target, target_corner[0]])
            preservation.append(dict(source_face=int(face), candidate_face=target,
                                     coordinates_and_raw_words_literal=bool(equal)))
        actual = audit(mesh, require_closed=True, expected_components=1,
                       expected_genus=[0], max_pairs=100000, max_tests=1000000)
        contacts = {}
        for name, context_path in contexts.items():
            _, summary = pair_contacts(mesh, Mesh.load(context_path), 0.,
                                       max_pairs=100000, max_tests=1000000)
            contacts[name] = summary
        with np.load(stage / f'HG-{side}-stored.npz', allow_pickle=False) as stored:
            decoded = stored['positions'].astype(float) + stored['geometryOrigin']
            storage_exact = (np.array_equal(decoded, mesh.vertices)
                             and np.array_equal(stored['faces'], mesh.faces)
                             and np.array_equal(stored['attribute_normals'], mesh.attributes['normals']))
        tri = mesh.triangles
        normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
        normal /= np.linalg.norm(normal, axis=1)[:, None]
        decoded_normals = words.astype(float)
        lengths = np.linalg.norm(decoded_normals, axis=2)
        decoded_normals /= lengths[:, :, None]
        dot = np.einsum('fki,fi->fk', decoded_normals, normal)
        preserved_faces = set(p['candidate_face'] for p in preservation if p['candidate_face'] is not None)
        generated_mask = ~np.isin(np.arange(len(tri)), list(preserved_faces))
        normal_result = dict(storage_dtype=str(mesh.attributes['normals'].dtype),
                             nonzero_corners=bool(np.all(lengths > 0)),
                             minimum_corner_face_dot=float(dot.min()),
                             nonpositive_corner_count=int(np.sum(dot <= 0)),
                             minimum_nonprotected_corner_face_dot=float(dot[generated_mask].min()),
                             protected_faces_with_nonpositive_corner=int(np.sum(np.any(dot[list(preserved_faces)] <= 0, axis=1))))
        unique = np.unique(mesh.vertices, axis=0)
        distance2, bone_faces, closest = igl.point_mesh_squared_distance(unique, hyoid.vertices, hyoid.faces)
        near = distance2 < 1e-12
        near_distances = np.sqrt(distance2[near])
        near_winding = np.abs(igl.winding_number(hyoid.vertices, hyoid.faces, unique[near]))
        adjacency = dict(actual_near_controls=int(np.sum(near)),
                         actual_gap_range_m=[float(near_distances.min()), float(near_distances.max())],
                         actual_hyoid_face_witnesses=bone_faces[near].tolist(),
                         maximum_absolute_winding=float(near_winding.max()),
                         within_declared_one_micrometre_tolerance=bool(near_distances.max() <= 1e-6),
                         histological_enthesis_established=False,
                         literal_shared_cortical_chart=False)
        old = Mesh.load(WORK / f'hyoglossus-origin-v4-replaced-lower/{side}-original-broad-origin-sheet.npz')
        sheet = Mesh.load(stage / f'{side}-original-broad-origin-sheet.npz')
        moving = np.concatenate(DOMAINS[side])
        stationary = ~np.isin(np.arange(len(old.faces)), moving)
        patch_proof = []
        for ids in DOMAINS[side]:
            old_ids = np.unique(old.faces[ids])
            patch_proof.append(dict(face_rows=ids,
                                    boundary_literal=boundary_signature(old.faces[ids]) == boundary_signature(sheet.faces[ids]),
                                    exact_same_controls=set(old_ids) == set(np.unique(sheet.faces[ids])),
                                    new_triangles_use_only_original_boundary_controls=bool(np.all(np.isin(sheet.faces[ids], old_ids))),
                                    unchanged_control_bounds_m=[old.vertices[old_ids].min(0).tolist(), old.vertices[old_ids].max(0).tolist()]))
        topology_proof = dict(all_coordinates_literal=np.array_equal(old.vertices, sheet.vertices),
                              all_unaffected_face_rows_literal=np.array_equal(old.faces[stationary], sheet.faces[stationary]),
                              changed_face_rows=int(np.sum(np.any(old.faces != sheet.faces, axis=1))),
                              patches=patch_proof,
                              old_self_intersecting_integral_not_occupied_volume=True,
                              old_algebraic_integral_m3=topology(old)['oriented_triangle_integral_m3'],
                              repaired_occupied_volume_m3=topology(sheet)['oriented_triangle_integral_m3'])
        passed = bool(actual['passed'] and all(clear(s) for s in contacts.values())
                      and all(p['coordinates_and_raw_words_literal'] for p in preservation)
                      and storage_exact and normal_result['nonzero_corners']
                      and normal_result['nonpositive_corner_count'] == 0
                      and adjacency['actual_near_controls'] == row['adjacency']['contactControls']
                      and adjacency['within_declared_one_micrometre_tolerance']
                      and adjacency['maximum_absolute_winding'] < .5
                      and topology_proof['all_coordinates_literal']
                      and topology_proof['all_unaffected_face_rows_literal']
                      and all(p['boundary_literal'] and p['exact_same_controls']
                              and p['new_triangles_use_only_original_boundary_controls'] for p in patch_proof))
        result = dict(side=side, passed=passed, mesh=str(path), sha256=digest_file(path),
                      geometry_sha256=mesh.geometry_sha, actual_audit=actual,
                      actual_context_contacts=contacts, protected_source_faces=len(protected),
                      protected_geometry_and_raw_words_literal=sum(p['coordinates_and_raw_words_literal'] for p in preservation),
                      protected_face_map=preservation, actual_storage_roundtrip_exact=storage_exact,
                      normal_field=normal_result, cortical_adjacency=adjacency,
                      rim_topology_only_proof=topology_proof,
                      source_footprint=str(footprint_path), source_footprint_sha256=digest_file(footprint_path),
                      footprint_course_scope='Source-bound hyoid body/greater-cornu footprint, positive affine narrowing to the native HG section; exact histological fibre trajectories are not established',
                      actual_normal_chart_vertices=len(mesh.vertices), actual_triangles=len(mesh.faces),
                      estimated_runtime_bytes=18*len(mesh.vertices)+12*len(mesh.faces))
        rows.append(result)
        atomic_json(output / f'{side}-independent-audit.json', result)
    sources_unchanged = (digest_file(catalogue) == catalogue_sha and
                         all(digest_file(p['path']) == p['sha256'] for p in context_bindings.values()))
    result = dict(schema='hyoglossus-v7-independent-actual-array-audit/v1',
                  candidate=str(candidate_path), candidate_sha256=digest_file(candidate_path),
                  source_catalogue=dict(path=str(catalogue), sha256=catalogue_sha),
                  source_audit=dict(path=str(source_audit), sha256=digest_file(source_audit)),
                  fixed_context=context_bindings, source_and_fixed_catalogues_unchanged=sources_unchanged,
                  rows=rows, passed=sources_unchanged and all(r['passed'] for r in rows),
                  anatomical_acceptance=False, installed=False,
                  elapsed_seconds=time.monotonic()-start)
    atomic_json(output / 'summary.json', result)
    print(json.dumps(dict(passed=result['passed'], elapsed_seconds=result['elapsed_seconds'],
                          sides=[dict(side=r['side'], passed=r['passed'], protected=r['protected_geometry_and_raw_words_literal'],
                                      normals=r['normal_field'], adjacency=r['cortical_adjacency']['actual_gap_range_m'],
                                      self_defects=r['actual_audit']['self_intersections']['defect_pairs'],
                                      context_defects=sum(s['crossing_face_pairs']+s['coplanar_overlap_face_pairs'] for s in r['actual_context_contacts'].values()))
                                 for r in rows])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    main(args.stage.resolve(), args.output.resolve())
