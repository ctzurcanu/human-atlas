"""Finite native-facial side openings using the lingual construction pipeline.

No overlapping cylinders or retained artificial proximal caps are accepted.
Ownership packing, actual Float32 replay, context and render review are separate
requirements; this bounded construction alone never permits installation.
"""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import manifold3d as md
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames, path_tangents
from geometry_tools.algorithms.tube_continuity import boundary_loops, loft_faces, annular_rim, finite_wall_cells, continuity_audit
from geometry_tools.algorithms.finite_tissue_boolean import make_solid, provenance, CommonPositionBank
from geometry_tools.audit import audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--route-version', default='v2')
    parser.add_argument('--output-version', default='v1')
    parser.add_argument('--max-pairs',type=int,default=250000)
    parser.add_argument('--max-tests',type=int,default=2000000)
    args = parser.parse_args()
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/circulatory-completion'
    route_file = base / ('male-left-facial-cervical-routes-' + args.route_version) / 'routes.json'
    plan = json.loads(route_file.read_bytes())
    assert digest_file(ROOT / plan['source_catalogue']) == plan['source_catalogue_sha256']
    directory = base / ('male-left-facial-cervical-finite-network-' + args.output_version)
    assert not directory.exists()
    directory.mkdir()
    parent = Mesh.load(ROOT / plan['fixed_parent']['mesh'])
    rings = np.asarray(plan['fixed_parent']['source_station_components'][0]['rings'])
    centers = parent.vertices[rings].mean(1)
    origin = centers.mean(0)
    native_inner = parent.vertices.copy()
    for ring, center in zip(rings, centers):
        delta = parent.vertices[ring] - center
        native_inner[ring] -= plan['illustrative_parent_wall_m'] * delta / np.linalg.norm(delta, axis=1)[:, None]
    fields = [np.array([parent.vertices, native_inner])]
    carriers = [parent.faces]
    all_rings = [rings]
    names = ['facial-left'] + [r['name'] for r in plan['routes']]
    terms = ['A12.2.05.020'] + [r['term'] for r in plan['routes']]
    guides = []
    for route in plan['routes']:
        guide = ROOT / route['guide']
        assert digest_file(guide) == route['guide_sha256']
        archive = np.load(guide)
        path = archive['points']
        outer_radii, inner_radii = archive['outer_radii'], archive['inner_radii']
        source_ring = rings[route['parent_station']]
        delta = parent.vertices[source_ring] - centers[route['parent_station']]
        _, _, vt = np.linalg.svd(delta, full_matrices=False)
        normal = vt[-1]
        if normal @ (centers[route['parent_station'] + 1] - centers[route['parent_station'] - 1]) < 0:
            normal = -normal
        reference = delta[0] - normal * (delta[0] @ normal)
        reference /= np.linalg.norm(reference)
        basis = np.c_[reference, np.cross(normal, reference), normal]
        profile = delta @ basis
        profile[:, 2] = 0
        profile /= np.linalg.norm(profile, axis=1).mean()
        if 'source_profile_template' in route:
            profile = np.asarray(route['source_profile_template'], float)
            assert profile.shape == delta.shape and np.all(profile[:, 2] == 0)
        tangent = archive['analytic_tangents'] if 'analytic_tangents' in archive else path_tangents(path)
        assert tangent.shape == path.shape and np.isfinite(tangent).all()
        frames = rotation_minimizing_frames(path, tangent, initial_reference=reference)
        section = np.einsum('kij,nj->kni', frames, profile)
        outer = path[:, None] + outer_radii[:, None, None] * section
        inner = path[:, None] + inner_radii[:, None, None] * section
        face = loft_faces(len(path), len(source_ring))
        field = np.array([outer.reshape(-1, 3), inner.reshape(-1, 3)])
        check = finite_wall_cells(field, face)
        atomic_json(directory / (route['name'] + '-paired-field.json'), check)
        atomic_npz(directory / (route['name'] + '-fields.npz'), fields=field, faces=face, guide_points=path, frame=np.array(parent.frame))
        assert check['finite_wall_cells_pass'], (route['name'], check['invalid_counts_by_interval'])
        fields.append(field)
        carriers.append(face)
        all_rings.append(np.arange(len(path) * len(source_ring)).reshape(-1, len(source_ring)))
        guides.append(path)
    ports = [(0, 0), (0, len(rings) - 1)] + [(i, len(all_rings[i]) - 1) for i in range(1, len(fields))]
    atomic_json(directory / 'construction-plan.json', dict(source_catalogue_sha256=plan['source_catalogue_sha256'], route_contract_sha256=digest_file(route_file),
                tissue_names=names, terms=terms, expected_external_ports=ports, native_parent_outer_coordinates_fixed=True,
                private_geometry_copied=False, installed=False, anatomical_acceptance=False))
    input_owners = {}
    sheets = []
    ancestries = []
    blood = None
    for boundary in [0, 1]:
        solids = []
        for component, (field, face, component_rings) in enumerate(zip(fields, carriers, all_rings)):
            vertices = field[boundary].copy()
            triangles = face.copy()
            caps = {}
            for loop in boundary_loops(triangles):
                polygon = vertices[loop]
                order = np.lexsort((polygon[:, 2], polygon[:, 1], polygon[:, 0]))
                center_index = len(vertices)
                vertices = np.r_[vertices, polygon[order].mean(0)[None]]
                cap = np.array([[center_index, int(loop[(j + 1) % len(loop)]), int(loop[j])] for j in range(len(loop))])
                external = any(component == ci and set(loop) == set(component_rings[station]) for ci, station in ports)
                caps.update({int(k): external for k in np.arange(len(triangles), len(triangles) + len(cap))})
                triangles = np.r_[triangles, cap]
            solid, identifier = make_solid((vertices - origin) * 1000, triangles)
            input_owners[identifier] = dict(component=component, caps=caps)
            solids.append(solid)
        joined = md.Manifold.batch_boolean(solids, md.OpType.Add)
        assert joined.status() == md.Error.NoError
        raw = joined.to_mesh64()
        vertices = np.asarray(raw.vert_properties)[:, :3] / 1000 + origin
        triangles = np.asarray(raw.tri_verts, int)
        owner, source_face = provenance(raw)
        external_cap = np.array([input_owners[int(o)]['caps'].get(int(f)) is True for o, f in zip(owner, source_face)])
        internal_cap = np.array([input_owners[int(o)]['caps'].get(int(f)) is False for o, f in zip(owner, source_face)])
        cap_report = dict(boundary=boundary, retained_internal_caps=int(internal_cap.sum()), retained_external_caps=int(external_cap.sum()),
                          internal_by_component={names[c]: sum(input_owners[int(o)]['component'] == c for o in owner[internal_cap]) for c in range(len(names))})
        atomic_json(directory / (('inner' if boundary else 'outer') + '-classifier-caps.json'), cap_report)
        assert not internal_cap.any(), ('An internal classifier cap survived', cap_report)
        bank = CommonPositionBank(np.concatenate([field[boundary] for field in fields]), tolerance=2e-12)
        mapping = np.array([bank.insert(point) for point in vertices])
        sheet = Mesh(bank.positions, mapping[triangles[~external_cap]], parent.frame)
        name = 'inner' if boundary else 'outer'
        sheet.save(directory / (name + '-sheet.npz'))
        check = continuity_audit(sheet.vertices, sheet.faces, expected_mouths=len(ports), expected_components=1)
        atomic_json(directory / (name + '-continuity.json'), check)
        assert check['topological_continuity_pass'], check
        sheets.append(sheet)
        ancestries.append(dict(component=np.array([input_owners[int(o)]['component'] for o in owner[~external_cap]]), source_face=source_face[~external_cap]))
        if boundary:
            blood = Mesh(bank.positions, mapping[triangles], parent.frame)
            blood.save(directory / 'closed-blood-classifier.npz')
    outer, inner = sheets
    offset = len(outer.vertices)
    vertices = np.r_[outer.vertices, inner.vertices]
    faces = [outer.faces, inner.faces[:, ::-1] + offset]
    outer_loops, inner_loops = boundary_loops(outer.faces), boundary_loops(inner.faces)
    port_records = []
    for component, station in ports:
        expected = fields[component][:, all_rings[component][station]]
        first = min(outer_loops, key=lambda loop: np.linalg.norm(outer.vertices[loop].mean(0) - expected[0].mean(0)))
        second = min(inner_loops, key=lambda loop: np.linalg.norm(inner.vertices[loop].mean(0) - expected[1].mean(0)))
        order = np.argmin(np.linalg.norm(expected[0][None] - outer.vertices[first][:, None], axis=2), axis=1)
        target = expected[1][order]
        permutation = np.argmin(np.linalg.norm(inner.vertices[second][None] - target[:, None], axis=2), axis=1)
        assert len(set(permutation)) == len(second) == 12
        assert np.linalg.norm(inner.vertices[second[permutation]] - target, axis=1).max() < 1e-11
        annulus = annular_rim(first, second[permutation] + offset)
        faces.append(annulus)
        port_records.append(dict(component=component, source_station=station, outer_indices=first.tolist(), inner_indices=(second[permutation] + offset).tolist(), recipient=plan['routes'][component - 1]['recipient'] if component else None))
    wall = Mesh(vertices, np.concatenate(faces), parent.frame)
    wall.save(directory / 'finite-wall.npz')
    wall_check = audit(wall, require_closed=True, expected_components=1, expected_genus=[len(ports) - 1],max_pairs=args.max_pairs,max_tests=args.max_tests)
    atomic_json(directory / 'wall-audit.json', wall_check)
    assert wall_check['passed'], wall_check['failures']
    blood_check = audit(blood, require_closed=True, expected_components=1, expected_genus=[0],max_pairs=args.max_pairs,max_tests=args.max_tests)
    atomic_json(directory / 'blood-audit.json', blood_check)
    assert blood_check['passed'], blood_check['failures']
    atomic_npz(directory / 'surfaces-and-ancestry.npz', vertices=vertices, outer_faces=outer.faces, inner_faces=inner.faces[:, ::-1] + offset,
               outer_component=ancestries[0]['component'], outer_source_face=ancestries[0]['source_face'], inner_component=ancestries[1]['component'], inner_source_face=ancestries[1]['source_face'], frame=np.array(parent.frame), **{f'annulus_{i}_faces': face for i, face in enumerate(faces[2:])})
    atomic_json(directory / 'ports.json', dict(ports=port_records, component_names=names, component_terms=terms))
    result = dict(schema='finite-left-facial-cervical-network/v1', source_catalogue_sha256=plan['source_catalogue_sha256'], geometry_sha256=wall.geometry_sha,
                  new_terms=sorted(set(terms[1:])), branch_tubes=len(guides), open_external_blood_ports=len(ports), wall_closed_genus=len(ports)-1,
                  wall_self_crossings=wall_check['self_intersections']['defect_pairs'], blood_self_crossings=blood_check['self_intersections']['defect_pairs'], internal_classifier_caps=0,
                  native_parent_outer_coordinates_retained=True, shared_owner_partition_pending=True, full_context_and_emitted_bank_verification_pending=True,
                  installed=False, installable=False, anatomical_acceptance=False, elapsed_seconds=time.monotonic()-began)
    atomic_json(directory / 'evidence.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
