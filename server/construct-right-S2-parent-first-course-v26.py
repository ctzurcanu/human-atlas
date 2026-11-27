"""One source-owned arterial forest construction, independent of airway ends.

Actual source parent stations own child prefixes. C2 suffixes use local
parenchymal/airway surface halfspaces; these are proposals, not clearance proof.
"""
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.spatial.transform import Rotation

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH
from geometry_tools.audit import winding_contains, topology
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.tube_sections import arclength, path_tangents, bend_diagnostics
from geometry_tools.algorithms.tube_c2_course import C2Course

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OUT = BASE / 'right-S2-parent-first-course-v26'
FRAME = 'model:male:assembled-metres'
S2 = 'TA98PULMONARY:male:right-superior:posterior'


def proper(a, b):
    a = a / np.linalg.norm(a, axis=-1)[:, None]
    b = b / np.linalg.norm(b, axis=-1)[:, None]
    dot = np.sum(a * b, axis=1)
    if np.any(dot < -1 + 1e-8):
        raise ValueError('Own-vessel tangent reversal cannot be hidden by profile transport')
    v = np.cross(a, b)
    r = np.tile(np.eye(3), (len(a), 1, 1))
    for i, (xyz, cosine) in enumerate(zip(v, dot)):
        x, y, z = xyz
        k = np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])
        r[i] += k + k @ k / (1 + cosine)
    return r


def main():
    began = time.monotonic(); OUT.mkdir(exist_ok=True)
    family = json.loads((BASE / 'family-contract-v2.json').read_text())
    assert digest_file(Path(family['source_catalogue'])) == family['source_catalogue_sha256']
    atlas = json.loads(Path(family['source_catalogue']).read_text()); parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    old = json.loads((BASE / 'right-S2-ordered-chart-contract-v21/contract.json').read_text())
    arrays = np.load(BASE / 'right-S2-ordered-chart-contract-v21/literal-ordered-guides-and-roots.npz', allow_pickle=False)
    selected = old['selected_source_curves']; bykey = {r['key']: r for r in selected}
    sources = {r['key']: arrays[r['source_sections_array']] for r in selected}
    centres = {k: q.mean(1) for k, q in sources.items()}
    roots = {r['key']: (len(centres[r['key']]) - 1 if r['larger_radius_endpoint'] else 0) for r in selected}
    owners = {}; pending = set(bykey); order = []
    for key in ['p7_c1', 'p31_c1']:
        owners[key] = None; order.append(key); pending.remove(key)
    while pending:
        progress = False
        for key in sorted(pending):
            p = centres[key][roots[key]]
            matches = []
            for parent in order:
                d = np.linalg.norm(centres[parent] - p, axis=1); station = int(d.argmin())
                if d[station] < 1e-6:
                    matches.append((float(d[station]), parent, station))
            if matches:
                distance, parent, station = min(matches)
                owners[key] = dict(parent=parent, parent_station=station, child_station=roots[key],
                                   literal_centre_gap_m=distance, finite_surface_interface_unproven=True)
                order.append(key); pending.remove(key); progress = True; break
        if not progress:
            raise ValueError('Source parent station ownership unresolved: ' + repr(sorted(pending)))
    territory = reader.read_part(atlas, parts[S2], frame=FRAME)
    assert topology(territory)['closed_oriented_manifold']
    territory_tree = BVH(territory)
    tv = territory.triangles
    face_normals = np.cross(tv[:, 1] - tv[:, 0], tv[:, 2] - tv[:, 0])
    face_normals /= np.linalg.norm(face_normals, axis=1)[:, None]
    sign = np.sign(np.einsum('ij,ij->i', tv[:, 0], np.cross(tv[:, 1], tv[:, 2])).sum())
    face_normals *= sign
    airway_ids = [r['id'] for r in json.loads((BASE / 'full-context-participation-v11/evidence.json').read_text())['records'] if r['retained'] and r['role'] == 'current pulmonary airway surface']
    # Role wording is not authority: use the explicitly registered 27-piece set.
    if not airway_ids:
        airway_ids = [p['id'] for p in atlas['parts'] if p['id'].startswith('LOCAL:male:pulmonary:bronchus:')]
    assert len(airway_ids) == 27
    vv = []; ff = []; offset = 0
    for pid in airway_ids:
        mesh = reader.read_part(atlas, parts[pid], frame=FRAME)
        vv.extend(mesh.vertices); ff.extend(mesh.faces + offset); offset += len(mesh.vertices)
    airway = Mesh(np.asarray(vv), np.asarray(ff), FRAME); airway_tree = BVH(airway)
    contract = dict(schema='right-S2-parent-first-source-owned-course-contract/v26',
                    source_catalogue_sha256=family['source_catalogue_sha256'], source_parent_station_owners=owners,
                    construction_order=order, source_identity_review_sha256=old['source_semantic_root_review_sha256'],
                    fixed_territory=dict(id=S2, geometry_sha256=territory.geometry_sha),
                    fixed_airway_ids=airway_ids, fixed_airway_geometry_sha256=airway.geometry_sha,
                    fixed_native_origins=['p7_c1:last', 'p31_c1:last'], protected_prefix_m=.004,
                    distal_parenchymal_constraint_after_root_m=.020, local_surface_clearance_margin_m=.0006,
                    suffix_minimum_bend_radius_factor_of_largest_profile=1/.75,
                    source_parameter_speed_lower_bound=.60, suffix_trust_radius_mm=12., solver_seconds_per_suffix=3.,
                    airway_truncation_endpoints_not_targets=True, source_calibre_preserved=True,
                    root_position_first_and_second_derivative_jets_fixed=True,
                    child_prefixes_carry_parent_station_proper_transform=True,
                    local_halfspaces_do_not_establish_global_containment_or_clearance=True,
                    current_all_54_family_sources_and_1524_context_owners_required=True,
                    finite_junctions_and_patency_unproven=True, catalogue_changed=False, anatomical_acceptance=False)
    (OUT / 'contract.json').write_text(json.dumps(contract, indent=2) + '\n')
    transformed = {}; diagnostics = []; saved = {}
    for key in order:
        source = sources[key]; cc = centres[key]; own = owners[key]
        root = roots[key]; reverse = root != 0
        if own:
            parent = transformed[own['parent']]; j = own['parent_station']
            rr = parent['rotations'][j]
            parent_old = centres[own['parent']][j]; parent_new = parent['centres'][j]
            rigid = parent_new + np.einsum('ij,svj->svi', rr, source - parent_old)
        else:
            rr = np.eye(3); rigid = source.copy()
        initial = rigid.mean(1); oriented = rigid[::-1] if reverse else rigid
        initial_oriented = oriented.mean(1); s = arclength(initial_oriented)
        prefix_last = min(int(np.searchsorted(s, .004)), len(s) - 4)
        u = (s[prefix_last:] - s[prefix_last]) * 1000
        initial_free = initial_oriented[prefix_last:]
        radius = np.linalg.norm(oriented - initial_oriented[:, None], axis=2).max(1)
        inside = winding_contains(territory, initial_free, distance_backend='bvh')['inside']
        targets = initial_free.copy(); planes = []; target_witnesses = []
        for j, (point, inside_point) in enumerate(zip(initial_free, inside)):
            full_station = prefix_last + j; original_distance = s[full_station]
            hit = territory_tree.closest_point(point); normal = face_normals[hit['face_id']]
            margin = radius[full_station] + .0006
            if original_distance > .020:
                if not inside_point or hit['distance_m'] < margin:
                    targets[j] = hit['point'] - normal * margin
                planes.append((float(u[j]), -normal, float((-normal @ hit['point'] + margin) * 1000), dict(owner=S2, face=int(hit['face_id']))))
            target_witnesses.append(dict(station=full_station, original_inside=bool(inside_point),
                                         source_to_boundary_m=float(hit['distance_m']), boundary_face=int(hit['face_id'])))
        targets[0] = initial_free[0]
        # Surface clearance proposals about the actual current airway, never
        # physiological endpoints or a copied donor bronchial tangent field.
        for j, point in enumerate(targets):
            full_station = prefix_last + j
            if j == 0: continue
            hit = airway_tree.closest_point(point)
            normal = point - hit['point']; distance = np.linalg.norm(normal)
            margin = radius[full_station] + .0006
            if distance < 1e-10:
                raise ValueError('No geometrically defined airway-clearance direction at ' + key)
            if distance < margin + .003:
                normal /= distance
                planes.append((float(u[j]), normal, float((normal @ hit['point'] + margin) * 1000), dict(owner='fixed airway surface', face=int(hit['face_id']))))
        count = min(20, max(10, len(u) // 3))
        control_u = np.linspace(0., u[-1], count)
        reference = CubicSpline(u, targets * 1000, axis=0)(control_u)
        source_spline = CubicSpline(s * 1000, initial_oriented * 1000, axis=0)
        root_parameter = s[prefix_last] * 1000
        course = C2Course(reference, float(u[-1]), count, float(radius.max() * 1000 / .75),
                          speed_lower=.60, trust_radius_mm=12., boundary_condition='not-a-knot')
        equality = [(0., 1, source_spline(root_parameter, 1), key + ':prefix tangent'),
                    (0., 2, source_spline(root_parameter, 2), key + ':prefix curvature')]
        # Retain the literal linearised problem, including failed proposals;
        # infeasibility is local and must be inspectable without another fit.
        problem_path = OUT / (key + '-literal-conic-problem.npz')
        if not problem_path.exists():
            np.savez_compressed(problem_path, native_source_sections=source,
                                owned_prefix_sections=oriented[:prefix_last + 1],
                                source_station_parameters_mm=u, proposed_reference_controls_mm=reference,
                                proposed_station_targets_m=targets,
                                plane_parameters_mm=np.asarray([p[0] for p in planes]),
                                plane_normals=np.asarray([p[1] for p in planes]),
                                plane_minimum_dot_mm=np.asarray([p[2] for p in planes]),
                                plane_owner_json=np.asarray(json.dumps([p[3] for p in planes])),
                                root_tangent_jet=equality[0][2], root_curvature_jet=equality[1][2])
        controls, report = course.solve(planes, time_limit=3., equalities=equality)
        diagnostic = dict(key=key, source_parent=own, prefix_last_oriented_station=prefix_last,
                          local_surface_planes=len(planes), target_boundary_witnesses=target_witnesses, solver=report,
                          literal_problem=str(problem_path.resolve()), literal_problem_sha256=digest_file(problem_path))
        diagnostics.append(diagnostic)
        if controls is None or not report['certificate']['continuous_curvature_bound_pass']:
            failed = dict(contract_sha256=digest_file(OUT / 'contract.json'), candidate_emitted=False,
                          failed_curve=key, diagnostics=diagnostics, installed=False, anatomical_acceptance=False,
                          elapsed_seconds=time.monotonic() - began)
            (OUT / 'construction-evidence.json').write_text(json.dumps(failed, indent=2) + '\n')
            print(json.dumps(dict(output=str(OUT), failed_curve=key, status=report['status'], elapsed_seconds=failed['elapsed_seconds'])))
            return
        new_centres = initial_oriented.copy(); new_centres[prefix_last:] = course.evaluate(controls, u) / 1000
        tangents = path_tangents(initial_oriented)
        new_tangents = tangents.copy(); new_tangents[prefix_last:] = course.evaluate(controls, u, 1)
        rotations = proper(tangents, new_tangents)
        sections = new_centres[:, None] + np.einsum('sij,svj->svi', rotations, oriented - initial_oriented[:, None])
        sections[:prefix_last + 1] = oriented[:prefix_last + 1]
        rotations[:prefix_last + 1] = np.eye(3)
        if reverse:
            sections = sections[::-1]; new_centres = new_centres[::-1]; rotations = rotations[::-1]
        rotations = rotations @ rr
        transformed[key] = dict(sections=sections, centres=sections.mean(1), rotations=rotations)
        saved[key + '_source_sections'] = source; saved[key + '_candidate_sections'] = sections
        saved[key + '_proper_rotations'] = rotations
        diagnostic.update(maximum_source_displacement_m=float(np.linalg.norm(source - sections, axis=2).max()),
                          source_length_m=float(arclength(cc)[-1]), candidate_length_m=float(arclength(sections.mean(1))[-1]),
                          sampled_bend_max_support=float(bend_diagnostics(sections.mean(1), sections[None])['curvature_support'].max()))
    # Actual literal fragments first; no whole-context work starts on a failed
    # cheap self gate. Every other whole-family source remains literal below.
    maps = defaultdict(dict); normal_maps = defaultdict(dict); fragment_paths = []; self_rows = []
    for row in selected:
        key = row['key']; source = sources[key]; target = transformed[key]['sections']; rotations = transformed[key]['rotations']
        whole = reader.read_part(atlas, parts[row['part_id']], frame=FRAME)
        ff = arrays[row['source_faces_array']]; used, inv = np.unique(whole.faces[ff], return_inverse=True)
        lookup = {tuple(p): (i, j) for i, ring in enumerate(source) for j, p in enumerate(ring)}
        sr = np.array([lookup[tuple(p)] for p in whole.vertices[used]], int)
        vv = target[sr[:, 0], sr[:, 1]]
        fragment = Mesh(vv, inv.reshape(-1, 3), FRAME, source=dict(source_part_id=row['part_id'], source_component=row['component'], source_face_ids=ff.tolist()))
        path = OUT / (key + '-candidate-fragment.npz'); fragment.save(path)
        summary, _ = self_contact_summary(fragment, max_pairs=100000, max_tests=500000, stop_on_defect=True)
        fragment_paths.append(dict(part_id=row['part_id'], component=row['component'], mesh=str(path.resolve()), geometry_sha256=fragment.geometry_sha))
        self_rows.append(dict(key=key, summary=summary))
        for si, (oldring, newring) in enumerate(zip(source, target)):
            for a, b in zip(oldring, newring): maps[row['part_id']][tuple(a)] = b; normal_maps[row['part_id']][tuple(a)] = rotations[si]
    whole_paths = []
    for pi, pid in enumerate(json.loads((BASE / 'takeoff-and-root-controls-v10/evidence.json').read_text())['whole_runtime_qualified_blood_source_ids']):
        mesh = reader.read_part(atlas, parts[pid], frame=FRAME); vv = mesh.vertices.copy(); attrs = {k:v.copy() for k,v in mesh.attributes.items()}
        if pid in maps:
            for vi, point in enumerate(mesh.vertices):
                key = tuple(point)
                if key not in maps[pid]: continue
                vv[vi] = maps[pid][key]
                if 'normals' in attrs: attrs['normals'][vi] = np.rint(np.clip(normal_maps[pid][key] @ np.asarray(mesh.attributes['normals'][vi], float), -32767, 32767)).astype(np.int16)
        whole = Mesh(vv, mesh.faces, FRAME, attrs, dict(part_id=pid, source_geometry_sha256=mesh.geometry_sha, finite_wall_acceptance=False))
        path = OUT / ('whole-family-part-%02d.npz' % pi); whole.save(path)
        whole_paths.append(dict(id=pid, mesh=str(path.resolve()), geometry_sha256=whole.geometry_sha, source_geometry_sha256=mesh.geometry_sha))
    np.savez_compressed(OUT / 'literal-parent-first-ring-transport.npz', **saved)
    evidence = dict(schema='right-S2-parent-first-proper-section-candidate/v26', source_catalogue_sha256=family['source_catalogue_sha256'],
                    contract_sha256=digest_file(OUT / 'contract.json'), candidate_meshes=whole_paths,
                    selected_fragment_meshes=fragment_paths, diagnostics=diagnostics, self_gate=self_rows,
                    self_gate_pass=all(r['summary']['complete_within_pair_budget'] and not r['summary']['defect_pairs'] for r in self_rows),
                    source_calibre_scaled=False, airway_truncation_endpoint_attraction_removed=True,
                    whole_context_clearance_pending=True, finite_wall_lumen_and_patency_pending=True,
                    catalogue_changed=False, installed=False, anatomical_acceptance=False, elapsed_seconds=time.monotonic() - began)
    (OUT / 'candidate-evidence.json').write_text(json.dumps(evidence, indent=2) + '\n')
    print(json.dumps(dict(output=str(OUT), evidence_sha256=digest_file(OUT / 'candidate-evidence.json'), self_gate_pass=evidence['self_gate_pass'],
                          self_defects=[dict(key=r['key'], positive_certificate=r['summary']['defect_pairs']) for r in self_rows], elapsed_seconds=evidence['elapsed_seconds'])))


if __name__ == '__main__': main()
