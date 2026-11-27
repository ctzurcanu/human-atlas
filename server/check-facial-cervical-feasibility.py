"""Independent whole-guide, tissue-mouth and fixed-context feasibility gates."""
from pathlib import Path
import argparse
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
import igl
import manifold3d as md
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.finite_tissue_boolean import make_solid


def occupancy_classifier(mesh, allow_conduit_caps=False):
    """Actual occupied solid; named conduit caps are temporary classifiers only."""
    gate = topology(mesh)
    if gate['closed_oriented_manifold']:
        return mesh, 'literal_closed_source'
    if not allow_conduit_caps:
        return None, 'source_is_not_a_closed_occupied_volume'
    if any(gate[key] for key in ('nonmanifold_edges', 'inconsistent_winding_edges', 'invalid_vertex_links', 'duplicate_faces', 'zero_area_faces')):
        return None, 'invalid_native_conduit_topology'
    vertices, faces = mesh.vertices.copy(), mesh.faces.copy()
    try:
        loops = boundary_loops(faces)
    except ValueError:
        return None, 'unresolved_native_conduit_boundary'
    for loop in loops:
        center = len(vertices)
        vertices = np.r_[vertices, vertices[loop].mean(0)[None]]
        faces = np.r_[faces, np.array([[center, int(loop[(j+1)%len(loop)]), int(loop[j])] for j in range(len(loop))])]
    closed = Mesh(vertices, faces, mesh.frame)
    if not topology(closed)['closed_oriented_manifold']:
        return None, 'temporary_conduit_classifier_is_not_closed'
    return closed, 'temporary_native_conduit_end_caps_never_rendered'


def solid(mesh, origin):
    return make_solid((mesh.vertices-origin)*1000, mesh.faces)[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-version', default='v1')
    parser.add_argument('--route-version', default='v2')
    parser.add_argument('--occupied-volume', action='store_true')
    parser.add_argument('--work-seconds',type=float,default=48.)
    args = parser.parse_args()
    assert 0<args.work_seconds<=48
    began = time.monotonic()
    base = ROOT / '.local-models/geometry-workbench/circulatory-completion'
    directory = base / ('male-left-facial-cervical-finite-network-' + args.network_version)
    out = directory / ('feasibility-v2' if args.occupied_volume else 'feasibility-v1')
    assert not out.exists()
    out.mkdir()
    plan = json.loads((base / ('male-left-facial-cervical-routes-' + args.route_version) / 'routes.json').read_bytes())
    inputs = json.loads((ROOT / plan['source_inputs']).read_bytes())
    wall = Mesh.load(directory / 'finite-wall.npz')
    blood = Mesh.load(directory / 'closed-blood-classifier.npz')
    data = np.load(directory / 'surfaces-and-ancestry.npz')
    ports = json.loads((directory / 'ports.json').read_bytes())
    sources = {r['role']: r for r in inputs['fixed_sources'] if not r.get('missing_native_source')}
    sources.update({r['role']: r for r in plan['additional_fixed_sources']})
    meshes = {role: Mesh.load(ROOT / row['mesh']) for role, row in sources.items()}
    origin = wall.bounds.mean(0)
    occupied = []
    if args.occupied_volume:
        parent_classifier, parent_mode = occupancy_classifier(Mesh.load(ROOT / plan['fixed_parent']['mesh']), True)
        assert parent_classifier is not None, parent_mode
        parent_solid = solid(parent_classifier, origin)
        for route in plan['routes']:
            fields = np.load(directory / (route['name']+'-fields.npz'))
            outer = Mesh(fields['fields'][0], fields['faces'], wall.frame)
            classifier, mode = occupancy_classifier(outer, True)
            assert classifier is not None, (route['name'], mode)
            new_volume = solid(classifier, origin) - parent_solid
            assert new_volume.status() == md.Error.NoError
            occupied.append(new_volume)
    components = []
    for component, route in enumerate(plan['routes'], 1):
        faces = np.concatenate([data['outer_faces'][data['outer_component'] == component], data['inner_faces'][data['inner_component'] == component], data[f'annulus_{component+1}_faces']])
        mesh = Mesh(data['vertices'], faces, wall.frame)
        mesh.save(out / (route['name'] + '-actual-wall.npz'))
        components.append(mesh)
    segments = []
    for component, route in enumerate(plan['routes'], 1):
        points = np.load(ROOT / route['guide'])['points']
        segments.extend((component, i, a, b) for i, (a, b) in enumerate(zip(points[:-1], points[1:])))
    parent_centers = np.array(plan['fixed_parent']['source_station_components'][0]['centers_m'])
    segments.extend((0, i, a, b) for i, (a, b) in enumerate(zip(parent_centers[:-1], parent_centers[1:])))
    midpoints = np.array([(a+b)/2 for _, _, a, b in segments])
    winding = np.abs(igl.winding_number(blood.vertices, blood.faces, midpoints))
    triangles = wall.triangles
    low, high = triangles.min(1), triangles.max(1)
    hits = []
    for component, interval, a, b in segments:
        direction = b-a
        selected = np.flatnonzero(np.all(high >= np.minimum(a,b), axis=1) & np.all(low <= np.maximum(a,b), axis=1))
        tri = triangles[selected]
        e1, e2 = tri[:,1]-tri[:,0], tri[:,2]-tri[:,0]
        h = np.cross(direction, e2)
        det = np.einsum('ij,ij->i', e1, h)
        valid = abs(det) > 1e-24
        selected, tri, e1, e2, h, det = selected[valid], tri[valid], e1[valid], e2[valid], h[valid], det[valid]
        s = a-tri[:,0]
        u = np.einsum('ij,ij->i', s, h)/det
        cross = np.cross(s, e1)
        v = np.einsum('j,ij->i', direction, cross)/det
        amount = np.einsum('ij,ij->i', e2, cross)/det
        failed = np.flatnonzero((u >= 0) & (v >= 0) & (u+v <= 1) & (amount >= 0) & (amount <= 1))
        hits.extend(dict(component=component, interval=interval, face=int(selected[i]), amount=float(amount[i])) for i in failed)
    atomic_npz(out / 'full-guide-segments.npz', segments=np.array([[a,b] for _,_,a,b in segments]), midpoint_winding=winding)
    terminal_rows = []
    for port in ports['ports'][2:]:
        component = port['component']
        route = plan['routes'][component-1]
        tissue = meshes[route['recipient']['recipient_role']]
        ring = wall.vertices[port['outer_indices']]
        center = ring.mean(0)
        radius = np.linalg.norm(ring-center, axis=1).max()
        distance = np.sqrt(igl.point_mesh_squared_distance(center[None], tissue.vertices, tissue.faces)[0][0])
        numbers = abs(igl.winding_number(tissue.vertices, tissue.faces, ring))
        terminal_rows.append(dict(name=route['name'], recipient=route['recipient']['recipient_role'], whole_annular_mouth_inside_tissue=bool(numbers.min() > .99 and distance > radius),
                                  complete_annulus_clearance_lower_m=float(distance-radius), tissue_geometry_sha256=tissue.geometry_sha))
    source = ROOT / plan['source_catalogue']
    assert digest_file(source) == plan['source_catalogue_sha256'], 'Rebase fixed-source plan before context acceptance'
    catalogue = json.loads(source.read_bytes())
    reader = AtlasReader()
    allowed = {
        'ascending-palatine': {'ZA:Palatopharyngeus muscle.l', 'ZA:Soft palate'},
        'tonsillar': {'ZA:Palatine tonsil.l', 'ZA:Superior pharyngeal constrictor.l'},
        'glandular-0': {'ZA:Submandibular gland.l'},
        'glandular-1': {'ZA:Submandibular gland.l'},
    }
    new_bounds = np.array([np.min([m.bounds[0] for m in components], axis=0)-.002, np.max([m.bounds[1] for m in components], axis=0)+.002])
    context_rows, crossings, occupied_crossings, unknown_occupancy = [], [], [], []
    stopped_early = False
    for part in catalogue['parts']:
        pid, bounds = part['id'], part.get('bounds')
        if part.get('suppressed') or bounds is None or not pid.startswith(('ZA:', 'TA98TONGUE:', 'TA98HYOID:', 'TA98ARTERIAL:')):
            continue
        if pid == 'ZA:Facial artery.l' or part.get('system') not in {'skeletal','muscular','nervous','arterial','venous','digestive','respiratory','integumentary'}:
            continue
        if not (np.all(np.asarray(bounds)[1] >= new_bounds[0]) and np.all(np.asarray(bounds)[0] <= new_bounds[1])):
            continue
        target = reader.read_part(catalogue, part, frame=wall.frame, catalogue_path=source, catalogue_sha=plan['source_catalogue_sha256'])
        target_classifier, target_mode = (occupancy_classifier(target, part.get('system') in {'arterial', 'venous'}) if args.occupied_volume else (None, None))
        target_solid = None
        if target_classifier is not None:
            try:
                target_solid = solid(target_classifier, origin)
            except ValueError as error:
                target_mode = 'closed_source_boolean_rejected: '+str(error)
        for ci, (component, route) in enumerate(zip(components, plan['routes'])):
            if time.monotonic()-began > args.work_seconds:
                stopped_early = True
                break
            if not (np.all(target.bounds[1] >= component.bounds[0]) and np.all(target.bounds[0] <= component.bounds[1])):
                continue
            contacts, summary = pair_contacts(component, target, 0., max_pairs=100000, max_tests=1000000)
            assert summary['complete_within_pair_budget'], (pid, route['name'])
            intentional = pid in allowed[route['name']]
            row = dict(part_id=pid, branch=route['name'], target_geometry_sha256=target.geometry_sha, intended_tissue_entry=intentional, summary=summary)
            if args.occupied_volume:
                row['occupied_volume_classifier'] = target_mode
                if target_solid is not None:
                    overlap = occupied[ci] ^ target_solid
                    assert overlap.status() == md.Error.NoError
                    row['new_conduit_occupied_overlap_mm3'] = float(overlap.volume())
                    row['new_conduit_occupied_intersection_empty'] = bool(overlap.is_empty())
                    if not intentional and not overlap.is_empty():
                        occupied_crossings.append(dict(part_id=pid, branch=route['name'], overlap_mm3=float(overlap.volume())))
                elif not intentional:
                    unknown_occupancy.append(dict(part_id=pid, branch=route['name'], reason=target_mode))
            context_rows.append(row)
            if summary['crossing_face_pairs']:
                index = len(context_rows)-1
                atomic_npz(out / (str(index)+'-contacts.npz'), **contacts)
                target.save(out / (str(index)+'-fixed-source.npz'))
                if not intentional:
                    crossings.append(dict(part_id=pid, branch=route['name'], crossing_pairs=summary['crossing_face_pairs']))
        if stopped_early:
            break
    passed = bool(not stopped_early and not hits and winding.min()>.99 and all(t['whole_annular_mouth_inside_tissue'] for t in terminal_rows) and not crossings and not occupied_crossings and not unknown_occupancy)
    result = dict(schema='source-bound-facial-cervical-feasibility/v1', source_catalogue_sha256=plan['source_catalogue_sha256'], candidate_geometry_sha256=wall.geometry_sha,
                  complete_guide_segments=len(segments), wall_guide_intersections=hits, minimum_midpoint_blood_winding=float(winding.min()), terminal_tissue_mouths=terminal_rows,
                  native_context_branch_pairs=len(context_rows), context=context_rows, forbidden_crossings=crossings,
                  occupied_volume_checked=args.occupied_volume, forbidden_occupied_volume_intersections=occupied_crossings, unverified_native_occupied_volumes=unknown_occupancy,
                  occupied_volume_scope='Full new outer conduit occupied volume minus the literal native parent occupied volume; blood plus tissue cannot occupy unrelated solid anatomy. Temporary source-conduit caps never enter rendered geometry.',
                  stopped_before_complete_native_context=stopped_early,
                  intentional_entry_is_not_complete_organ_perfusion=True, context_frame='Only installed native ZA and explicitly installed same-frame tongue, hyoid and lingual families; unregistered donor alternatives cannot certify packing',
                  oral_cavity_closed_volume_available=False, oral_lumen_exclusion_full_proof_pending=True, shared_owner_bank_pending=True,
                  numerical_feasibility_pass=passed, installable=False, installed=False, anatomical_acceptance=False, elapsed_seconds=time.monotonic()-began)
    atomic_json(out / 'evidence.json', result)
    print(json.dumps({k:result[k] for k in ['complete_guide_segments','minimum_midpoint_blood_winding','terminal_tissue_mouths','native_context_branch_pairs','forbidden_crossings','forbidden_occupied_volume_intersections','unverified_native_occupied_volumes','stopped_before_complete_native_context','numerical_feasibility_pass','elapsed_seconds']}))


if __name__ == '__main__':
    main()
