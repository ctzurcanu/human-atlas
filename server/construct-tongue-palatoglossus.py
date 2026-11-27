"""Bounded original palatoglossus proposal with witnessed native attachments.

Run `--stage plan` before `--stage construct`. The adult primary reference
establishes a flattened faucial belly, fan at the palate and tapered lingual
termination, but supplies no measured dimensions for this particular atlas.
All chosen dimensions remain explicit modeling estimates. No source vendor
muscle, texture or UV data is copied, and this script never installs anatomy.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'server'))

import numpy as np

from geometry_tools.algorithms.tube_sections import (
    centripetal_curve, resample_path, rotation_minimizing_frames)
from geometry_tools.audit import audit, winding_contains
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.contact import pair_contacts
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH
from tongue_whole_input import CATALOGUE, SOURCE_AUDIT, WORK, fixed_bone_sources, whole_input

FRAME = 'model:local-ta98:assembled-metres'
REFERENCE = 'https://pubmed.ncbi.nlm.nih.gov/281278/'
WHOLE = WORK / 'basal-whole-protected-v1/whole-tongue.npz'
FIT = WORK / 'basal-whole-protected-v1/fit-result.json'


def require(value, message):
    if not value:
        raise ValueError(message)


def unit(value):
    value = np.asarray(value, float)
    return value / np.linalg.norm(value, axis=-1, keepdims=True)


def bound(path):
    path = Path(path).resolve()
    return dict(path=str(path), sha256=digest_file(path))


def load_bound(row):
    path = Path(row['path'])
    require(path.is_relative_to(WORK) and digest_file(path) == row['sha256'],
            'Changed or nonlocal construction input')
    return Mesh.load(path)


def oral_palate_origin(mesh, x, z):
    """Exact vertical intersection with the outward oral palatal underside."""
    hits = []
    for face, tri in enumerate(mesh.triangles):
        edges = (tri[1:] - tri[0])[:, [0, 2]].T
        if abs(np.linalg.det(edges)) < 1e-18:
            continue
        uv = np.linalg.solve(edges, np.array([x, z]) - tri[0, [0, 2]])
        weights = np.r_[1-uv.sum(), uv]
        if weights.min() < -1e-12:
            continue
        normal = unit(np.cross(tri[1]-tri[0], tri[2]-tri[0]))
        if normal[1] >= -.1:
            continue
        point = weights @ tri
        hits.append((point[1], face, point, weights, normal))
    require(hits, 'No actual oral palatal underside at the declared fan station')
    _, face, point, weights, normal = min(hits, key=lambda row: row[0])
    return dict(part_id='ZA:Soft palate', face=int(face), point=point.tolist(),
                barycentric=weights.tolist(), outward_normal=normal.tolist(),
                geometry_sha256=mesh.geometry_sha,
                method='literal vertical-ray/source-triangle intersection; outward oral underside')


def lingual_insertion(mesh, side):
    target = np.array([side*.018, 1.633, .048])
    witness = BVH(mesh).closest_point(target)
    tri = mesh.triangles[witness['face_id']]
    normal = unit(np.cross(tri[1]-tri[0], tri[2]-tri[0]))
    require(side*normal[0] > .25,
            'The proposed endpoint did not witness the actual lateral tongue')
    return dict(part_id='ZA:Tongue', face=int(witness['face_id']),
                point=witness['point'].tolist(), barycentric=witness['barycentric'].tolist(),
                outward_normal=normal.tolist(), geometry_sha256=mesh.geometry_sha,
                guide_point=target.tolist(), guide_distance_m=witness['distance_m'],
                method='closest actual lateral facet to an explicitly planned posterior-tongue guide')


def palatal_depth(mesh, witness):
    """First actual opposite tissue intersection along the inward source normal."""
    point = np.asarray(witness['point'])
    direction = -np.asarray(witness['outward_normal'])
    hits = []
    for face, tri in enumerate(mesh.triangles):
        matrix = np.column_stack([direction, -(tri[1]-tri[0]), -(tri[2]-tri[0])])
        if abs(np.linalg.det(matrix)) < 1e-20:
            continue
        distance, u, v = np.linalg.solve(matrix, tri[0]-point)
        if distance > 1e-7 and u >= -1e-10 and v >= -1e-10 and u+v <= 1+1e-10:
            hits.append((distance, face))
    require(hits, 'No finite opposite palatal tissue boundary; cannot invent origin thickness')
    depth, opposite = min(hits)
    center = point+.3*depth*direction
    require(winding_contains(mesh, center[None], boundary_tolerance_m=1e-8)['inside'][0],
            'Actual measured palatal material interval has no interior origin witness')
    return float(depth), int(opposite), center


def plan(output):
    require(output.is_relative_to(WORK) and not output.exists(), 'Use a new local proposal directory')
    tongue, tongue_path, binding = whole_input(WHOLE, FIT)
    source = json.loads(SOURCE_AUDIT.read_text())
    rows = {row['partId']: row for row in source['contextParts']}
    palate_row = rows['ZA:Soft palate']
    palate = Mesh.load(palate_row['meshSavePath'])
    require(digest_file(palate_row['meshSavePath']) == palate_row['meshSaveSha256'], 'Palate source changed')
    output.mkdir(parents=True)
    atlas = json.loads(CATALOGUE.read_text())
    by_id = {part['id']: part for part in atlas['parts']}
    reader = AtlasReader()
    contexts = []
    ids = ['ZA:Palatine tonsil.l', 'ZA:Palatine tonsil.r',
           'ZA:Palatopharyngeus muscle.l', 'ZA:Palatopharyngeus muscle.r']
    for index, part_id in enumerate(ids):
        require(part_id in by_id, 'The explicit adjoining context is absent: '+part_id)
        mesh = reader.read_part(atlas, by_id[part_id], frame=FRAME,
                                catalogue_path=CATALOGUE, catalogue_sha=binding['source_catalogue_sha256'])
        path = output / f'fixed-adjoining-{index}.npz'
        mesh.save(path)
        contexts.append(dict(part_id=part_id, geometry_sha256=mesh.geometry_sha, **bound(path)))
    sides = {}
    for side, label in [(1, 'left'), (-1, 'right')]:
        origin = oral_palate_origin(palate, side*.013, .040)
        insertion = lingual_insertion(tongue, side)
        op, on = np.array(origin['point']), np.array(origin['outward_normal'])
        ip, inn = np.array(insertion['point']), np.array(insertion['outward_normal'])
        depth, opposite, center = palatal_depth(palate, origin)
        # One cubic retains the actual fan's normal departure. Piecewise
        # interpolation through an arbitrary lateral waypoint made the wide
        # native roof fan cross its first two stations in v2.
        end = ip-.002*inn
        control = np.vstack([center, op+.009*on,
                             ip+.006*inn+np.array([0,.007,0]), end])
        sides[label] = dict(side=side, palate_origin=origin, tongue_insertion=insertion,
                            control_points_m=control.tolist(),
                            palate_origin_inset_m=.3*depth, tongue_insertion_inset_m=.002,
                            measured_palatal_normal_depth_m=depth, opposite_palatal_face=opposite)
    report = dict(schema='tongue-palatoglossus-construction-plan/v1',
        inventory=bound(ROOT/'reports/ta98-tongue-inventory.json'),
        authority_term='A05.1.04.110', canonical_parent='A04.1.06.001',
        whole_organ='Tongue, male first; bilateral palatoglossus in the native oral assembly',
        representation='paired finite macroscopic muscle envelopes; not exclusive microscopic fibre ownership',
        reference=dict(url=REFERENCE, study='1978 adult cadaver investigation, 25 specimens',
                       shape='flattened faucial belly; fan-shaped palatal attachment; tapered lingual termination',
                       exact_source_dimensions_established=False),
        source_catalogue=bound(CATALOGUE), source_audit=bound(SOURCE_AUDIT),
        whole=bound(tongue_path), whole_input=binding,
        soft_palate=bound(palate_row['meshSavePath']),
        fixed_bones=[dict(part_id=row['partId'], **bound(row['meshSavePath'])) for row in fixed_bone_sources()],
        fixed_adjoining=contexts, sides=sides,
        dimensions=dict(modeling_estimates=True, fan_full_width_m=.0054,
                        belly_full_width_m=.0032, belly_full_thickness_m=.0008,
                        distal_full_width_m=.0007, distal_full_thickness_m=.0003),
        adjoining_interfaces=[
            dict(host='native soft palate', relation='palatal aponeurotic origin; fan remains within/continuous with palate',
                 endpoint_inside_required=True, microscopic_shared_fibres_established=False),
            dict(host='corrected native lateral posterior tongue', relation='intralingual tapered insertion',
                 endpoint_inside_required=True, muscle_envelope_overlap_qualified=True),
            dict(host='faucial-pillar mucosa', relation='covering must be constructed/reviewed before completed mucosal assembly',
                 represented=False, pending=True),
            dict(host='mandible, temporal bones and hyoid', relation='fixed disjoint bone context', movement_authorized=False),
            dict(host='palatine tonsils and palatopharyngeus', relation='fixed adjoining tissue; disjoint course required', movement_authorized=False)],
        sections=25, ring_controls=12, max_triangles_per_side=600,
        max_seconds=45, source_and_skeleton_modified=False,
        course_algorithm='single-cubic-bezier-with-witnessed-outward-palatal-tangent',
        anatomical_acceptance=False, install_authorized_by_this_script=False)
    atomic_json(output/'plan.json', report)
    print(json.dumps(dict(plan=str(output/'plan.json'), sides=sides)))


def normals(vertices, faces):
    tri = vertices[faces]
    area = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    values = np.zeros_like(vertices)
    for corner in range(3):
        np.add.at(values, faces[:, corner], area)
    require(np.all(np.linalg.norm(values, axis=1) > 0), 'Degenerate averaged normal')
    return np.rint(unit(values)*32767).clip(-32767, 32767).astype('<i2')


def normal_charts(vertices, faces):
    """Preserve every physical triangle and give sharp tissue rims local charts."""
    tri = vertices[faces]
    areas = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
    directions = unit(areas)
    owners = [[] for _ in vertices]
    for face, ids in enumerate(faces):
        for vertex in ids:
            owners[vertex].append(face)
    words = np.empty((*faces.shape, 3), dtype='<i2')
    for face, ids in enumerate(faces):
        for corner, vertex in enumerate(ids):
            neighbors = np.asarray(owners[vertex])
            nearby = neighbors[directions[neighbors]@directions[face] >= np.cos(np.deg2rad(45))]
            words[face, corner] = np.rint(unit(areas[nearby].sum(0))*32767).astype('<i2')
    keys = np.c_[faces.reshape(-1), words.reshape(-1,3)].astype(np.int64)
    charts, inverse = np.unique(keys, axis=0, return_inverse=True)
    positions, indices = vertices[charts[:,0]], inverse.reshape(-1,3)
    require(np.array_equal(positions[indices], vertices[faces]), 'A normal chart changed physical triangles')
    return positions, indices, charts[:,1:].astype('<i2')


def construct(output):
    began = time.monotonic()
    data = json.loads((output/'plan.json').read_text())
    require(data['schema'] == 'tongue-palatoglossus-construction-plan/v1', 'Unknown explicit plan')
    require(not (output/'candidate.json').exists(), 'Candidate geometry is immutable; use a new planned stage')
    require(digest_file(CATALOGUE) == data['source_catalogue']['sha256'], 'Installed source catalogue changed')
    tongue, palate = load_bound(data['whole']), load_bound(data['soft_palate'])
    bones = [(row['part_id'], load_bound(row)) for row in data['fixed_bones']]
    adjoining = [(row['part_id'], load_bound(row)) for row in data['fixed_adjoining']]
    results, meshes = {}, []
    for label, side_row in data['sides'].items():
        control = np.asarray(side_row['control_points_m'])
        sample = np.linspace(0, 1, 257)[:,None]
        polyline = (1-sample)**3*control[0]+3*(1-sample)**2*sample*control[1] \
                   +3*(1-sample)*sample**2*control[2]+sample**3*control[3]
        tangents = unit(3*(1-sample)**2*(control[1]-control[0]) \
                       +6*(1-sample)*sample*(control[2]-control[1]) \
                       +3*sample**2*(control[3]-control[2]))
        centers, tangents = resample_path(polyline, data['sections'], tangents)
        frames = rotation_minimizing_frames(centers, tangents, initial_reference=[1, 0, 0])
        t = np.linspace(0, 1, len(centers))
        # Tangential fan -> flattened faucial ribbon, using a smooth first-quarter roll.
        blend = np.clip((t-.2)/.5, 0, 1)
        roll = .5*np.pi*blend*blend*(3-2*blend)
        n, b = frames[:, :, 0].copy(), frames[:, :, 1].copy()
        major = np.cos(roll)[:, None]*n + np.sin(roll)[:, None]*b
        minor = -np.sin(roll)[:, None]*n + np.cos(roll)[:, None]*b
        a = np.interp(t, [0, .2, .75, 1], [.0027, .0016, .0016, .00035])
        b = np.interp(t, [0, .2, .75, 1], [.00035, .0004, .0004, .00015])
        angles = 2*np.pi*np.arange(data['ring_controls'])/data['ring_controls']
        rings = centers[:, None, :] + a[:, None, None]*np.cos(angles)[None, :, None]*major[:, None, :] \
                + b[:, None, None]*np.sin(angles)[None, :, None]*minor[:, None, :]
        # The fan starts INSIDE the actual finite palate, not .8 mm beyond a
        # thin source roof. Each fan control has its own material-ray witness.
        fan_witnesses = []
        origin_point = np.asarray(side_row['palate_origin']['point'])
        for corner, angle in enumerate(angles):
            witness = oral_palate_origin(palate, origin_point[0]+a[0]*np.cos(angle),
                                        origin_point[2]+b[0]*np.sin(angle))
            depth, opposite, position = palatal_depth(palate, witness)
            rings[0, corner] = position
            fan_witnesses.append(dict(**witness, opposite_face=opposite,
                                      normal_material_depth_m=depth, interior_point=position.tolist()))
        vertices = np.vstack([rings.reshape(-1, 3), centers[0], centers[-1]])
        count, sides = len(centers), data['ring_controls']
        faces = []
        for station in range(count-1):
            for corner in range(sides):
                i, j = station*sides+corner, station*sides+(corner+1)%sides
                k, ell = i+sides, j+sides
                faces.extend([[i, j, ell], [i, ell, k]])
        for corner in range(sides):
            faces.append([count*sides, (corner+1)%sides, corner])
            faces.append([count*sides+1, (count-1)*sides+corner, (count-1)*sides+(corner+1)%sides])
        faces = np.asarray(faces, np.int64)
        # Test exactly the intended dyadic-origin Float32 runtime coordinates.
        origin = np.array([0, 1.625, 0])
        vertices = (vertices-origin).astype('<f4').astype(float)+origin
        fan_inside = winding_contains(palate, vertices[:sides], boundary_tolerance_m=1e-8)
        positions, indices, words = normal_charts(vertices, faces)
        mesh = Mesh(positions, indices, FRAME, {'normals':words},
                    dict(ta98Term='A05.1.04.110', side=label,
                         construction='original parameterized flattened fan-ribbon with witnessed native endpoint facets',
                         source_catalogue_sha256=data['source_catalogue']['sha256'],
                         reference=REFERENCE, donor_dimensions=False, mucosal_pillar_pending=True,
                         source_vendor_geometry_or_texture_copied=False,
                         geometry_origin_for_exact_float32=origin.tolist()))
        path = output/f'palatoglossus-{label}.npz'
        mesh.save(path)
        numerical = audit(mesh, max_pairs=30000, max_tests=200000,
                          expected_components=1, expected_genus=[0], require_closed=True)
        checks = []
        for part_id, context in bones+adjoining:
            require(time.monotonic()-began < data['max_seconds'], 'Native construction deadline reached')
            arrays, summary = pair_contacts(mesh, context, 0, max_pairs=30000, max_tests=200000)
            witness = output/(label+'-'+part_id.replace(':', '_').replace(' ', '_')+'.npz')
            atomic_npz(witness, **arrays)
            checks.append(dict(part_id=part_id, summary=summary,
                               witness=str(witness), witness_sha256=digest_file(witness)))
        inside_tongue = winding_contains(tongue, centers, boundary_tolerance_m=1e-8)
        inside_palate = winding_contains(palate, centers, boundary_tolerance_m=1e-8)
        tongue_states = inside_tongue['inside']
        palate_states = inside_palate['inside']
        entry = next((i for i, yes in enumerate(tongue_states) if yes), None)
        exit_palate = next((i for i, yes in enumerate(palate_states) if not yes), None)
        enters_once = entry is not None and all(tongue_states[entry:])
        leaves_once = palate_states[0] and exit_palate is not None and not any(palate_states[exit_palate:])
        actual_normals = mesh.attributes['normals'][mesh.faces].astype(float)
        tri = mesh.triangles
        area = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
        minimum_dot = float(np.min(np.einsum('fci,fi->fc', unit(actual_normals), unit(area))))
        result = dict(mesh=str(path), sha256=digest_file(path), geometry_sha256=mesh.geometry_sha,
            numerical=numerical, fixed_context=checks, inside_tongue=inside_tongue,
            inside_palate=inside_palate, centerline_enters_tongue_once=bool(enters_once),
            centerline_leaves_palate_once=bool(leaves_once), first_tongue_station=entry,
            first_free_palate_station=exit_palate, minimum_corner_face_normal_dot=minimum_dot,
            palatal_fan_material_witnesses=fan_witnesses, palatal_fan_controls_inside=fan_inside,
            centerline_m=centers.tolist(), section_semiaxes_m=np.c_[a,b].tolist(),
            physical_triangles=len(faces), approximate_course_length_m=float(np.linalg.norm(np.diff(centers,axis=0),axis=1).sum()),
            estimated_packed_bytes=len(positions)*18+len(faces)*12)
        results[label] = result
        meshes.append(mesh)
    arrays, paired = pair_contacts(*meshes, 0, max_pairs=30000, max_tests=200000)
    atomic_npz(output/'paired-witnesses.npz', **arrays)
    passed = all(row['numerical']['passed'] and row['minimum_corner_face_normal_dot'] > 0 and
        row['centerline_enters_tongue_once'] and row['centerline_leaves_palate_once'] and
        all(row['palatal_fan_controls_inside']['inside']) and
        all(c['summary']['complete_within_pair_budget'] and
            not c['summary']['crossing_face_pairs'] and not c['summary']['coplanar_overlap_face_pairs']
            for c in row['fixed_context']) for row in results.values()) and paired['complete_within_pair_budget'] \
            and not paired['crossing_face_pairs'] and not paired['coplanar_overlap_face_pairs']
    candidate = dict(schema='tongue-palatoglossus-candidate/v1', plan=bound(output/'plan.json'),
                     sides=results, paired=paired, numerical_candidate_pass=bool(passed),
                     anatomical_acceptance=False, installed=False,
                     mucosal_pillar_continuity_pending=True,
                     elapsed_seconds=time.monotonic()-began)
    atomic_json(output/'candidate.json', candidate)
    spec = dict(frame=FRAME, resolution=640, wall_seconds=60,
                meshes=[dict(mesh=data['whole']['path'], name='Corrected whole tongue', color=[.8,.48,.5,1]),
                        dict(mesh=data['soft_palate']['path'], name='Native soft palate', color=[.66,.49,.67,1]),
                        *[dict(mesh=row['mesh'], name='Palatoglossus '+label, color=[.72,.25,.29,1]) for label,row in results.items()]],
                views=[dict(name='palatoglossus-front', direction=[0,0,1]),
                       dict(name='palatoglossus-right', direction=[1,0,0]),
                       dict(name='palatoglossus-posterior', direction=[0,0,-1]),
                       dict(name='palatoglossus-sagittal', direction=[1,0,0],
                            section=dict(normal=[1,0,0], point=[0,1.64,.045], clear_outer=True))])
    atomic_json(output/'review-spec.json', spec)
    print(json.dumps(dict(output=str(output), numerical_candidate_pass=passed,
        elapsed_seconds=candidate['elapsed_seconds'], sides={label:dict(
            audit=row['numerical']['passed'], normal_dot=row['minimum_corner_face_normal_dot'],
            context_crossings=sum(c['summary']['crossing_face_pairs'] for c in row['fixed_context']),
            tongue_entry=row['centerline_enters_tongue_once'], palate_exit=row['centerline_leaves_palate_once'])
            for label,row in results.items()})))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stage', choices=['plan', 'construct'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    (plan if args.stage == 'plan' else construct)(args.output.resolve())
