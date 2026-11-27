"""Closed computational root/cavity cutters; no fluid is a rendered organ."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl
from scipy.spatial import ConvexHull
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.algorithms.finite_tissue_boolean import make_solid, provenance, classifier_caps_outside
from geometry_tools.algorithms.annular_pipe_envelopes import close_with_shared_disks, boundary_loops

ap = argparse.ArgumentParser()
ap.add_argument('--bank', required=True)
ap.add_argument('--fit', required=True)
ap.add_argument('--stage', required=True)
ap.add_argument('--actual-inner-charts', action='store_true',
                help='Use the exported root inner carrier, including its shared cut charts')
ap.add_argument('--prove-cap-exclusion-in-carve', action='store_true',
                help='Defer touching-cap rejection to complete Boolean output ancestry')
ap.add_argument('--clearance-micrometers', type=float, default=0.,
                help='Computational cutter clearance within the preserved root wall (0 to 20 um)')
args = ap.parse_args()
assert 0. <= args.clearance_micrometers <= 20.
assert not args.clearance_micrometers or args.actual_inner_charts
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
bank = h / args.bank
emission = json.loads((bank / 'emission.json').read_text())
assert digest_file(bank / 'cells.bin') == emission['bank_sha256']
atlas_path = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas = json.loads(atlas_path.read_text())
wall_ids = [pid for pid in atlas['finiteCardiacIntegration']['wallPartIds'] if pid in emission['parts']]
reader = AtlasReader()
mini = dict(chunks=[dict(url='/local-models/' + str((bank / 'cells.bin').relative_to(ROOT / '.local-models')))])
walls = [reader.read_part(mini, dict(emission['parts'][pid], chunk=0, geometryOrigin=[0., 0., 0.]), frame='heart-native-local-metres') for pid in wall_ids]
offset = np.r_[0, np.cumsum([len(m.vertices) for m in walls])]
v, inv = np.unique(np.concatenate([m.vertices for m in walls]), axis=0, return_inverse=True)
f = np.concatenate([inv[m.faces + offset[i]] for i, m in enumerate(walls)])
_, first, count = np.unique(np.sort(f, axis=1), axis=0, return_index=True, return_counts=True)
f = f[first[count == 1]]
host = Mesh(v, f, walls[0].frame)
host.save(out / 'host-boundary.npz')
origin = np.array(emission['geometry_origin'])
placement = json.loads((h / args.fit / 'placement.json').read_text())
if args.actual_inner_charts:
    source_directory = h / emission['source_stage']
    source_metadata = json.loads((source_directory / 'ownership.json').read_text())
    source_arrays = np.load(source_directory / 'arrangement.npz')
    source_positions = source_arrays['positions'].astype('<f4').astype(float)
    surface_definitions = {r['id']: r for r in source_metadata['surfaces']}
records = []
for index, (source, fluid_name, matrix, root_id) in enumerate(zip(
        ['aorta-rounded-course-v276', 'pulmonary-spatial-course-v278'],
        ['fluid-boundary-4.npz', 'fluid-boundary-3.npz'], placement['transforms'],
        ['ZA:Ascending aorta', 'ZA:Pulmonary trunk'])):
    if time.monotonic() - start > 25:
        raise TimeoutError('Outflow computational field construction deadline')
    if args.actual_inner_charts:
        root_uses = source_metadata['cells'][root_id]
        inner_uses = [u for u in root_uses if
                      surface_definitions[u['surface']].get('role') == 'arterial-inner:' + root_id]
        assert inner_uses, ('No explicit inner carrier role', root_id)
        inner_references = [dict(surface=u['surface'], source_face=i)
            for u in inner_uses for i in range(len(source_arrays[surface_definitions[u['surface']]['faces']]))]
        faces = np.concatenate([
            source_arrays[surface_definitions[u['surface']]['faces']][:, [0, 2, 1]]
            if u['side'] > 0 else source_arrays[surface_definitions[u['surface']]['faces']]
            for u in inner_uses])
        used, faces = np.unique(faces, return_inverse=True)
        faces = faces.reshape(-1, 3)
        inner = source_positions[used]
        source = emission['source_stage']
    else:
        native = np.load(h / source / 'source-walls.npz')
        # Literal Float32 source sheet controls used by the tissue bank.
        inner = (native['inner_positions'] + np.array(matrix)[:3, 3] - origin).astype('<f4').astype(float)
        faces = native['source_faces']
    clearance_proof = None
    if args.clearance_micrometers:
        # The cutter may enter the finite root wall by a tiny explicit margin,
        # but that wall is protected and is never cut. This avoids a coincident
        # Boolean chart being mistaken for a new inward-facing host fragment.
        fit_meta = json.loads((h / args.fit / 'ownership.json').read_text())
        root_stage = fit_meta['root_stages'][index]
        native = np.load(h / root_stage / 'source-walls.npz')
        native_inner = native['inner_positions'] + np.array(matrix)[:3, 3] - origin
        native_outer = native['outer_positions'] + np.array(matrix)[:3, 3] - origin
        native_faces = native['source_faces']
        distance2, parents, closest = igl.point_mesh_squared_distance(inner, native_inner, native_faces)
        tri = native_inner[native_faces[parents]]
        bary = igl.barycentric_coordinates(closest, tri[:, 0], tri[:, 1], tri[:, 2])
        delta = np.einsum('ij,ijk->ik', bary, (native_outer-native_inner)[native_faces[parents]])
        width = np.linalg.norm(delta, axis=1)
        assert np.sqrt(distance2.max()) < 5e-8, ('Inner chart/source mismatch', float(np.sqrt(distance2.max())))
        margin = args.clearance_micrometers * 1e-6
        assert width.min() > 50 * margin, ('Insufficient protected wall thickness', float(width.min()))
        inner = inner + delta * (margin / width[:, None])
        clearance_proof = dict(computational_offset_m=margin,
            minimum_corresponding_wall_width_m=float(width.min()),
            maximum_source_projection_error_m=float(np.sqrt(distance2.max())),
            protected_physical_root_unchanged=True, source_stage=root_stage)
    loops = boundary_loops(faces)
    assert len(loops) == 2, ('An arterial inner carrier must have two open rims', root_id, len(loops))
    chamber = Mesh.load(h / 'blood-spaces-v130' / fluid_name)
    chamber = Mesh(chamber.vertices - origin, chamber.faces, host.frame)
    centers = np.array([inner[loop].mean(0) for loop in loops])
    proximal = int(np.argmin(centers[:, 1]))
    mouth = centers[proximal]
    delta = np.array(np.meshgrid(np.linspace(-.025, .025, 17), np.linspace(-.035, -.006, 13), np.linspace(-.025, .025, 17), indexing='ij')).reshape(3, -1).T
    points = mouth + delta
    inside = abs(igl.winding_number(chamber.vertices, chamber.faces, points)) > .99
    points = points[inside]
    assert len(points), 'No ventricular cavity witness in the declared local search bounds'
    distance2 = igl.point_mesh_squared_distance(points, host.vertices, host.faces)[0]
    eligible = np.flatnonzero(distance2 > .003 ** 2)
    assert len(eligible), 'No cavity witness has 3 mm clearance from the fixed tissue'
    order = eligible[np.argsort(np.linalg.norm(points[eligible] - mouth, axis=1))]
    target = None
    for candidate in order[:30]:
        point = points[candidate]
        if abs(igl.winding_number(host.vertices, host.faces, point[None])[0]) < .5:
            target = point
            radius = min(.003, float(np.sqrt(distance2[candidate])) * .45)
            break
    assert target is not None, 'The source blood-space witness lies in actual tissue'
    # A closed convex cavity connector has no flat end disk in the myocardium.
    directions = np.array([q for q in np.ndindex(3, 3, 3) if q != (1, 1, 1)], float) - 1
    directions /= np.linalg.norm(directions, axis=1)[:, None]
    if args.actual_inner_charts:
        # One adjacent triangle band gives the exact exported rim direction.
        rim = loops[proximal]
        adjacent_faces = faces[np.any(np.isin(faces, rim), axis=1)]
        neighbor_points = inner[np.setdiff1d(np.unique(adjacent_faces), rim)]
    else:
        levels = native['station_indices']
        level = int(levels[loops[proximal][0]])
        neighbor = 1 if level == 0 else int(levels.max()) - 1
        neighbor_points = inner[levels == neighbor]
    points = np.r_[inner[loops[proximal]], neighbor_points, target + radius * directions]
    hull = ConvexHull(points)
    hf = hull.simplices.copy()
    tri = points[hf]
    reverse = np.einsum('ti,ti->t', np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), hull.equations[:, :3]) < 0
    hf[reverse] = hf[reverse][:, [0, 2, 1]]
    connector = make_solid(points, hf)[0]
    vv, ff, caps = close_with_shared_disks(inner, faces)
    fluid, source_id = make_solid(vv, ff)
    final = fluid + connector
    mm = final.to_mesh64()
    q, fs = mm.vert_properties[:, :3], np.asarray(mm.tri_verts)
    original, parent = provenance(mm)
    source_caps = np.flatnonzero((original == source_id) & (parent >= len(faces)))
    if args.actual_inner_charts:
        final_references = [inner_references[int(i)] if not args.clearance_micrometers and family == source_id and i < len(faces)
                            else None for family, i in zip(original, parent)]
        atomic_json(out / f'fluid-{index}-source-charts.json', final_references)
    atomic_npz(out / f'fluid-{index}-diagnostic.npz', positions=q, faces=fs,
               cap_face_ids=source_caps, original_inner=inner, original_faces=faces)
    # Every surviving computational disk must be completely outside tissue.
    try:
        cap_proof = classifier_caps_outside(q, fs, source_caps, v, f)
    except ValueError as exc:
        if not args.prove_cap_exclusion_in_carve:
            raise
        cap_proof = dict(strict_separation_passed=False, reason=str(exc),
                         acceptance_pending=True,
                         required_certificate='No Boolean host output face may descend from a cap')
    mesh = Mesh(q, fs, host.frame)
    mesh.save(out / f'fluid-{index}.npz')
    atomic_npz(out / f'fluid-{index}-caps.npz', cap_face_ids=source_caps)
    records.append(dict(source_stage=source, root_part_id=root_id,
        actual_exported_inner_carrier=args.actual_inner_charts,
        target_point_local_m=target.tolist(), target_radius_m=radius,
        proximal_source_mouth=proximal, original_sheet_controls_literal=not bool(args.clearance_micrometers),
        computational_clearance=clearance_proof,
        surviving_computational_caps=len(source_caps), cap_outside_tissue_proof=cap_proof,
        geometry_sha256=mesh.geometry_sha, computational_only=True))
    atomic_json(out / 'progress.json', dict(records=records, elapsed_seconds=time.monotonic() - start), overwrite=True)
    print(source, 'ventricular connector;',
          'cap exclusion pending output ancestry' if cap_proof.get('acceptance_pending')
          else 'caps outside tissue', flush=True)
atomic_json(out / 'fields.json', dict(passed=True, source_bank_sha256=emission['bank_sha256'], wall_part_ids=wall_ids,
    prove_cap_exclusion_by_output_ancestry=args.prove_cap_exclusion_in_carve,
    source_catalogue_sha256=digest_file(atlas_path), records=records,
    computational_only=True, geometry_and_anatomical_acceptance_pending=True, elapsed_seconds=time.monotonic() - start))
