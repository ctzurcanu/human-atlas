"""Original local styloglossus envelopes from explicit native attachment guides.

No proprietary reference vertices, faces or textures are copied. The fixed
styloid cortex supplies a literal origin patch. The intralingual course is a
low-control original construction guided by actual tongue cross-sections.
This is an uninstalled anatomical proposal, not accepted muscle microanatomy.
"""
from pathlib import Path
from collections import defaultdict
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import numpy as np
from scipy.interpolate import CubicSpline

from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.audit import audit, winding_contains
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.solid_attachment import bridge_patches
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames

SOURCE = ROOT / '.local-models/geometry-workbench/tongue-completion/source-v1'
OUT = ROOT / '.local-models/geometry-workbench/tongue-completion/styloglossus-v2'
FRAME = 'model:local-ta98:assembled-metres'
START = time.monotonic()
DEADLINE = START + 40


def remaining():
    if time.monotonic() >= DEADLINE:
        raise TimeoutError('Styloglossus stage exceeded its explicit 40 second native budget')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def unit(v):
    return np.asarray(v, float) / np.linalg.norm(v)


def face_normals(mesh):
    t = mesh.triangles
    return np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])


def source_mesh(row):
    path = Path(row['meshSavePath'])
    assert digest(path) == row['meshSaveSha256']
    result = Mesh.load(path)
    assert result.frame == FRAME
    return result


def origin_patch(bone, ligament_guide, side):
    """An explicit anterolateral styloid patch, retaining exact source faces."""
    bvh = BVH(bone)
    closest = bvh.closest_point(ligament_guide)
    centroids = bone.triangles.mean(1)
    normals = face_normals(bone)
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    distances = np.linalg.norm(centroids-closest['point'], axis=1)
    # The ligament tip establishes which native bony projection is being used.
    # Facing tests then choose its anterior-lateral cortical aspect, not a
    # random closest point on temporal bone or a bounding-box landmark.
    candidates = np.flatnonzero((distances < .004) &
        (side*normals[:, 0] > .05) & (normals[:, 2] > .6))
    if not len(candidates):
        raise ValueError('No witnessed anterolateral styloid cortex within 4 mm of the ligament tip')
    seed = int(candidates[np.argmin(distances[candidates])])
    _, inv = np.unique(bone.vertices, axis=0, return_inverse=True)
    faces = inv[bone.faces]
    edge_faces = defaultdict(list)
    for i, tri in enumerate(faces):
        for j in range(3):
            edge_faces[tuple(sorted((int(tri[j]), int(tri[(j+1) % 3]))))].append(i)
    adjacency = [set() for _ in faces]
    for owners in edge_faces.values():
        for i in owners:
            adjacency[i].update(j for j in owners if j != i)
    chosen, todo = {seed}, [seed]
    while todo:
        i = todo.pop()
        for j in sorted(adjacency[i]):
            if j in chosen or np.linalg.norm(centroids[j]-centroids[seed]) > .00115 \
                    or normals[j]@normals[seed] < .65:
                continue
            chosen.add(j)
            todo.append(j)
    ids = sorted(chosen)
    patch = Mesh(bone.vertices, bone.faces[ids], FRAME)
    # Reject holes or non-disk patches rather than making an invented cap.
    vertices, inverse = np.unique(patch.triangles.reshape(-1, 3), axis=0, return_inverse=True)
    local_faces = inverse.reshape(-1, 3)
    loops = boundary_loops(local_faces)
    if len(loops) != 1:
        raise ValueError('The selected bone footprint is not a single boundary disk')
    area_vectors = face_normals(patch)
    areas = np.linalg.norm(area_vectors, axis=1)
    center = np.einsum('i,ij->j', areas, patch.triangles.mean(1)) / areas.sum()
    return patch, center, dict(
        host_part_id=bone.source['record']['id'], host_geometry_sha256=bone.geometry_sha,
        seed_face=seed, source_faces=ids, ligament_guide=np.asarray(ligament_guide).tolist(),
        ligament_to_bone_distance_m=closest['distance_m'], closest_bone_face=int(closest['face_id']),
        closest_bone_point=closest['point'].tolist(), closest_bone_barycentric=closest['barycentric'].tolist(),
        seed_face_normal=normals[seed].tolist(), seed_face_centroid=centroids[seed].tolist(),
        area_weighted_patch_centroid=center.tolist(), patch_area_m2=float(areas.sum()/2),
        source_face_winding_and_coordinates_literal=True,
        anatomical_choice='Anterolateral native styloid cortex next to the actual stylohyoid ligament superior end; exact histological footprint unverified')


def tongue_surface_guide(tongue, z, side):
    """Exact lower-lateral source contour ray witness at one axial station."""
    segments = []
    for face, tri in enumerate(tongue.triangles):
        points = []
        for i in range(3):
            j = (i+1) % 3
            d0, d1 = tri[i, 2]-z, tri[j, 2]-z
            if d0*d1 < 0:
                points.append(tri[i] + (-d0)/(d1-d0)*(tri[j]-tri[i]))
        if len(points) == 2:
            segments.append((face, np.asarray(points)))
    all_points = np.concatenate([s for _, s in segments])
    center = (all_points[:, :2].min(0)+all_points[:, :2].max(0))*.5
    direction = unit([side*.9, -.6])
    hits = []
    for face, segment in segments:
        matrix = np.column_stack([direction, -(segment[1, :2]-segment[0, :2])])
        if abs(np.linalg.det(matrix)) < 1e-16:
            continue
        distance, segment_parameter = np.linalg.solve(matrix, segment[0, :2]-center)
        if distance > 0 and 0 <= segment_parameter <= 1:
            hits.append((distance, face, np.r_[center+distance*direction, z]))
    if not hits:
        raise ValueError('No actual inferior-lateral contour intersection at the declared station')
    _, face, point = min(hits, key=lambda row: row[0])
    normal = unit(face_normals(tongue)[face])
    witness = BVH(tongue).closest_point(point, seed_face=face)
    assert witness['distance_m'] < 1e-10
    return point, normal, dict(source_face=face, point=point.tolist(), normal=normal.tolist(),
        barycentric=witness['barycentric'].tolist(), station_z_m=z,
        method='actual source-plane intersection, lower-lateral ray; not a box-only target')


def smooth_normals(vertices, faces):
    normals = np.zeros_like(vertices)
    tri = vertices[faces]
    vectors = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    for corner in range(3):
        np.add.at(normals, faces[:, corner], vectors)
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    return np.rint(normals*32767).clip(-32767, 32767).astype('<i2')


def construct(tongue, bone, guide, side, fixed_bones):
    remaining()
    patch, origin, origin_evidence = origin_patch(bone, guide, side)
    witnesses, internal = [], []
    for z in [.0367, .0437, .0507, .0577, .0647, .0717, .0777]:
        surface, normal, evidence = tongue_surface_guide(tongue, z, side)
        inset = float(np.interp(z, [.0367, .0647, .0777], [.0028, .0028, .0012]))
        internal.append(surface-inset*normal)
        witnesses.append(dict(**evidence, inward_offset_m=inset))
    internal = np.asarray(internal)
    departure = unit(internal[0]-origin)
    outward = unit(face_normals(patch).sum(0))
    # Leave the actual cortex into free soft-tissue space before turning
    # medially. A direct chord to tongue would cut through this bony projection.
    departure = unit(departure+outward*max(0., .45-departure@outward))
    origin_evidence['initial_departure_direction'] = departure.tolist()
    origin_evidence['initial_departure_outward_normal_dot'] = float(departure@outward)
    control = np.vstack([origin, origin+departure*.004,
        .55*origin+.45*internal[0]+np.array([-side*.003, .002, 0]), internal])
    parameter = np.r_[0, np.cumsum(np.linalg.norm(np.diff(control, axis=0), axis=1))]
    spline = CubicSpline(parameter, control, axis=0,
        bc_type=((1, departure), (1, unit(control[-1]-control[-2]))))
    stations = np.linspace(0, parameter[-1], 25)
    centers = spline(stations)
    centers[0], centers[-1] = origin, internal[-1]
    t = stations/stations[-1]
    radii = np.column_stack([
        np.interp(t, [0, .18, .65, .90, 1], [.0007, .0011, .00105, .00065, .00025]),
        np.interp(t, [0, .18, .65, .90, 1], [.0010, .0023, .0021, .0011, .00035])])
    frames = rotation_minimizing_frames(centers)
    n, b = frames[-1, :, 0], frames[-1, :, 1]
    angles = 2*np.pi*np.arange(12)/12
    ring = centers[-1]+radii[-1, 0]*np.cos(angles)[:, None]*n+radii[-1, 1]*np.sin(angles)[:, None]*b
    # The distal closure is a construction boundary within the tongue muscle
    # volume, not an invented mucosal or separate histological insertion disk.
    vertices = np.vstack([centers[-1], ring])
    end_faces = np.array([[0, i+1, (i+1) % 12+1] for i in range(12)], dtype=np.int64)
    direction = unit(centers[-1]-centers[-2])
    if np.cross(vertices[1]-vertices[0], vertices[2]-vertices[0])@direction < 0:
        end_faces = end_faces[:, [0, 2, 1]]
    end = Mesh(vertices, end_faces, FRAME)
    positions, surfaces, evidence = bridge_patches(patch, end, centers, radii,
        first_side=-1, second_side=1, section_controls=12)
    faces = np.concatenate(list(surfaces.values())).astype(np.int64)
    # Verify the actual Float32 representation that would be delivered.
    positions = positions.astype('<f4').astype(float)
    result = Mesh(positions, faces, FRAME, {'normals':smooth_normals(positions, faces)}, dict(
        construction='Original slender curved parametric reconstruction; literal Z-Anatomy cortical origin only',
        ta98Term='A05.1.04.105', side='left' if side==1 else 'right',
        source_catalogue_sha256=tongue.source['catalogue_sha256'],
        references=['https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/'],
        anatomy='Styloid origin, descending/anterior course and inferior-lateral intralingual continuation',
        limitation='Exact donor footprint, fibre interdigitation and small vessels/nerves are unverified. No separate source mucosal attachment is claimed.',
        vendor_vertices_faces_or_textures_copied=False))
    numerical = audit(result, max_pairs=25000, max_tests=150000,
        expected_components=1, expected_genus=[0], require_closed=True)
    remaining()
    # Only free wall is tested against cortex: the literal reversed cortical
    # footprint is a declared contact, not an accidental penetration.
    free = Mesh(positions, surfaces['free_wall'], FRAME)
    contacts, summary = pair_contacts(free, bone, 0., max_pairs=20000, max_tests=150000)
    fixed_checks = []
    for host in fixed_bones:
        remaining()
        _, crossing = pair_contacts(result, host, 0., max_pairs=20000, max_tests=150000)
        fixed_checks.append(dict(part_id=host.source['record']['id'], contact=crossing))
    center_inside = winding_contains(tongue, centers, boundary_tolerance_m=1e-8)
    entry = next((i for i, inside in enumerate(center_inside['inside']) if inside), None)
    monotonic_entry = entry is not None and all(center_inside['inside'][entry:])
    return result, dict(origin=origin_evidence, contour_guides=witnesses,
        course_control_points_m=control.tolist(), course_stations_m=centers.tolist(),
        section_semiaxes_m=radii.tolist(), construction=evidence,
        numerical=numerical, free_wall_to_origin_bone=summary,
        other_fixed_bones=fixed_checks,
        centerline_inside_tongue=center_inside, first_intralingual_station=entry,
        centerline_enters_once=monotonic_entry,
        physical_origin_patch_faces_literal=True,
        exact_intralingual_fibre_interfaces_established=False,
        anatomical_acceptance=False)


def main():
    source_audit_path = SOURCE/'source-audit.json'
    source_audit = json.loads(source_audit_path.read_text())
    inventory = ROOT/'reports/ta98-tongue-inventory.json'
    assert inventory.exists() and (ROOT/'TA98-TONGUE.md').exists()
    catalogue = Path(source_audit['sourceCatalogues']['male']['path'])
    assert digest(catalogue) == source_audit['sourceCatalogues']['male']['sha256']
    rows = {row['partId']:row for row in source_audit['contextParts']}
    tongue = source_mesh(rows['ZA:Tongue'])
    fixed_bones = [source_mesh(row) for pid,row in rows.items()
        if pid=='ZA:Mandible:Bone' or pid.startswith('TA98HYOID:')]
    OUT.mkdir(parents=True, exist_ok=True)
    results, render_entries, total = {}, [], 0
    for side, label, suffix in [(1,'left','l'), (-1,'right','r')]:
        bone = source_mesh(rows[f'ZA:Temporal bone.{suffix}:Bone-1'])
        guide = source_audit['guides'][f'styloid-end-{label}']['coordinate']
        mesh, evidence = construct(tongue, bone, guide, side, fixed_bones)
        remaining()
        path = OUT/f'styloglossus-{label}.npz'
        mesh.save(path)
        total += len(mesh.faces)
        results[label] = dict(mesh=str(path), sha256=digest(path), geometry_sha256=mesh.geometry_sha, **evidence)
        render_entries.append(dict(mesh=str(path), name=f'Original styloglossus ({label})', color=[.74,.32,.34,1]))
    assert total <= 2000
    assert digest(catalogue) == source_audit['sourceCatalogues']['male']['sha256']
    review = dict(frame=FRAME, resolution=640, wall_seconds=60,
        meshes=[dict(mesh=rows['ZA:Tongue']['meshSavePath'], name='Native tongue',color=[.81,.52,.55,1]),
                *render_entries,
                dict(mesh=rows['ZA:Soft palate']['meshSavePath'],name='Native soft palate',color=[.66,.49,.67,1])],
        views=[dict(name='styloglossus-front',direction=[0,0,1]),
               dict(name='styloglossus-right',direction=[1,0,0]),
               dict(name='styloglossus-top',direction=[0,1,0],up=[0,0,-1]),
               dict(name='styloglossus-mid-sagittal',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.63,.05],clear_outer=True))])
    (OUT/'review-spec.json').write_text(json.dumps(review,indent=2)+'\n')
    passed = all(r['numerical']['passed'] and r['centerline_enters_once']
        and r['free_wall_to_origin_bone']['complete_within_pair_budget']
        and r['free_wall_to_origin_bone']['crossing_face_pairs']==0
        and all(c['contact']['complete_within_pair_budget'] and
                c['contact']['crossing_face_pairs']==0 for c in r['other_fixed_bones'])
        for r in results.values())
    evidence = dict(schema='original-styloglossus-candidate/v1', source_audit_sha256=digest(source_audit_path),
        inventory_path=str(inventory), inventory_sha256=digest(inventory),
        source_catalogue_path=str(catalogue), source_catalogue_sha256=digest(catalogue),
        sides=results, total_triangles=total, textures_added=0,
        skeleton_and_source_catalogues_modified=False, proprietary_geometry_copied=False,
        numerical_candidate_pass=passed, status='uninstalled candidate; anatomical surface/section review pending',
        elapsed_seconds=time.monotonic()-START)
    (OUT/'candidate.json').write_text(json.dumps(evidence,indent=2)+'\n')
    print(json.dumps(dict(output=str(OUT),total_triangles=total,numerical_candidate_pass=passed,
        elapsed_seconds=evidence['elapsed_seconds'], sides={side:dict(
            geometry_pass=row['numerical']['passed'], failures=row['numerical']['failures'],
            bone_crossings=row['free_wall_to_origin_bone']['crossing_face_pairs'],
            centerline_enters_once=row['centerline_enters_once']) for side,row in results.items()})))


if __name__ == '__main__':
    main()
