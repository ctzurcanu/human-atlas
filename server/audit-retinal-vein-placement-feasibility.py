"""Stage 2/2: one licensed donor registration and literal outflow-port audit.

This is a feasibility artifact, never an install proposal.  In particular,
contact with a closed donor envelope cannot establish an open venous junction.
"""
from pathlib import Path
import json
import sqlite3
import sys
import time

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, digest_file
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH


START = time.monotonic()
BASE = ROOT / '.local-models/geometry-workbench/circulatory-completion'
SOURCE = BASE / 'male-retinal-vein-source-feasibility-v2/inputs.json'
OUT = BASE / 'male-retinal-vein-placement-feasibility-v1'
OUT.mkdir(exist_ok=False)
data = json.loads(SOURCE.read_text())
INVENTORY = ROOT / 'reports/ta98-blood-circulation-inventory-installed-v8.json'
inventory = json.loads(INVENTORY.read_text())
inventory_terms = {row['id']: row for row in inventory['terms']}
assert inventory['authority_cache_sha256'] == data['authority_sha256']
assert sum(row['male_applicable'] for row in inventory['terms']) == 1149
assert digest_file(data['native_catalogue']) == data['native_catalogue_sha256']
assert digest_file(data['donor_catalogue']) == data['donor_catalogue_sha256']
assert digest_file(data['authority_path']) == data['authority_sha256']
by_role = {(row['side'], row['role']): row for row in data['sources']}


def load(side, role, welded=False):
    row = by_role[side, role]
    assert digest_file(row['mesh']) == row['mesh_sha256']
    return Mesh.load(row['exact_coordinate_welded_analysis'] if welded else row['mesh'])


def sphere_landmarks(sclera, nerve):
    """Fit actual scleral samples; preserve residuals and source witnesses."""
    points = np.unique(sclera.vertices, axis=0)
    origin = points.mean(axis=0)
    local = points - origin
    solution = np.linalg.lstsq(
        np.c_[2 * local, np.ones(len(local))],
        np.einsum('ij,ij->i', local, local), rcond=None)[0]
    center = origin + solution[:3]
    distances = np.linalg.norm(points - center, axis=1)
    radius = np.sqrt(solution[3] + np.dot(solution[:3], solution[:3]))
    nerve_points = np.unique(nerve.vertices, axis=0)
    # This is a witnessed placement cue, not a named optic-disc boundary.
    count = max(8, int(np.ceil(len(nerve_points) * .06)))
    ids = np.argsort(np.linalg.norm(nerve_points - center, axis=1), kind='stable')[:count]
    entry = nerve_points[ids].mean(axis=0)
    posterior = entry - center
    posterior /= np.linalg.norm(posterior)
    up = np.array([0., 1., 0.])
    up -= np.dot(up, posterior) * posterior
    up /= np.linalg.norm(up)
    right = np.cross(up, posterior)
    frame = np.column_stack([right, up, posterior])
    assert np.allclose(frame.T @ frame, np.eye(3), atol=1e-12)
    assert np.linalg.det(frame) > .999999999
    return dict(center=center, radius=radius, frame=frame,
                radial_fit_residual_rms_m=np.sqrt(np.mean((distances-radius)**2)),
                radial_fit_residual_max_m=np.max(np.abs(distances-radius)),
                nerve_entry_patch_center=entry,
                nerve_entry_unique_vertex_ids=ids,
                nerve_entry_source_vertices=nerve_points[ids],
                qualification='Least-squares scleral globe and nearest six-percent optic-nerve patch; not an anatomical registration acceptance or named optic-disc boundary')


def exact_shared_positions(first, second):
    a = np.unique(first.vertices, axis=0)
    b = np.unique(second.vertices, axis=0)
    dtype = np.dtype((np.void, 3 * np.dtype('<f8').itemsize))
    return int(len(np.intersect1d(np.ascontiguousarray(a).view(dtype).ravel(),
                                 np.ascontiguousarray(b).view(dtype).ravel())))


def actual_mouth_distances(candidate, host, host_row):
    result = []
    bvh = BVH(candidate)
    for number, loop in enumerate(host_row['boundary_polygons']):
        vertices = host.vertices[loop['indices']]
        queries, proof = bvh.closest_points(vertices, backend='bvh')
        center = vertices.mean(axis=0)
        result.append(dict(
            actual_host_boundary_number=number,
            actual_host_boundary_vertex_ids=loop['indices'],
            actual_host_boundary_vertices_m=vertices,
            actual_host_boundary_center_m=center,
            boundary_to_closed_donor_surface_min_m=queries['distance_m'].min(),
            boundary_to_closed_donor_surface_max_m=queries['distance_m'].max(),
            per_vertex_source_face_witnesses=queries,
            center_to_closed_donor_surface=bvh.closest_point(center),
            proof=proof,
            matched_donor_lumen_mouth=False,
            qualification='Distances are to literal closed envelope triangles, not proof of a communicating lumen or a paired annulus'))
    return result


# Parent closure includes these A15 codes in the existing 1,149 male blood
# terms. Code prefixes must never be used as an inventory-membership test.
# Equal FMA/source identities are recorded instead of silently double-counted.
scope = []
with sqlite3.connect('file:' + data['authority_path'] + '?mode=ro', uri=True) as db:
    db.row_factory = sqlite3.Row
    for term in data['terms']:
        row = term['authority']
        actual_inventory_row = inventory_terms[row['id']]
        assert actual_inventory_row['parent_id'] == row['parent_id']
        assert actual_inventory_row['scope'] == 'blood-circulation'
        assert actual_inventory_row['male_applicable']
        same_source = [dict(other) for other in db.execute(
            'SELECT id,name_en,parent_id,source_id,fma_id FROM ta98 WHERE source_id=? AND id<>? ORDER BY id',
            (row['source_id'], row['id']))]
        same_fma = [dict(other) for other in db.execute(
            "SELECT id,name_en,parent_id,source_id,fma_id FROM ta98 WHERE fma_id=? AND fma_id IS NOT NULL AND fma_id<>'' AND id<>? ORDER BY id",
            (row['fma_id'], row['id']))]
        scope.append(dict(
            id=row['id'], name=row['name_en'], parent_id=row['parent_id'],
            source_id=row['source_id'], fma_id=row['fma_id'],
            representation='small venous vessel' if row['id'].startswith('A15.') else term['representation_type'],
            authority_code_prefix=row['id'].split('.')[0],
            inside_existing_1149_male_blood_closure=True,
            additional_to_existing_denominator=False,
            inventory_row=actual_inventory_row,
            same_source_identity_other_authority_records=same_source,
            same_fma_identity_other_authority_records=same_fma,
            duplicates_canonical_A12_identity=any(other['id'].startswith('A12.') for other in same_source + same_fma),
            direct_distinct_geometry_available=False, represented_or_installed=False))
atomic_json(OUT / 'scope-and-authority.json', dict(
    schema='retinal-cross-chapter-blood-scope/v2', male_blood_parent_closure_denominator=1149,
    denominator_follows_authority_parent_closure_across_code_chapters=True,
    A15_retinal_venules_already_inside_existing_denominator=7,
    additional_terms_to_add_to_denominator=0,
    A15_retinal_venule_status='Missing distinct geometry; included in the existing denominator, neither available nor installed',
    inventory_path=str(INVENTORY), inventory_sha256=digest_file(INVENTORY),
    inventory_source_catalogue_sha256=inventory['source_catalogue_sha256'],
    inventory_scope_policy=inventory['scope_policy'],
    authority_sha256=data['authority_sha256'], terms=scope))

sides = []
for side in ['left', 'right']:
    donor_sclera = load(side, 'donor-sclera', True)
    donor_nerve = load(side, 'donor-optic-nerve', True)
    native_sclera = load(side, 'native-sclera', True)
    native_nerve = load(side, 'native-optic-nerve', True)
    donor = sphere_landmarks(donor_sclera, donor_nerve)
    native = sphere_landmarks(native_sclera, native_nerve)
    rotation = native['frame'] @ donor['frame'].T
    scale = native['radius'] / donor['radius']
    translation = native['center'] - scale * rotation @ donor['center']
    transform = np.eye(4)
    transform[:3, :3] = scale * rotation
    transform[:3, 3] = translation
    assert scale > 0 and np.linalg.det(rotation) > .999999999
    vessel = load(side, 'donor-vein', True)
    candidate = Mesh(vessel.vertices @ rotation.T * scale + translation,
                     vessel.faces, by_role[side, 'native-sclera']['frame'],
                     source=dict(kind='licensed-donor-feasibility-only',
                                 source_mesh=by_role[side, 'donor-vein']['mesh'],
                                 source_mesh_sha256=by_role[side, 'donor-vein']['mesh_sha256'],
                                 positive_similarity=transform.tolist(),
                                 anatomical_registration_verified=False))
    path = OUT / (side + '-provisional-registered-outer-envelope.npz')
    candidate.save(path)
    candidate_bvh = BVH(candidate)
    donor_host = load(side, 'donor-receiving-vein', True)
    native_host = load(side, 'native-receiving-vein', True)
    source_contact = BVH(vessel).closest_mesh(BVH(donor_host))
    native_contact = candidate_bvh.closest_mesh(BVH(native_host))
    # Sampled point-to-triangle distances qualify placement only.  They are
    # not an exhaustive intersection/containment acceptance of a new vessel.
    sampled_context = []
    for role in ['native-retina', 'native-sclera', 'native-optic-nerve']:
        context = load(side, role, True)
        samples = candidate.vertices[np.linspace(0, len(candidate.vertices)-1, 24, dtype=int)]
        distances, proof = BVH(context).closest_points(samples, backend='bvh')
        sampled_context.append(dict(role=role, samples_m=samples, witnesses=distances,
                                    proof=proof, exhaustive_crossing_check=False))
    artery = load(side, 'native-retinal-artery')
    sides.append(dict(
        side=side, donor_landmarks=donor, native_landmarks=native,
        proper_rotation=rotation, positive_uniform_scale=scale,
        rotation_degrees=np.degrees(np.arccos(np.clip((np.trace(rotation)-1)/2, -1, 1))),
        transform=transform, candidate_mesh=str(path), candidate_mesh_sha256=digest_file(path),
        source_vertex_and_face_order_preserved=True,
        donor_source_contact_with_own_SOV=source_contact,
        donor_exact_shared_positions_with_own_SOV=exact_shared_positions(vessel, donor_host),
        candidate_contact_with_native_SOV=native_contact,
        donor_source_boundary_loops=by_role[side, 'donor-vein']['boundary_polygons'],
        native_host_actual_mouths=actual_mouth_distances(candidate, native_host, by_role[side, 'native-receiving-vein']),
        sampled_fixed_eye_context=sampled_context,
        native_retinal_artery_fixed=dict(part_id=by_role[side, 'native-retinal-artery']['part_id'],
                                       bounds_m=artery.bounds, exhaustive_crossing_check=False),
        represented_receiving_host=by_role[side, 'native-receiving-vein']['part_id'],
        actual_open_donor_mouth_available=False, finite_inner_wall_available=False,
        intraocular_extraocular_partition_source_witnessed=False,
        install_eligible=False,
        obstacle='Closed donor outer envelope has no actual open outflow annulus or inner wall to pair with either native SOV boundary. Nearest contact cannot establish a patent venous junction.'))

assert digest_file(data['native_catalogue']) == data['native_catalogue_sha256']
assert digest_file(data['donor_catalogue']) == data['donor_catalogue_sha256']
evidence = dict(
    schema='retinal-vein-two-stage-feasibility/v1', stage1_inputs=str(SOURCE),
    stage1_inputs_sha256=digest_file(SOURCE), scope_file=str(OUT/'scope-and-authority.json'),
    scope_sha256=digest_file(OUT/'scope-and-authority.json'), sides=sides,
    inventory_path=str(INVENTORY), inventory_sha256=digest_file(INVENTORY),
    scope_note='All seven A15 retinal venules are already included in the existing 1,149 male blood parent-closure terms and remain unrepresented.',
    fixed_catalogue_sha256=data['native_catalogue_sha256'],
    donor_catalogue_sha256=data['donor_catalogue_sha256'],
    native_geometry_or_catalogue_changed=False, existing_heart_eye_skeleton_vessels_fixed=True,
    registration_is_provisional=True, installation_eligible=False,
    stopped_at_interface_evidence=True,
    next_defensible_stage='If separately authorized, construct an independently modeled finite central retinal vein wall with source-calibre optic-nerve passage and a real cut paired SOV side-ostium; do not insert a closed donor tip or count the seven absent retinal venules.',
    elapsed_seconds=time.monotonic()-START)
atomic_json(OUT/'evidence.json', evidence)
print(json.dumps(dict(output=str(OUT), sides=[dict(
    side=row['side'], scale=row['positive_uniform_scale'], rotation_deg=row['rotation_degrees'],
    donor_SOV_surface_gap_mm=float(row['donor_source_contact_with_own_SOV']['distance_m']*1000),
    registered_native_SOV_surface_gap_mm=float(row['candidate_contact_with_native_SOV']['distance_m']*1000),
    native_mouth_to_donor_surface_min_mm=[float(x['boundary_to_closed_donor_surface_min_m']*1000) for x in row['native_host_actual_mouths']])for row in sides],
    A15_missing=7, install_eligible=False, elapsed_seconds=evidence['elapsed_seconds'])))
