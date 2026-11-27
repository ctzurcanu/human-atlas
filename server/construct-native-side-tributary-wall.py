"""Source-bound side-tributary finite-wall producer; no atlas installation.

Stages: prepare, outer, inner, assemble. Each stage must be supervised by
geometry_tools.runtime (normally 60 s). Source NPZs contain vertices, faces,
profiles and profile_vertex_indices. Profiles index EVERY native vertex once.
The supplied plan and its SHA256 are mandatory; this program invents neither
branch identity nor anatomical permission. Native controls are not deformed.
"""
from pathlib import Path
import argparse
import json
import sys
import time

import manifold3d as md
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.finite_tissue_boolean import (  # noqa: E402
    CommonPositionBank, make_solid, oriented_key, provenance,
)
from geometry_tools.algorithms.tube_continuity import (  # noqa: E402
    annular_rim, boundary_loops, continuity_audit,
)
from geometry_tools.algorithms.tube_sections import (  # noqa: E402
    arclength, bend_diagnostics, rotation_minimizing_frames, unit,
)
from geometry_tools.audit import audit  # noqa: E402
from geometry_tools.common import atomic_json, atomic_npz, digest_file  # noqa: E402
from geometry_tools.model import Mesh  # noqa: E402


SCHEMA = 'native-side-tributary-finite-wall/v1'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def bound_file(binding, base):
    path = Path(binding['path'])
    path = path if path.is_absolute() else base / path
    path = path.resolve()
    require(path.is_file(), f'Missing declared source: {path}')
    require(digest_file(path) == binding['sha256'], f'Stale declared source: {path}')
    return path


def load_plan(args):
    path = args.plan.resolve()
    require(digest_file(path) == args.plan_sha256, 'Plan SHA256 mismatch')
    plan = json.loads(path.read_text())
    require(plan.get('schema') == SCHEMA, 'Unsupported plan schema')
    require(plan.get('plan_validated_for_construction') is True,
            'Root must explicitly validate the source-bound construction plan')
    require(plan.get('anatomical_acceptance') is False,
            'Construction inputs must not claim anatomical acceptance')
    require(isinstance(plan.get('frame'), str) and plan['frame'], 'Explicit frame required')
    bound_file(plan['inventory'], path.parent)
    require(bool(plan.get('guards')), 'At least one current catalogue/context guard required')
    for binding in plan['guards']:
        bound_file(binding, path.parent)
    source_paths = [bound_file(plan[k], path.parent) for k in ('parent', 'tributary')]
    require(source_paths[0] != source_paths[1], 'Two distinct indexed native sources required')
    require(plan['branch_mouth'] in (0, 1), 'branch_mouth must be 0 or 1')
    fraction = float(plan['inner_radial_fraction'])
    require(0 < fraction < 1, 'Positive, subunit illustrative inner radial fraction required')
    require(plan['mode'] in ('native-overlap', 'extend'), 'Explicit native-overlap or extend mode required')
    if plan['mode'] == 'extend':
        extension = plan['extension']
        require(0 < float(extension['max_length_m']) <= .005,
                'Declared extension maximum must be positive and at most 5 mm')
        for key in ('spacing_m', 'handle_m', 'max_station_turn_rad', 'curvature_support_limit'):
            require(np.isfinite(float(extension[key])) and float(extension[key]) > 0,
                    f'Positive finite extension {key} required')
        require(float(extension['spacing_m']) <= .001, 'Extension stations must be <= 1 mm apart')
        require(float(extension['max_station_turn_rad']) < np.pi / 2, 'Explicit turn guard too loose')
        require(float(extension['curvature_support_limit']) < 1, 'Curvature support must remain below 1')
    else:
        require(0 < float(plan['max_axis_residual_m']) <= 1e-7,
                'Native overlap needs an explicit <= 0.1-micron axis residual guard')
        require('extension' not in plan, 'Native-overlap mode must not declare an extension')
    require(0 < float(plan['junction_support_radius_m']) <= .01,
            'Explicit local Boolean support radius, at most 10 mm, required')
    origin = np.asarray(plan['storage_origin_m'], float)
    require(origin.shape == (3,) and np.isfinite(origin).all(), 'Finite 3-vector storage origin required')
    return path, plan, source_paths


def source_arrays(path, frame):
    with np.load(path, allow_pickle=False) as z:
        arrays = {k: z[k].copy() for k in ('vertices', 'faces', 'profiles', 'profile_vertex_indices')}
        if 'frame' in z:
            require(str(z['frame']) == frame, f'Source frame mismatch: {path}')
    v, f, p, indices = (arrays[k] for k in ('vertices', 'faces', 'profiles', 'profile_vertex_indices'))
    require(v.ndim == 2 and v.shape[1] == 3 and np.isfinite(v).all(), 'Invalid native vertices')
    require(f.ndim == 2 and f.shape[1] == 3 and np.issubdtype(f.dtype, np.integer), 'Indexed integer native triangles required')
    require(f.min() >= 0 and f.max() < len(v), 'Native face index out of range')
    require(p.ndim == 3 and p.shape[0] >= 2 and p.shape[1] >= 3 and p.shape[2] == 3,
            'Two ordered native profiles required')
    require(indices.shape == p.shape[:2] and np.issubdtype(indices.dtype, np.integer), 'Profile index shape/type mismatch')
    require(np.array_equal(np.sort(indices.ravel()), np.arange(len(v))),
            'Profile indices must cover every unique native vertex exactly once')
    require(np.array_equal(v[indices], p), 'Profile controls must be literal indexed source coordinates')
    station = np.empty(len(v), dtype=np.int64)
    station[indices] = np.arange(len(p))[:, None]
    require(np.all(np.ptp(station[f], axis=1) == 1), 'Native faces must join adjacent station rings only')
    check = continuity_audit(v, f, expected_mouths=2, expected_components=1, expected_euler=0)
    require(check['topological_continuity_pass'], 'Source is not a single two-mouth native side sheet')
    loops = boundary_loops(f)
    require({frozenset(q) for q in loops} == {frozenset(indices[0]), frozenset(indices[-1])},
            'Native mouths do not match the declared ordered profile endpoints')
    return arrays


def parent_axis_anchor(plan, parent):
    anchor = plan['parent_anchor']
    j, fraction = int(anchor['station']), float(anchor.get('fraction', 0.))
    require(0 <= j < len(parent['profiles']) - 1 and 0 <= fraction <= 1,
            'Parent anchor must identify an actual source axis interval')
    pc = parent['profiles'].mean(1)
    return (1 - fraction) * pc[j] + fraction * pc[j + 1]


def extension_profiles(plan, parent, branch):
    """Append full spatial contours through proper, rotation-minimizing frames."""
    mouth = 0 if plan['branch_mouth'] == 0 else -1
    adjacent = 1 if mouth == 0 else -2
    p = branch['profiles']
    start = p[mouth].mean(0)
    outward = unit(start - p[adjacent].mean(0))
    target = parent_axis_anchor(plan, parent)
    displacement = target - start
    direct = float(np.linalg.norm(displacement))
    e = plan['extension']
    require(1e-7 < direct <= float(e['max_length_m']), 'Anchor beyond declared extension bound, or coincident')
    final = unit(np.asarray(e.get('arrival_direction_m', displacement), float))
    require(final.shape == (3,) and np.dot(final, displacement) > 0,
            'Arrival direction must enter the actual parent axis')
    handle = float(e['handle_m'])
    require(handle <= direct, 'Bezier handle may not exceed endpoint distance')
    controls = np.array([start, start + handle * outward, target - handle * final, target])
    # Fixed parameter grid is a proposal, not an adaptive retry or radius reduction.
    count = int(e.get('stations', max(9, np.ceil(3 * direct / float(e['spacing_m'])) + 1)))
    require(3 <= count <= 257, 'Extension station count must be 3..257')
    t = np.linspace(0., 1., count)
    guide = ((1-t)[:, None]**3 * controls[0] + 3*((1-t)**2*t)[:, None]*controls[1]
             + 3*((1-t)*t*t)[:, None]*controls[2] + t[:, None]**3*controls[3])
    tangent = (3*(1-t)[:, None]**2*(controls[1]-controls[0])
               + 6*((1-t)*t)[:, None]*(controls[2]-controls[1])
               + 3*t[:, None]**2*(controls[3]-controls[2]))
    guide[[0, -1]] = [start, target]
    tangent = unit(tangent)
    length = float(arclength(guide)[-1])
    # A Bezier curve lies in its control polygon; polygon length upper-bounds
    # the continuous arclength. Do not accept only a short sampled polyline.
    upper_length = float(np.linalg.norm(np.diff(controls, axis=0), axis=1).sum())
    require(upper_length <= float(e['max_length_m']),
            'Continuous extension length upper bound exceeds declared maximum')
    require(np.linalg.norm(np.diff(guide, axis=0), axis=1).max() <= float(e['spacing_m']),
            'Prepared guide exceeds explicit station spacing')
    frames = rotation_minimizing_frames(guide, tangent, initial_reference=p[mouth, 0] - start)
    rotations = frames @ frames[0].T
    require(np.max(abs(rotations.transpose(0, 2, 1) @ rotations - np.eye(3))) < 1e-12,
            'Nonorthogonal extension station transform')
    require(np.linalg.det(rotations).min() > 1 - 1e-12, 'Improper extension station transform')
    extended = np.einsum('sij,rj->sri', rotations, p[mouth] - start) + guide[:, None]
    extended[0] = p[mouth]  # Literal source boundary, whose frame transform is identity.
    before = np.linalg.norm(p[mouth, :, None] - p[mouth, None, :], axis=2)
    after = np.linalg.norm(extended[:, :, None] - extended[:, None, :], axis=3)
    isometry_error = float(np.max(abs(after - before)))
    require(isometry_error < 2e-12, 'Extension profile shape changed')
    bends = bend_diagnostics(guide, extended[None], tangent,
                             float(e['curvature_support_limit']))
    require(not len(bends['failed_stations']), 'Extension too tight for literal source profile')
    require(float(bends['station_turn_rad'].max()) <= float(e['max_station_turn_rad']),
            'Extension exceeds explicit station turn bound')
    return extended, target, dict(method='cubic Bezier with double-reflection proper station transport',
        sampled_length_m=length, continuous_length_upper_bound_m=upper_length,
        maximum_profile_shape_error_m=isometry_error, radii_scaled=False,
        maximum_curvature_support=float(bends['curvature_support'].max()),
        maximum_station_turn_rad=float(bends['station_turn_rad'].max()),
        sampled_bend_check_not_global_injectivity_proof=True, stations=count)


def prepare(out, plan, paths, began):
    out.mkdir(exist_ok=False)
    arrays = [source_arrays(p, plan['frame']) for p in paths]
    if plan['mode'] == 'extend':
        extension, anchor, proof = extension_profiles(plan, *arrays)
    else:
        anchor = parent_axis_anchor(plan, arrays[0])
        mouth = 0 if plan['branch_mouth'] == 0 else -1
        residual = float(np.linalg.norm(arrays[1]['profiles'][mouth].mean(0) - anchor))
        require(residual <= float(plan['max_axis_residual_m']),
                'Native tributary mouth is not at the declared parent axis anchor')
        extension = arrays[1]['profiles'][mouth][None].copy()
        proof = dict(method='literal native sidewall overlap; no extension',
                     axis_residual_m=residual, radii_scaled=False, stations=0)
    # Preserve every original branch triangle; orient the appended loft from
    # the actual directed open boundary rather than guessing ring winding.
    branch = arrays[1]
    ring = branch['profile_vertex_indices'][0 if plan['branch_mouth'] == 0 else -1]
    loop = next(q for q in boundary_loops(branch['faces']) if set(q) == set(ring))
    indices = np.empty((len(extension), len(ring)), dtype=np.int64)
    indices[0] = ring
    indices[1:] = np.arange(len(branch['vertices']), len(branch['vertices']) + (len(extension)-1)*len(ring)).reshape(-1, len(ring))
    position = {int(v): j for j, v in enumerate(ring)}
    order = np.array([position[int(v)] for v in loop])
    strips = []
    for a, b in zip(indices[:-1], indices[1:]):
        for j in range(len(ring)):
            u, v = int(order[j]), int(order[(j+1) % len(ring)])
            # Oppose the existing directed boundary on the source-facing edge.
            strips.extend([[int(a[v]), int(a[u]), int(b[u])],
                           [int(a[v]), int(b[u]), int(b[v])]])
    branch_vertices = np.r_[branch['vertices'], extension[1:].reshape(-1, 3)]
    branch_faces = np.r_[branch['faces'], np.asarray(strips, dtype=np.int64).reshape(-1, 3)]
    near = 0 if plan['branch_mouth'] == 0 else -1
    far = -1 if near == 0 else 0
    alpha = float(plan['inner_radial_fraction'])
    inner = []
    for a in arrays:
        centers = a['profiles'].mean(1)
        q = np.empty_like(a['vertices'])
        q[a['profile_vertex_indices']] = centers[:, None] + alpha*(a['profiles']-centers[:, None])
        inner.append(q)
    ex_center = extension.mean(1)
    inner[1] = np.r_[inner[1], (ex_center[:, None] + alpha*(extension-ex_center[:, None]))[1:].reshape(-1, 3)]
    minimum_thickness = min(float(np.linalg.norm(a-b, axis=1).min())
                            for a, b in zip([arrays[0]['vertices'], branch_vertices], inner))
    require(minimum_thickness > 1e-8, 'A source profile collapses its illustrative finite wall')
    mouths = dict(parent=[arrays[0]['profile_vertex_indices'][0].tolist(), arrays[0]['profile_vertex_indices'][-1].tolist()],
                  tributary=[branch['profile_vertex_indices'][far].tolist(), indices[-1].tolist()])
    check = continuity_audit(branch_vertices, branch_faces, expected_mouths=2, expected_components=1, expected_euler=0)
    require(check['topological_continuity_pass'], 'Appended branch sheet fails literal continuity')
    require(np.array_equal(branch_vertices[:len(branch['vertices'])], branch['vertices']), 'Native branch coordinates moved')
    require(np.array_equal(branch_faces[:len(branch['faces'])], branch['faces']), 'Native branch faces changed before Boolean')
    atomic_npz(out/'prepared.npz', parent_outer=arrays[0]['vertices'], parent_inner=inner[0],
        parent_faces=arrays[0]['faces'], tributary_outer=branch_vertices, tributary_inner=inner[1],
        tributary_faces=branch_faces, anchor=anchor, extension_profiles=extension,
        parent_native_face_count=np.array(len(arrays[0]['faces'])),
        tributary_native_face_count=np.array(len(branch['faces'])))
    return dict(prepared_sha256=digest_file(out/'prepared.npz'), mouths=mouths,
        extension=proof, native_coordinates_and_faces_literal_before_boolean=True,
        minimum_corresponding_wall_thickness_m=minimum_thickness,
        illustrative_radial_wall_fraction=alpha, not_a_TA98_wall_coat=True,
        anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)


def stage_inputs(out):
    prep = json.loads((out/'prepare.json').read_text())
    require(digest_file(out/'prepared.npz') == prep['prepared_sha256'], 'Prepared arrays changed')
    with np.load(out/'prepared.npz', allow_pickle=False) as z:
        arrays = {k: z[k].copy() for k in z.files}
    return prep, arrays


def sheet(out, plan, level, began):
    prep, arrays = stage_inputs(out)
    origin = np.asarray(plan['storage_origin_m'], float)
    inputs, solids, source_points = {}, [], []
    for component, name in enumerate(('parent', 'tributary')):
        vertices = arrays[f'{name}_{level}'].copy()
        source_points.append(vertices.copy())
        faces = arrays[f'{name}_faces'].copy()
        native_count = int(arrays[f'{name}_native_face_count'])
        external_loops = prep['mouths'][name][:2 if name == 'parent' else 1]
        caps = {}
        for loop in boundary_loops(faces):
            external = any(set(loop) == set(q) for q in external_loops)
            tip = len(vertices)
            vertices = np.r_[vertices, vertices[loop].mean(0)[None]]
            start = len(faces)
            faces = np.r_[faces, np.array([[tip, int(loop[(j+1) % len(loop)]), int(loop[j])]
                                           for j in range(len(loop))])]
            caps.update({k: external for k in range(start, len(faces))})
        solid, oid = make_solid((vertices-origin)*1000, faces)
        require(solid.status() == md.Error.NoError, f'Invalid {name} Boolean classifier')
        inputs[oid] = dict(component=component, native_faces=native_count, caps=caps)
        solids.append(solid)
    joined = md.Manifold.batch_boolean(solids, md.OpType.Add)
    require(joined.status() == md.Error.NoError, f'{level} Boolean union failed')
    raw = joined.to_mesh64()
    positions = np.asarray(raw.vert_properties)[:, :3]/1000 + origin
    faces = np.asarray(raw.tri_verts, dtype=np.int64)
    owner, source_face = provenance(raw)
    external = np.array([inputs[int(o)]['caps'].get(int(f)) is True for o, f in zip(owner, source_face)])
    internal = np.array([inputs[int(o)]['caps'].get(int(f)) is False for o, f in zip(owner, source_face)])
    require(not internal.any(), 'Internal classifier disk survives: junction is not patent')
    classifier_components = np.array([inputs[int(o)]['component'] for o in owner], dtype=np.int64)
    classifier_source_faces = source_face.copy()
    bank = CommonPositionBank(np.concatenate(source_points), tolerance=2e-12)
    mapping = np.array([bank.insert(p) for p in positions], dtype=np.int64)
    candidate = Mesh(bank.positions, mapping[faces[~external]], plan['frame'])
    continuity = continuity_audit(candidate.vertices, candidate.faces,
                                 expected_mouths=3, expected_components=1, expected_euler=-1)
    require(continuity['topological_continuity_pass'], 'Union is not one three-mouth sheet')
    components = np.array([inputs[int(o)]['component'] for o in owner[~external]], dtype=np.int64)
    source_face = source_face[~external]
    # A cut outside the declared support is forbidden. Source triangles cut by
    # the Boolean may be subdivided locally; all remote native triangles remain
    # literal, with exact oriented coordinates from CommonPositionBank.
    emitted = {oriented_key(face) for face in candidate.faces}
    preservation = []
    for component, name in enumerate(('parent', 'tributary')):
        vertices = arrays[f'{name}_{level}']
        native_faces = arrays[f'{name}_faces'][:int(arrays[f'{name}_native_face_count'])]
        triangles = vertices[native_faces]
        literal_indices = np.array([bank.insert(p) for p in vertices], dtype=np.int64)
        retained = np.array([oriented_key(face) in emitted for face in literal_indices[native_faces]])
        changed = np.flatnonzero(~retained)
        radius = float(plan['junction_support_radius_m'])
        if len(changed):
            extent = np.linalg.norm(triangles[changed] - arrays['anchor'], axis=2).max()
            require(extent <= radius + 2e-12,
                    f'Boolean changed a native {name} triangle outside declared support')
        preservation.append(dict(source=name, native_faces=len(native_faces),
            literal_retained_faces=int(retained.sum()), locally_cut_native_face_ids=changed.tolist(),
            changed_native_triangles_within_declared_support=True))
    candidate.save(out/f'{level}-sheet.npz')
    atomic_npz(out/f'{level}-ancestry.npz', component=components, source_face=source_face)
    classifier = Mesh(bank.positions, mapping[faces], plan['frame'])
    classifier.save(out/f'{level}-classifier.npz')
    # These rows correspond one-for-one to the actual classifier's indexed
    # triangles. A subsequent Boolean difference must retain this distinction
    # between virtual closure fans and actual native/generated tube walls.
    atomic_npz(out/f'{level}-classifier-provenance.npz',
        external_virtual_cap=external, internal_virtual_cap=internal,
        component=classifier_components, source_face=classifier_source_faces)
    return dict(sheet_sha256=digest_file(out/f'{level}-sheet.npz'),
        ancestry_sha256=digest_file(out/f'{level}-ancestry.npz'),
        classifier_sha256=digest_file(out/f'{level}-classifier.npz'),
        classifier_provenance_sha256=digest_file(out/f'{level}-classifier-provenance.npz'),
        internal_classifier_caps=0, external_classifier_faces_removed=int(external.sum()),
        external_blood_mouths=3, continuity=continuity, source_preservation=preservation,
        anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)


def assemble(out, plan, began):
    prep, arrays = stage_inputs(out)
    sheets = []
    ancestry = []
    classifier_provenance = []
    for level in ('outer', 'inner'):
        record = json.loads((out/f'{level}.json').read_text())
        for field, suffix in [('sheet_sha256', 'sheet.npz'), ('ancestry_sha256', 'ancestry.npz'),
                              ('classifier_sha256', 'classifier.npz'),
                              ('classifier_provenance_sha256', 'classifier-provenance.npz')]:
            require(digest_file(out/f'{level}-{suffix}') == record[field], f'{level} stage changed')
        sheets.append(Mesh.load(out/f'{level}-sheet.npz'))
        with np.load(out/f'{level}-ancestry.npz', allow_pickle=False) as z:
            ancestry.append({k: z[k].copy() for k in z.files})
        with np.load(out/f'{level}-classifier-provenance.npz', allow_pickle=False) as z:
            classifier_provenance.append({k: z[k].copy() for k in z.files})
    outer, inner = sheets
    origin = np.asarray(plan['storage_origin_m'], float)
    classifiers = [Mesh.load(out/f'{level}-classifier.npz') for level in ('outer', 'inner')]
    solid_rows = [make_solid((m.vertices-origin)*1000, m.faces) for m in classifiers]
    solids = [row[0] for row in solid_rows]
    classifier_by_original_id = {int(row[1]): i for i, row in enumerate(solid_rows)}
    require(len(classifier_by_original_id) == 2, 'Classifier source IDs must remain distinct')
    for m, p in zip(classifiers, classifier_provenance):
        require(all(len(p[k]) == len(m.faces) for k in
                    ('external_virtual_cap', 'internal_virtual_cap', 'component', 'source_face')),
                'Classifier ancestry does not cover every actual indexed triangle')
        require(p['external_virtual_cap'].dtype == np.bool_ and p['internal_virtual_cap'].dtype == np.bool_,
                'Classifier cap provenance must use literal Boolean masks')
        require(not p['internal_virtual_cap'].any(), 'An internal virtual cap entered a saved classifier')
    outside = solids[1] - solids[0]
    require(outside.status() == md.Error.NoError, 'Inner-minus-outer Boolean classification failed')
    strict_empty = outside.num_tri() == 0
    outside_raw = outside.to_mesh64()
    outside_positions = np.asarray(outside_raw.vert_properties)[:, :3].reshape(-1, 3)/1000 + origin
    outside_faces = np.asarray(outside_raw.tri_verts, dtype=np.int64).reshape(-1, 3)
    outside_owner, outside_face_id = provenance(outside_raw)
    outside_level = np.full(len(outside_faces), -1, dtype=np.int64)
    outside_external_cap = np.zeros(len(outside_faces), dtype=bool)
    outside_component = np.full(len(outside_faces), -1, dtype=np.int64)
    outside_native_face = np.full(len(outside_faces), -1, dtype=np.int64)
    for j, (oid, face_id) in enumerate(zip(outside_owner, outside_face_id)):
        require(int(oid) in classifier_by_original_id, 'Outside facet has unknown Boolean ancestry')
        level = classifier_by_original_id[int(oid)]
        rows = classifier_provenance[level]
        require(0 <= int(face_id) < len(rows['external_virtual_cap']),
                'Outside facet refers to a nonexistent classifier triangle')
        outside_level[j] = level
        outside_external_cap[j] = bool(rows['external_virtual_cap'][int(face_id)])
        outside_component[j] = int(rows['component'][int(face_id)])
        outside_native_face[j] = int(rows['source_face'][int(face_id)])
    used = np.concatenate([m.vertices[np.unique(m.faces)] for m in classifiers])
    max_extent_mm = float(np.ptp(used, axis=0).max()*1000)
    roundoff_budget_m3 = float(64*np.finfo(np.float64).eps*max_extent_mm**3/1e9)
    outside_volume_m3 = float(outside.volume()/1e9)
    cap_only = bool(outside_external_cap.all())
    volume_bounded = bool(np.isfinite(outside_volume_m3) and
                          0 <= outside_volume_m3 <= roundoff_budget_m3)
    artifact = out/'inner-minus-outer-diagnostic.npz'
    atomic_npz(artifact, vertices=outside_positions, faces=outside_faces,
        classifier_original_id=outside_owner, classifier_face_id=outside_face_id,
        classifier_level=outside_level, external_virtual_cap=outside_external_cap,
        original_component=outside_component, original_input_face=outside_native_face,
        frame=np.array(plan['frame']))
    containment = dict(schema='finite-wall-classifier-containment/v1',
        strict_empty=strict_empty, exact_full_containment=strict_empty,
        contained_away_from_external_classifier_caps=cap_only,
        nonempty_difference_is_external_cap_roundoff_only=bool(not strict_empty and cap_only and volume_bounded),
        outside_faces=len(outside_faces), outside_vertices=len(outside_positions),
        outside_actual_wall_source_faces=int((~outside_external_cap).sum()),
        outside_volume_m3=outside_volume_m3, maximum_classifier_extent_mm=max_extent_mm,
        roundoff_budget_m3=roundoff_budget_m3,
        roundoff_budget_formula='64 * Float64 epsilon * maximum_classifier_extent_mm^3 / 1e9',
        finite_roundoff_volume_bound_passed=volume_bounded,
        outside_artifact=str(artifact), outside_artifact_sha256=digest_file(artifact),
        classifier_input_original_ids={str(oid): ('outer' if level == 0 else 'inner')
                                       for oid, level in classifier_by_original_id.items()},
        source_classifier_sha256={level: digest_file(out/f'{level}-classifier.npz')
                                  for level in ('outer', 'inner')},
        source_classifier_provenance_sha256={level: digest_file(out/f'{level}-classifier-provenance.npz')
                                             for level in ('outer', 'inner')},
        virtual_external_caps_are_not_rendered_tissue=True, anatomical_acceptance=False,
        passed=bool(strict_empty or cap_only and volume_bounded))
    atomic_json(out/'inner-containment.json', containment)
    require(containment['passed'],
            'Inner classifier leaks through an actual wall, or virtual-cap difference exceeds roundoff budget')
    tissue_volume = float((solids[0].volume() - solids[1].volume()) / 1e9)
    require(tissue_volume > 0, 'Nonpositive finite tissue volume')
    offset = len(outer.vertices)
    vertices = np.r_[outer.vertices, inner.vertices]
    pieces = [outer.faces, inner.faces[:, ::-1] + offset]
    ports = []
    ol, il = boundary_loops(outer.faces), boundary_loops(inner.faces)
    used_o, used_i = set(), set()
    for name, mouth in [('parent', 0), ('parent', 1), ('tributary', 0)]:
        indices = np.array(prep['mouths'][name][mouth])
        expected_o = arrays[f'{name}_outer'][indices]
        expected_i = arrays[f'{name}_inner'][indices]
        center = expected_o.mean(0)
        oi = min(range(len(ol)), key=lambda j: np.linalg.norm(outer.vertices[ol[j]].mean(0)-center))
        ii = min(range(len(il)), key=lambda j: np.linalg.norm(inner.vertices[il[j]].mean(0)-center))
        require(oi not in used_o and ii not in used_i, 'Two external mouth identities select the same boundary')
        used_o.add(oi)
        used_i.add(ii)
        a, b = ol[oi], il[ii]
        require(len(a) == len(b) == len(indices), 'Boolean altered an external source profile')
        ao = np.argmin(np.linalg.norm(expected_o[None] - outer.vertices[a][:, None], axis=2), axis=1)
        require(len(set(ao)) == len(indices), 'External outer mouth correspondence not bijective')
        require(np.linalg.norm(outer.vertices[a] - expected_o[ao], axis=1).max() < 2e-12,
                'External outer profile changed')
        bi = np.argmin(np.linalg.norm(inner.vertices[b][None] - expected_i[ao][:, None], axis=2), axis=1)
        require(len(set(bi)) == len(indices), 'External inner mouth correspondence not bijective')
        require(np.linalg.norm(inner.vertices[b[bi]] - expected_i[ao], axis=1).max() < 2e-12,
                'External inner profile changed')
        rim = annular_rim(a, b[bi] + offset)
        pieces.append(rim)
        ports.append(dict(source=name, mouth=mouth, outer_indices=a.tolist(),
                          inner_indices=(b[bi]+offset).tolist()))
    wall = Mesh(vertices, np.concatenate(pieces), plan['frame'])
    wall.save(out/'finite-wall.npz')
    wall_check = audit(wall, require_closed=True, expected_components=1, expected_genus=[2],
                       max_pairs=250000, max_tests=2000000)
    atomic_json(out/'wall-audit.json', wall_check)
    require(wall_check['passed'], f'Finite wall audit failed: {wall_check["failures"]}')
    blood = Mesh.load(out/'inner-classifier.npz')
    blood_check = audit(blood, require_closed=True, expected_components=1, expected_genus=[0],
                        max_pairs=250000, max_tests=2000000)
    atomic_json(out/'blood-audit.json', blood_check)
    require(blood_check['passed'], f'Blood classifier audit failed: {blood_check["failures"]}')
    atomic_npz(out/'surfaces-and-ancestry.npz', vertices=vertices, outer_faces=outer.faces,
        inner_faces=inner.faces[:, ::-1]+offset, outer_component=ancestry[0]['component'],
        outer_source_face=ancestry[0]['source_face'], inner_component=ancestry[1]['component'],
        inner_source_face=ancestry[1]['source_face'],
        **{f'annulus_{i}': faces for i, faces in enumerate(pieces[2:])})
    atomic_json(out/'ports.json', dict(ports=ports))
    return dict(passed=True, finite_wall_sha256=digest_file(out/'finite-wall.npz'),
        geometry_sha256=wall.geometry_sha, wall_faces=len(wall.faces), external_blood_mouths=3,
        closed_finite_wall_genus=2, closed_blood_classifier_genus=0, internal_classifier_caps=0,
        inner_blood_classifier_strictly_contained=strict_empty,
        contained_away_from_external_classifier_caps=cap_only,
        inner_containment_sha256=digest_file(out/'inner-containment.json'), tissue_volume_m3=tissue_volume,
        lumen_disks=0, finite_wall_not_anatomical_layer_claim=True,
        source_attribute_transport_and_packing_pending=True, organ_context_and_render_review_pending=True,
        anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan', type=Path, required=True)
    parser.add_argument('--plan-sha256', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', choices=('prepare', 'outer', 'inner', 'assemble'), required=True)
    args = parser.parse_args()
    began = time.monotonic()
    plan_path, plan, paths = load_plan(args)
    out = args.output.resolve()
    code_sha = digest_file(Path(__file__).resolve())
    if args.stage == 'prepare':
        require(not out.exists(), 'Prepare output must be a new directory')
    else:
        receipt = json.loads((out/'prepare.json').read_text())
        require(receipt['plan_sha256'] == args.plan_sha256 and receipt['producer_sha256'] == code_sha,
                'Prepared plan/producer changed; make a distinct new candidate')
        require(not (out/f'{args.stage}.json').exists(), 'Stage receipt already exists; no silent retries')
    try:
        if args.stage == 'prepare':
            result = prepare(out, plan, paths, began)
        elif args.stage in ('outer', 'inner'):
            result = sheet(out, plan, args.stage, began)
        else:
            result = assemble(out, plan, began)
        result.update(schema=SCHEMA, stage=args.stage, plan=str(plan_path),
                      plan_sha256=args.plan_sha256, producer_sha256=code_sha)
        atomic_json(out/f'{args.stage}.json', result)
        print(json.dumps(result, default=str))
    except Exception as error:
        if out.exists() and not (out/f'{args.stage}-failure.json').exists():
            atomic_json(out/f'{args.stage}-failure.json', dict(stage=args.stage,
                error_type=type(error).__name__, error=str(error), plan_sha256=args.plan_sha256,
                producer_sha256=code_sha, installed=False, anatomical_acceptance=False,
                elapsed_seconds=time.monotonic()-began))
        raise


if __name__ == '__main__':
    main()
