"""Source-bound bladder apical cord candidate; never installs an atlas.

The terminal fascia is an existing open reference sheet. Its contact is not
silently promoted to a reciprocal interface between two finite tissue cells.
"""
from pathlib import Path
import json
import numpy as np
from scipy.interpolate import CubicSpline

from .algorithms.annular_pipe_envelopes import boundary_loops
from .algorithms.solid_attachment import bridge_patches
from .algorithms.surface_normals import grouped_normal_words
from .audit import audit
from .bladder_features import back_graph
from .common import atomic_json, atomic_npz, digest_file
from .contact import pair_contacts
from .model import AtlasReader, Mesh
from .tissue import pack


def patch_centroid(mesh):
    t = mesh.triangles
    w = np.linalg.norm(np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), axis=1)
    return np.einsum('i,ij->j', w, t.mean(1))/w.sum()


def construct(plan_path, output):
    plan_path, output = Path(plan_path).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    plan = json.loads(plan_path.read_text()); source = Path(plan['source_dir'])
    manifest_path = source/'sources-and-assembly.json'
    if digest_file(manifest_path) != plan['source_manifest_sha256']:
        raise ValueError('Apical source manifest changed')
    source_info = json.loads(manifest_path.read_text())
    if digest_file(source_info['bladder_exterior']) != source_info['bladder_exterior_sha256']:
        raise ValueError('The fixed bladder exterior changed')
    cat = Path(source_info['source_catalogue'])
    if digest_file(cat) != source_info['source_catalogue_sha256']:
        raise ValueError('Reprepare from the actual current bladder revision')
    for r in source_info['context']:
        if digest_file(r['path']) != r['sha256']:
            raise ValueError('A fixed context mesh changed')
    a = json.loads(cat.read_text()); byid = {p['id']: p for p in a['parts']}
    record = byid[plan['first_part_id']]; reader = AtlasReader()
    host = reader.read_part(a, record, frame=source_info['frame'], catalogue_path=cat)
    ids = np.asarray(plan['first_faces'], int)
    if len(set(ids)) != len(ids) or np.any(ids < 0) or np.any(ids >= len(host.faces)):
        raise ValueError('Explicit unique source contact faces are required')
    if not all(any(not g.get('neighbor') and g['start']//3 <= i < (g['start']+g['count'])//3
                   for g in record['tissueDrawGroups']) for i in ids):
        raise ValueError('The apex contact must remain on actual exposed outer tissue')
    first = Mesh(host.vertices, host.faces[ids], host.frame, host.attributes)
    fascia_row = next(r for r in source_info['context'] if r['part_id'] == plan['second_part_id'])
    fascia = Mesh.load(fascia_row['path']); face = int(plan['second_face'])
    second = Mesh(fascia.vertices, fascia.faces[[face]], fascia.frame, fascia.attributes)
    c0, c1 = patch_centroid(first), patch_centroid(second)
    normal = lambda m: np.cross(m.triangles[:, 1]-m.triangles[:, 0],
                               m.triangles[:, 2]-m.triangles[:, 0]).sum(0)
    n0, n1 = normal(first), normal(second)
    n0 /= np.linalg.norm(n0); n1 /= np.linalg.norm(n1)
    if c1[1]-c0[1] < .03:
        raise ValueError('The explicitly selected umbilical attachment must be superior to the apex')
    h0, h1 = plan['first_transition_length_m'], plan['final_transition_length_m']
    start_y, end_y = c0[1]+h0, c1[1]-h1
    waypoint_step = plan.get('waypoint_step_m')
    if waypoint_step is not None and (not np.isfinite(waypoint_step) or waypoint_step <= 0):
        raise ValueError('Fascial sampling step must be finite and positive')
    waypoint_count = max(10, int(np.ceil((end_y-start_y)/waypoint_step))+1) if waypoint_step else 10
    ys = np.linspace(start_y, end_y, waypoint_count)
    xy = np.c_[np.zeros(len(ys)), ys]
    z, parents = back_graph(fascia.vertices*[1, 1, -1], fascia.faces, xy)
    if not np.isfinite(z).all():
        raise ValueError('The proposed course leaves the actual anterior fascial sheet')
    radii = np.asarray(plan['body_radii_m'], float)
    if radii.shape != (2,) or np.any(radii <= 0):
        raise ValueError('Two explicit physical section radii are required')
    waypoints = np.vstack([c0, np.c_[xy, -z-radii.max()-plan['anterior_wall_gap_m']], c1])
    arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(waypoints, axis=0), axis=1))]
    spline = CubicSpline(arc, waypoints, bc_type=((1, n0), (1, n1)))
    count = int(np.ceil(arc[-1]/plan['maximum_step_m']))+1
    t = np.linspace(0, arc[-1], count); centers = spline(t)
    # Gentle taper. End contacts themselves retain their actual source shape.
    distance = np.minimum(t, arc[-1]-t)
    taper = np.clip(distance/.009, 0, 1); taper = taper*taper*(3-2*taper)
    section_radii = radii[None, :]*(1.15-.15*taper[:, None])
    positions, surfaces, evidence = bridge_patches(first, second, centers, section_radii,
        first_side=-1, second_side=1, section_controls=plan['sides'])
    words, normal_evidence = grouped_normal_words(positions, {'free_wall': surfaces['free_wall']},
        {'free_wall': np.full(len(surfaces['free_wall']), 'fibromuscular-cord', object)})
    def source_words(mesh):
        values = mesh.attributes['normals'].astype('int32')
        if 'normalSigns' in mesh.attributes:
            values *= mesh.attributes['normalSigns'][:, None]
        return values[mesh.faces]
    words['first_contact'] = (-source_words(first)[:, [0, 2, 1]]).astype('<i2')
    words['second_contact'] = source_words(second).astype('<i2')
    output.mkdir(parents=True)
    keys = {name: 'faces_'+str(i) for i, name in enumerate(surfaces)}
    normal_keys = {name: 'words_'+str(i) for i, name in enumerate(surfaces)}
    atomic_npz(output/'construction.npz', positions=positions, centers=centers, radii=section_radii,
               **{keys[s]: f for s, f in surfaces.items()}, **{normal_keys[s]: w for s, w in words.items()})
    pid = 'BLADDER:male:median-umbilical-ligament'
    uses = [dict(surface=s, side=1, neighbor=None) for s in surfaces]
    packed = pack(positions, surfaces, {pid: uses}, {pid: pid}, host.frame,
        output/'packed', origin=record['geometryOrigin'], corner_words=words,
        max_pairs=1000000, construction=evidence)
    actual = Mesh.load(output/'packed'/(pid.replace(':', '_')+'.npz'))
    offsets = {}; start = 0
    for name, f in surfaces.items():
        offsets[name] = np.arange(start, start+len(f)); start += len(f)
    first_literal = np.array_equal(actual.triangles[offsets['first_contact']], first.triangles[:, [0, 2, 1]])
    second_literal = np.array_equal(actual.triangles[offsets['second_contact']], second.triangles)
    own_check = audit(actual, require_closed=True, expected_components=1, expected_genus=[0], max_pairs=1000000)
    exterior = Mesh.load(source_info['bladder_exterior'])
    crossing_arrays, crossing = pair_contacts(actual, exterior, 0, max_pairs=1000000)
    fascia_arrays, fascia_contact = pair_contacts(actual, fascia, 0, max_pairs=1000000)
    source_parents = json.loads(Path(source_info['face_parents']).read_text())
    overlap_pairs = crossing_arrays['pairs'][crossing_arrays['relation'] == 'coplanar-overlap']
    bladder_patch_only = (len(overlap_pairs) == len(ids)
        and set(overlap_pairs[:, 0]) == set(offsets['first_contact'])
        and {(source_parents[i]['part_id'], source_parents[i]['face_id']) for i in overlap_pairs[:, 1]}
            == {(plan['first_part_id'], int(i)) for i in ids})
    fascia_overlap = fascia_arrays['pairs'][fascia_arrays['relation'] == 'coplanar-overlap']
    fascia_patch_only = (len(fascia_overlap) == 1
        and np.array_equal(fascia_overlap[0], [offsets['second_contact'][0], face]))
    lining_vertices, lining_faces, count = [], [], 0
    for p in a['parts']:
        if p.get('suppressed') or p.get('urinaryOrgan') != 'bladder' \
                or not p.get('tissueWallLayer', '').startswith('mucosa'):
            continue
        m = reader.read_part(a, p, frame=host.frame, catalogue_path=cat)
        free = np.concatenate([np.arange(g['start']//3, (g['start']+g['count'])//3)
                               for g in p['tissueDrawGroups'] if not g.get('neighbor')])
        lining_vertices.append(m.vertices); lining_faces.append(m.faces[free]+count); count += len(m.vertices)
    lining = Mesh(np.concatenate(lining_vertices), np.concatenate(lining_faces), host.frame)
    _, lining_contact = pair_contacts(actual, lining, 0, max_pairs=1000000)
    passed = bool(packed['passed'] and own_check['passed'] and first_literal and second_literal
                  and not crossing['truncated'] and not crossing['crossing_face_pairs']
                  and bladder_patch_only and fascia_patch_only
                  and not fascia_contact['truncated'] and not fascia_contact['crossing_face_pairs']
                  and not lining_contact['truncated'] and not lining_contact['returned_triangle_pairs'])
    atomic_npz(output/'bladder-contacts.npz', **crossing_arrays)
    atomic_npz(output/'fascia-contacts.npz', **fascia_arrays)
    report = dict(schema='bladder-apical-cord-candidate/v1', plan=str(plan_path), plan_sha256=digest_file(plan_path),
        catalogue_sha256=digest_file(cat), packed=packed, own_check=own_check,
        first_source_contact_literal=first_literal, second_source_contact_literal=second_literal,
        bladder_contact=crossing, fascia_contact=fascia_contact,
        bladder_contact_only_declared_source_patch=bladder_patch_only,
        fascia_contact_only_declared_source_patch=fascia_patch_only,
        actual_bladder_lining_contact=lining_contact,
        normal_evidence=normal_evidence, source_boundaries=evidence,
        first_contact_source=dict(part_id=plan['first_part_id'], face_ids=ids.tolist()),
        second_contact_source=dict(part_id=plan['second_part_id'], face_ids=[face],
                                   kind='open-reference-sheet-contact; not a finite reciprocal tissue interface'),
        apex_point=c0.tolist(), umbilical_fascial_point=c1.tolist(),
        approximate_course_length_m=float(np.linalg.norm(np.diff(centers, axis=0), axis=1).sum()),
        source_organ_and_skeleton_geometry_modified=False, passed=passed, installed=False,
        anatomical_acceptance=False, limits=['Illustrative cord thickness; not donor histometry.',
            'Deep umbilical ring has no independent source mesh; contact uses native fascia beneath its skin landmark.',
            'Full-body packing and actual surface/section review remain required.'])
    atomic_json(output/'candidate.json', report)
    return report


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('plan'); parser.add_argument('--output', required=True)
    args = parser.parse_args(); result = construct(args.plan, args.output)
    print(json.dumps({k: result[k] for k in ['passed', 'first_source_contact_literal',
        'second_source_contact_literal', 'approximate_course_length_m', 'bladder_contact', 'fascia_contact']}, indent=2))
    raise SystemExit(0 if result['passed'] else 1)
