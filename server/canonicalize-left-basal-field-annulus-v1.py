"""Replay the declared canonical basal interface from fresh Boolean ancestry.

Exactly one ancestral proximal-cap/first-inner-strip patch is replaced. No
course fitting, radius change, Boolean retry or live atlas mutation occurs.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.tube_continuity import boundary_loops
from geometry_tools.algorithms.finite_tissue_boolean import boundary_signature


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=lambda x: x.tolist() if hasattr(x, 'tolist') else x.item()) + '\n')


def run(directory, resume_attribute_order_only=False):
    began = time.monotonic()
    contract = json.loads((directory / 'contract-before-construction.json').read_text())
    for key in ['fields', 'frozen_inputs', 'fixed_parent', 'field_proof']:
        assert sha(contract[key]['path']) == contract[key]['sha256']
    assert sha(ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json') == contract['source_catalogue_sha256']
    old = directory / 'candidate-stage-v1'
    proof = json.loads((old / 'proof.json').read_text())
    assert proof['complete'] and proof['passed']
    assert sha(old / 'finite-basal-wall.npz') == proof['wall_npz_sha256']
    assert sha(old / 'source-primitives.npz') == proof['source_primitives_sha256']
    out = directory / 'candidate-stage-v2'
    if resume_attribute_order_only:
        assert out.exists() and not (out / 'proof.json').exists()
        assert (out / 'finite-basal-wall.npz').exists() and (out / 'patch-and-ancestry.npz').exists()
    else:
        assert not out.exists(), 'One canonical patch replay only.'
    x = np.load(directory / 'construction-inputs.npz', allow_pickle=False)
    assert sha(directory / 'construction-inputs.npz') == contract['construction_inputs_sha256']
    m = Mesh.load(old / 'finite-basal-wall.npz')
    a = np.load(old / 'finite-basal-wall-ancestry.npz', allow_pickle=False)
    patch = np.zeros(len(m.faces), bool)
    ancestry_selection = []
    for primitive in proof['primitives']:
        if primitive['curve'] != 131:
            continue
        port = next(q for q in primitive['mouths'] if q['mouth'] == 0)
        ids = np.asarray(port['newFaceIds'], int)
        if primitive['role'] == 'fluid':
            # close_branch emits two numerical guard-side facets immediately
            # before each classifier disk. Native first strip is24faces.
            ids = np.unique(np.r_[ids, ids-1, ids-2, np.arange(24)])
        take = (a['primitive_ids'] == primitive['oid']) & np.isin(a['source_face_ids'], ids)
        patch |= take
        ancestry_selection.append(dict(curve=131, role=primitive['role'], primitive=primitive['oid'], source_faces=ids.tolist(), output_faces=np.flatnonzero(take).tolist()))
    assert patch.any()
    physical, inverse = np.unique(m.vertices, axis=0, return_inverse=True)
    loops = boundary_loops(inverse[m.faces[patch]])
    assert sorted(map(len, loops)) == [12, 12], 'Ancestral patch must have exactly two12control boundaries; no branch mouth may be consumed.'
    # Preserve every existing row and all outside oriented faces. Only missing
    # canonical inner controls may be appended; old controls never move.
    vertices = m.vertices.tolist()
    exact = {tuple(p): i for i, p in enumerate(m.vertices)}
    new_ids = []
    nearest_new = []
    tree = cKDTree(m.vertices)
    def literal_id(point):
        key = tuple(point)
        if key not in exact:
            error, _ = tree.query(point)
            assert error <= 4.348e-9, 'Canonical correction left its already declared bounded root neighborhood.'
            nearest_new.append(float(error))
            exact[key] = len(vertices)
            new_ids.append(len(vertices))
            vertices.append(np.asarray(point).tolist())
        return exact[key]
    cap_corners = x['parent_annulus_oriented_corners'][:, [0, 2, 1]]
    cap_faces = np.asarray([[literal_id(p) for p in tri] for tri in cap_corners], int)
    inner = x['curve_131_inner'].reshape(-1, 12, 3)
    root_ids = np.asarray([literal_id(p) for p in inner[0]], int)
    profile1_error, profile1_ids = tree.query(inner[1])
    assert float(profile1_error.max()) <= 1e-12
    # This keeps actual pre-existing boundary controls literal, including the
    # measured sub-picometre Boolean arithmetic at nonshared native terminals.
    strip = x['curve_131_faces']
    strip = strip[np.all(strip < 24, axis=1)]
    assert len(strip) == 24
    strip_ids = np.r_[root_ids, profile1_ids]
    strip_faces = strip_ids[strip][:, [0, 2, 1]]
    positions = np.asarray(vertices, float)
    combined = np.r_[cap_faces, strip_faces]
    # Boundary comparison uses exact physical coordinates; no tolerance weld,
    # fabricated fan cap or erased branch opening is permitted here.
    old_triangles = m.vertices[m.faces[patch]]
    new_triangles = positions[combined]
    both, rev = np.unique(np.r_[old_triangles.reshape(-1, 3), new_triangles.reshape(-1, 3)], axis=0, return_inverse=True)
    old_f = rev[:old_triangles.size//3].reshape(-1, 3)
    new_f = rev[old_triangles.size//3:].reshape(-1, 3)
    assert boundary_signature(old_f) == boundary_signature(new_f), 'The patch has a changed physical boundary; stop.'
    remaining = m.faces[~patch]
    faces = np.r_[remaining, combined]
    mesh = Mesh(positions, faces, m.frame)
    n = x['parent_distal_normal']
    cap_area = np.cross(cap_corners[:, 1]-cap_corners[:, 0], cap_corners[:, 2]-cap_corners[:, 0])
    assert np.all(cap_area @ n < 0)
    t = positions[strip_faces]
    centres = inner[:2].mean(1)
    radial = positions[strip_ids] - np.repeat(centres, 12, axis=0)
    dots = np.einsum('ij,ij->i', np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), -radial[strip[:, [0, 2, 1]]].mean(1))
    assert np.all(dots > 0)
    if resume_attribute_order_only:
        saved = Mesh.load(out / 'finite-basal-wall.npz')
        assert np.array_equal(saved.vertices, mesh.vertices) and np.array_equal(saved.faces, mesh.faces)
    else:
        out.mkdir()
        mesh.save(out / 'finite-basal-wall.npz')
        np.savez_compressed(out / 'patch-and-ancestry.npz', removed_v1_face_ids=np.flatnonzero(patch),
                            removed_v1_oriented_triangles=old_triangles, outside_v1_face_ids=np.flatnonzero(~patch),
                            new_annulus_faces=cap_faces, new_inner_strip_faces=strip_faces, new_vertex_ids=np.asarray(new_ids),
                            old_vertex_count=len(m.vertices), source_primitive_ids=np.r_[a['primitive_ids'][~patch], np.full(48, -1)],
                            source_face_ids=np.r_[a['source_face_ids'][~patch], np.full(48, -1)],
                            retained_parent_triangles=a['parent_triangles'][~patch], retained_corner_barycentric=a['corner_barycentric'][~patch])
    frozen = Path(contract['frozen_inputs']['path']).parent
    packet = frozen / 'candidate-stage-v2/canonical-annulus-corner-fields.npz'
    fields = {k: v.copy() for k, v in np.load(packet, allow_pickle=False).items()}
    # Different triangle exporters may choose a cyclic first corner. Map the
    # exact same physical source corners; words/UV/signs are never recomputed.
    original_positions = fields['world_positions']
    permutations = np.asarray([[int(np.flatnonzero(np.all(have == point, axis=1))[0]) for point in want]
                               for have, want in zip(original_positions, cap_corners)])
    assert all(tuple(p) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)] for p in permutations)
    fields = {k: np.take_along_axis(v, permutations[(...,)+(None,)*(v.ndim-2)], axis=1) for k, v in fields.items()}
    assert np.array_equal(fields['world_positions'], cap_corners)
    write(out / 'attribute-packet-order-diagnostic.json', dict(source_packet_sha256=sha(packet),
          exact_source_triangle_rows=True, cyclic_corner_permutations=permutations.tolist(),
          geometry_changed=False, raw_normal_words_UV_sign_values_recomputed=False))
    np.savez_compressed(out / 'canonical-annulus-corner-fields.npz', **fields)
    check = audit(mesh, require_closed=True, expected_components=1, expected_genus=[8])
    result = dict(schema='left-basal-new-field-canonical-annulus/v1', complete=True, passed=check['passed'],
                  source_catalogue_sha256=contract['source_catalogue_sha256'], field_sha256=contract['fields']['sha256'],
                  previous_candidate_sha256=sha(old / 'finite-basal-wall.npz'), candidate_sha256=sha(out / 'finite-basal-wall.npz'),
                  fresh_ancestry_selection=ancestry_selection, removed_patch_faces=int(patch.sum()), new_annulus_faces=24,
                  new_inner_strip_faces=24, boundary_loop_lengths=[12, 12], source_positions_moved=0,
                  appended_canonical_controls=len(new_ids), appended_nearest_old_max_m=max(nearest_new, default=0),
                  outside_patch_oriented_triangles_literal=np.array_equal(mesh.triangles[:len(remaining)], m.triangles[~patch]),
                  parent_annulus_exact=True, parent_packet_sha256=sha(packet),
                  blood_classifier_path=str(old / 'blood-classifier.npz'), blood_classifier_sha256=sha(old / 'blood-classifier.npz'),
                  wall_audit=check, actual_nine_port_replay_pending=True, fixed_context_pending=True,
                  anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic()-began)
    result['resumed_attribute_corner_order_only'] = resume_attribute_order_only
    write(out / 'proof.json', result)
    print(json.dumps(dict(output=str(out), passed=result['passed'], elapsed_seconds=result['elapsed_seconds'])))
    assert result['passed'], 'Stop on positive geometry defect; no cleanup loop.'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--resume-attribute-order-only', action='store_true')
    args = parser.parse_args()
    run(args.directory.resolve(), args.resume_attribute_order_only)
