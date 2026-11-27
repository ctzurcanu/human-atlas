"""One bounded, uninstalled native genial-origin patch proposal.

Source upper muscle coordinates and the skeleton remain literal. The small
new origin footprint lies on an actual mandibular face, with opposed cortex
winding. This is a review candidate, never an anatomical acceptance claim.
"""
from pathlib import Path
import hashlib
import json
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
from geometry_tools.algorithms.tube_rim_bridge import bridge_boundaries

SOURCE = ROOT / '.local-models/geometry-workbench/tongue-completion/source-v1'
OUT = ROOT / '.local-models/geometry-workbench/tongue-completion/native-attachment-candidate-v1'
START = time.monotonic()
DEADLINE = START + 50


def check_time():
    if time.monotonic() >= DEADLINE:
        raise TimeoutError('Native attachment candidate exceeded its 50 second native budget')


def load(i):
    return Mesh.load(SOURCE / 'meshes' / f'native-{i:02d}.npz')


def normals(v, f):
    raw = np.zeros_like(v)
    t = v[f]
    vectors = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
    for j in range(3):
        np.add.at(raw, f[:, j], vectors)
    lens = np.linalg.norm(raw, axis=1)
    raw[lens > 0] /= lens[lens > 0, None]
    return np.rint(np.clip(raw, -1, 1)*32767).astype('<i2')


def make_gg(native, bone, bone_face, side):
    # The actual closest source pair is GG face 2170, opposed to this native
    # mandibular facet. The source tip triangle and cortex are not bound boxes.
    removed_face = 2170
    cortex = bone.triangles[bone_face]
    patch_ids = native.faces[removed_face]
    # The general annulus helper rejected the first skew footprint before
    # geometry generation. Instead, replace just this existing origin disk's
    # three controls by their strictly interior barycentric projections onto
    # the witnessed cortex face. No new fan, annulus or muscle solid is added.
    basis = np.column_stack([cortex[0]-cortex[2], cortex[1]-cortex[2]])
    ab = np.linalg.lstsq(basis, (native.vertices[patch_ids]-cortex[2]).T, rcond=None)[0].T
    barycentric = np.column_stack([ab, 1-ab.sum(1)])
    barycentric = np.maximum(barycentric, .03)
    barycentric /= barycentric.sum(1)[:, None]
    footprint = barycentric @ cortex
    v = np.array(native.vertices)
    v[patch_ids] = footprint
    f = native.faces
    start = len(native.vertices)
    new_normals = normals(v, f)
    new_normals[:start] = native.attributes['normals']
    # Only the three retained patch-boundary controls get updated normals;
    # all untouched source normal words remain literally unchanged.
    affected_faces = np.flatnonzero(np.any(np.isin(f, patch_ids), axis=1))
    boundary_ids = np.unique(f[affected_faces])
    updated = normals(v, f)
    new_normals[boundary_ids] = updated[boundary_ids]
    bone_n = np.cross(cortex[1]-cortex[0], cortex[2]-cortex[0])
    bone_n /= np.linalg.norm(bone_n)
    new_normals[patch_ids] = np.rint(-bone_n*32767).astype('<i2')
    uvs = native.attributes['uvs']
    result = Mesh(v, f, native.frame, {'normals': new_normals, 'uvs': uvs}, {
        'record': native.source['record'],
        'proposal': 'Literal native course with a short opposed mandibular origin rim',
        'parent_geometry_sha256': native.geometry_sha,
        'fixed_bone_geometry_sha256': bone.geometry_sha,
        'source_catalogue_sha256': 'cbed9740aca9cd20ca12bf43a051eb7b36f60f90dc1919b64818dde8848f963f',
        'anatomical_acceptance': False,
    })
    evidence = dict(side=side, native_part_id=native.source['record']['id'],
        replaced_source_patch_face_ids=[removed_face], preserved_native_faces=len(f)-len(affected_faces),
        untouched_source_vertex_coordinates_literal=True,
        changed_source_vertex_coordinate_count=len(patch_ids),
        changed_source_vertex_ids=patch_ids.tolist(),
        vertex_displacements_m=np.linalg.norm(v[patch_ids]-native.vertices[patch_ids], axis=1).tolist(),
        new_origin_bone_part=bone.source['record']['id'], new_origin_bone_face=bone_face,
        bone_face_coordinates=cortex.tolist(), bone_face_outward_normal=bone_n.tolist(),
        footprint_barycentric_controls=barycentric.tolist(),
        footprint_coordinates=footprint.tolist(), muscle_contact_normal=(-bone_n).tolist(),
        original_annulus_helper_failure='No nondegenerate monotone annulus exists; no geometry produced by that attempt.',
        topology=topology(result),
        limitations=['Small finite footprint is a source-constrained proposal, not established histological attachment extent.',
                    'Separate native GG hyoid-body penetration remains unmodified in this genial-origin-only proposal.'])
    return result, evidence


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT/'candidate.json').exists():
        raise FileExistsError('Do not replace the one first candidate for this pair')
    bone = load(5)
    rows = []
    for i, face, side in [(1, 11088, 'left'), (2, 3912, 'right')]:
        check_time()
        native = load(i)
        candidate, evidence = make_gg(native, bone, face, side)
        path = OUT / f'genioglossus-{side}.npz'
        candidate.save(path)
        evidence['mesh'] = str(path)
        evidence['mesh_sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        evidence['open3d'] = diagnostics(candidate, self_intersections=True)
        check_time()
        contacts, summary = pair_contacts(candidate, bone, distance_m=0.0,
            max_pairs=25000, max_tests=200000)
        evidence['mandible_contact'] = summary
        np.savez_compressed(OUT / f'genioglossus-{side}-mandible.npz', **contacts)
        evidence['source_triangle_addition'] = len(candidate.faces)-len(native.faces)
        rows.append(evidence)
    report = dict(schema='native-tongue-attachment-first-candidate/v1',
        elapsed_seconds=time.monotonic()-START, native_deadline_seconds=50,
        source_catalogue_sha256='cbed9740aca9cd20ca12bf43a051eb7b36f60f90dc1919b64818dde8848f963f',
        skeleton_modified=False, installed=False, anatomical_acceptance=False, candidates=rows)
    (OUT / 'candidate.json').write_text(json.dumps(report, indent=2)+'\n')
    spec = json.loads((SOURCE.parent / 'native-attachments-v1/review-spec.json').read_text())
    for row in spec['meshes']:
        for side, idx in [('left', 1), ('right', 2)]:
            if row['mesh'].endswith(f'native-{idx:02d}.npz'):
                row['mesh'] = str(OUT / f'genioglossus-{side}.npz')
                row['name'] = f'Genioglossus genial-origin candidate ({side})'
    (OUT / 'review-spec.json').write_text(json.dumps(spec, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
