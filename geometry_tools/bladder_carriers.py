"""Recover and verify a three-mouth finite bladder wall for construction.

The radial sign is a source-selection proposal only. Independent closed-cell,
source-chain, occupied-volume and actual section review must accept the result.
Computational fluid disks never belong to emitted anatomy.
"""
from collections import Counter
from pathlib import Path

import numpy as np

from .algorithms.annular_pipe_envelopes import boundary_loops, face_components
from .algorithms.constrained_boundary_disk import constrained_boundary_disk, cancel_opposed_faces
from .algorithms.finite_tissue_boolean import make_solid
from .audit import audit
from .common import atomic_json, atomic_npz
from .model import Mesh


def _chain(vertices, faces):
    result = Counter()
    for face in faces:
        t = [tuple(p) for p in vertices[face]]
        key = min(tuple(t[k:] + t[:k]) for k in range(3))
        reverse = min(tuple((t[:1] + t[:0:-1])[k:] + (t[:1] + t[:0:-1])[:k]) for k in range(3))
        if key < reverse:
            result[key] += 1
        else:
            result[reverse] -= 1
    return {key: value for key, value in result.items() if value}


def recover(mesh, output, *, inner_face_ids=None, mouth_labels=None, cap_forbid_source_diagonals=False):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    v, inv = np.unique(mesh.vertices, axis=0, return_inverse=True)
    f = inv[mesh.faces]; t = v[f]
    normal = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    if np.any(np.linalg.norm(normal, axis=1) == 0):
        raise ValueError('Repair zero-area source triangles before carrier recovery')
    if inner_face_ids is None:
        center = mesh.bounds.mean(0)
        negative = np.flatnonzero(np.einsum('ij,ij->i', normal, t.mean(1) - center) < 0)
        count, labels = face_components(f[negative])
        if count < 1:
            raise ValueError('No inward-facing source-sheet proposal')
        component = int(np.argmax(np.bincount(labels)))
        inner_ids = negative[labels == component]
    else:
        inner_ids = np.asarray(inner_face_ids)
        if inner_ids.ndim != 1 or not np.issubdtype(inner_ids.dtype, np.integer) or not len(inner_ids) or len(np.unique(inner_ids)) != len(inner_ids) or inner_ids.min() < 0 or inner_ids.max() >= len(f):
            raise ValueError('Explicit distinct valid lining face IDs are required')
    inner = f[inner_ids]
    loops = boundary_loops(inner)
    expected = 3 if mouth_labels is None else len(mouth_labels)
    if len(loops) != expected:
        raise ValueError('The source does not expose the explicitly declared number of physiological mouths')
    # Anatomical labels are explicit source-frame proposals, reviewed below.
    if mouth_labels is None:
        neck = int(np.argmin([v[q].mean(0)[1] for q in loops]))
        names = {neck: 'urethral-neck'}
        rest = [i for i in range(3) if i != neck]
        names[min(rest, key=lambda i: v[loops[i]].mean(0)[0])] = 'right-ureteric'
        names[max(rest, key=lambda i: v[loops[i]].mean(0)[0])] = 'left-ureteric'
    else:
        if len(set(mouth_labels)) != expected:
            raise ValueError('Distinct explicit mouth labels are required')
        order = sorted(range(expected), key=lambda i: v[loops[i]].mean(0)[1])
        names = dict(zip(order, mouth_labels))
    fluid_side = inner[:, [0, 2, 1]]
    disks, disk_evidence = [], []
    for loop in boundary_loops(fluid_side):
        ring_edges = {tuple(sorted((int(a),int(b)))) for a,b in zip(loop,np.roll(loop,-1))}
        source_edges = {tuple(sorted((int(a),int(b)))) for tri in f for a,b in zip(tri,np.roll(tri,-1))}
        forbidden = source_edges-ring_edges if cap_forbid_source_diagonals else ()
        cap, proof = constrained_boundary_disk(v, loop, f, forbidden_edges=forbidden)
        disks.append(cap); disk_evidence.append(proof)
    disks = np.concatenate(disks)
    remaining = np.setdiff1d(np.arange(len(f)), inner_ids)
    outer = np.r_[f[remaining], disks]
    fluid = np.r_[fluid_side, disks]
    outer, outer_cancellation = cancel_opposed_faces(outer)
    fluid, fluid_cancellation = cancel_opposed_faces(fluid)
    wall = Mesh(v, f, mesh.frame)
    outer_mesh, fluid_mesh = Mesh(v, outer, mesh.frame), Mesh(v, fluid, mesh.frame)
    meshes = dict(wall=wall, outer=outer_mesh, fluid=fluid_mesh)
    checks = {name: audit(m, require_closed=True, expected_components=1,
                          expected_genus=[expected - 1] if name == 'wall' else [0], max_pairs=1000000)
              for name, m in meshes.items()}
    if not all(c['passed'] for c in checks.values()):
        atomic_json(output / 'failed-checks.json', checks)
        raise ValueError('A recovered source carrier fails the independent full-geometry checks')
    volumes = {}; solids = {}
    for name, m in meshes.items():
        solids[name], _ = make_solid(m.vertices, m.faces)
        volumes[name] = float(solids[name].volume())
    intersection = float((solids['wall'] ^ solids['fluid']).volume())
    residual = abs(volumes['outer'] - volumes['wall'] - volumes['fluid'])
    reconstructed = np.r_[outer, fluid[:, [0, 2, 1]]]
    exact_chain = _chain(v, f) == _chain(v, reconstructed)
    bound = 1e-12 * max(volumes.values())
    if not exact_chain or any(q <= 0 for q in volumes.values()) or intersection > bound or residual > bound:
        raise ValueError('Recovered occupied/fluid cells do not replay the original finite wall')
    atomic_npz(output / 'carriers.npz', vertices=v, wall_faces=f, outer_faces=outer,
               fluid_faces=fluid, inner_side_faces=inner, outer_and_rim_faces=f[remaining],
               computational_disks=disks, inner_source_face_ids=inner_ids,
               frame=np.array(mesh.frame), **{name.replace('-', '_') + '_ring': loops[i] for i, name in names.items()})
    for name, m in meshes.items():
        m.save(output / (name + '.npz'))
    result = dict(schema='bladder-three-port-carriers/v1' if mouth_labels is None else 'explicit-finite-wall-carriers/v1', source_geometry_sha256=mesh.geometry_sha,
                  frame=mesh.frame, source_faces_literal=True, source_positions_literal=True,
                  inner_source_faces=len(inner), mouth_count=expected,
                  mouths=[dict(name=names[i], controls=len(loop), center_m=v[loop].mean(0).tolist())
                          for i, loop in enumerate(loops)], checks=checks, volumes_m3=volumes,
                  fluid_wall_overlap_m3=intersection, volume_replay_residual_m3=residual,
                  source_boundary_chain_exact=exact_chain, computational_disk_evidence=disk_evidence,
                  opposed_boundary_cancellation=dict(outer=outer_cancellation, fluid=fluid_cancellation),
                  computational_disks_for_renderer=False,
                  scope='Recovered source carriers only. Source registration, shape, wall thickness, named coat domains and adjoining interfaces remain unaccepted.',
                  anatomical_acceptance=False, atlas_modified=False)
    atomic_json(output / 'carriers.json', result)
    return result


if __name__ == '__main__':
    import argparse
    import json
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mesh'); p.add_argument('--output', required=True)
    args = p.parse_args(); result = recover(Mesh.load(args.mesh), args.output)
    print(json.dumps({key: result[key] for key in ['mouth_count', 'volumes_m3', 'fluid_wall_overlap_m3', 'source_boundary_chain_exact']}, indent=2))
