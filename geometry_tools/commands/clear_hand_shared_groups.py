"""Clear a connected imported hand with one sparse shared-control solve.

All source section and mouth controls remain in their literal shared closure.
This is a bounded proposal; it does not authorize runtime installation.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.source_profile_io import load_source_pairing
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.algorithms.sparse_group_clearance import sparse_group_clearance
from geometry_tools.commands.check_hand_network import combine


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'candidate', 'profiles', 'bones', 'skin-reference', 'output', 'side']:
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--maximum-motion-mm', type=float, default=5.)
    parser.add_argument('--coherent-bone-sections',action='store_true')
    args = parser.parse_args(); started = time.monotonic()
    if args.side not in ['left', 'right']:
        raise ValueError('Explicit hand side required')
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve()
             for k in ['source', 'candidate', 'profiles', 'bones', 'skin-reference']}
    data = {k: json.loads(p.read_text()) for k, p in paths.items()}
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new shared-clearance stage')
    selected = [r for r in data['candidate']['rows'] if r['side'] == args.side]
    source = {r['partId']: r for r in data['source']['rows']}
    meshes, originals = [], []
    for row in selected:
        sr = source[row['partId']]
        if sha(row['path']) != row['sha256'] or sha(sr['path']) != sr['sha256']:
            raise ValueError('Source or proposal changed')
        m, s = Mesh.load(row['path']), Mesh.load(sr['path'])
        if not np.array_equal(m.faces, s.faces):
            raise ValueError('Source face ancestry changed')
        meshes.append(m); originals.append(s)
    physical, representative, inv = np.unique(np.concatenate([m.vertices for m in originals]),
                                              axis=0, return_index=True, return_inverse=True)
    raw = np.concatenate([m.vertices for m in meshes]); points = raw[representative]
    if not np.array_equal(raw, points[inv]):
        raise ValueError('A declared source physical alias separated')
    sections = set(); offset = 0; face_bank = []
    for row, m, s in zip(selected, meshes, originals):
        bank = next((b for b in data['profiles']['rows'] if b['partId'] == row['partId']), None)
        pairing=load_source_pairing(bank,s.faces)
        graph = source_profile_graph(s.vertices, s.faces, triangle_pairing=pairing, extend_annuli=True)
        for cycle in graph['profiles']+graph['boundary_raw']:
            sections.add(tuple(sorted(set(inv[offset+cycle].tolist()))))
        face_bank.append(inv[offset+m.faces]); offset += len(m.vertices)
    groups = []
    for group in map(set, sorted(sections)):
        touching = [i for i, old in enumerate(groups) if old.intersection(group)]
        for i in reversed(touching):
            group.update(groups.pop(i))
        groups.append(group)
    groups = [np.array(sorted(g), int) for g in groups]
    bones = [r for r in data['bones']['rows'] if r['partId'].endswith('.l' if args.side == 'left' else '.r')]
    skin = data['skin-reference']
    for row in [*bones, skin]:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Fixed context changed')
    if not skin['defaultHiddenExcluded']:
        raise ValueError('Actual active skin required')
    bone = SurfaceIndex(combine([Mesh.load(r['path']) for r in bones]))
    skin_index = SurfaceIndex(Mesh.load(skin['path']))
    faces = np.concatenate(face_bank)
    q, proof = sparse_group_clearance(points, faces, groups, bone, skin_index, wall_seconds=10.,
                                      maximum_motion_m=args.maximum_motion_mm*.001,
                                      coherent_bone_sections=args.coherent_bone_sections)
    gauge = max((float(np.max(np.abs(np.linalg.norm(q[g][:, None]-q[g][None], axis=2)-
                                    np.linalg.norm(points[g][:, None]-points[g][None], axis=2))))
                 for g in groups), default=0.)
    if gauge > 1e-12:
        raise ValueError('Rigid source section changed its calibre')
    if time.monotonic()-started >= 20:
        raise TimeoutError('20-second sparse whole-hand clearance stage')
    normals = np.zeros_like(points); tri = q[faces]
    fn = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    np.add.at(normals, faces.ravel(), np.repeat(fn, 3, axis=0))
    norms = np.linalg.norm(normals, axis=1)
    if (norms <= 1e-20).any():
        raise ValueError('Undefined source normal')
    normals = np.rint(normals/norms[:, None]*32767).astype('<i2')
    out.mkdir(); changed = {}; offset = 0
    for row, m in zip(selected, meshes):
        ids = inv[offset:offset+len(m.vertices)]; offset += len(m.vertices)
        attrs = dict(m.attributes, normals=normals[ids], normalSigns=np.ones(len(ids), 'i1'))
        final = Mesh(q[ids], m.faces, m.frame, attrs, m.source)
        path = out/f'part-{len(changed):02}.npz'; final.save(path)
        bank = next((b for b in data['profiles']['rows'] if b['partId'] == row['partId']), None)
        record = dict(row, path=str(path), sha256=sha(path), geometrySha256=final.geometry_sha)
        if bank:
            record['sourceTrianglePairing'] = dict(path=bank['path'], sha256=bank['sha256'],
                                                   array='retained_source_triangle_pairs')
        changed[row['partId']] = record
    report = copy.deepcopy(data['candidate'])
    for key in ['wholeAssembly', 'trials', 'catalogueSha256', 'bankSha256', 'sourcePackedCandidateSha256']:
        report.pop(key, None)
    report.update(rows=[changed.get(r['partId'], r) for r in report['rows']],
                  sparseSharedClearance=dict(side=args.side, sectionGaugeError_m=gauge, proof=proof),
                  localRepairComplete=proof['complete'], complete=False, installed=False,
                  installationAllowed=False, skeletonChanged=False,
                  inputs={k: dict(path=str(p), sha256=sha(p)) for k, p in paths.items()},
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(side=args.side, proof=proof, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
