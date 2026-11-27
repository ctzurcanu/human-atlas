"""Fit one joined source graph to explicit reviewed distal pose proposals.

Unlike independent part fits, literal shared artery controls are one graph node.
This bounded numerical proposal still needs section, collision and render gates.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.algorithms.profile_network_arap import fit_profile_network


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'left', 'right', 'left-pose', 'right-pose', 'output']:
        p.add_argument('--'+name, required=True)
    args = p.parse_args()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve()
             for k in ['source', 'left', 'right', 'left-pose', 'right-pose']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new saved diagnostic stage')
    out.mkdir()
    started = time.monotonic()
    rows, checks = [], []
    for side in ['left', 'right']:
        selected = [r for r in data['source']['rows'] if r['side'] == side]
        targets = {r['partId']: r for r in data[side]['rows']}
        palm = next(r for r in data[side+'-pose']['rows'] if 'common_palmar' in r['partId'])
        matrix = np.array(palm['componentTransforms'][0]['matrix'])
        meshes, proposal, faces, offset = [], [], [], 0
        for row in selected:
            if sha(row['path']) != row['sha256']:
                raise ValueError('Packed source changed')
            m = Mesh.load(row['path'])
            meshes.append(m)
            if row['partId'] in targets:
                r = targets[row['partId']]
                if sha(r['path']) != r['sha256'] or sha(r['sourcePath']) != r['sourceSha256']:
                    raise ValueError('Declared distal source/target changed')
                a, b = Mesh.load(r['sourcePath']), Mesh.load(r['path'])
                distance, ids = cKDTree(a.vertices).query(m.vertices)
                if distance.max() > .0002:
                    raise ValueError('Packed controls do not replay the declared distal source')
                target = m.vertices + (b.vertices-a.vertices)[ids]
            else:
                q = m.vertices@matrix[:3, :3].T + matrix[:3, 3]
                t = np.clip((1.07-m.vertices[:, 1])/.13, 0, 1)
                w = t*t*t*(10-15*t+6*t*t)
                target = m.vertices + w[:, None]*(q-m.vertices)
            proposal.append(target)
            faces.append(m.faces+offset)
            offset += len(m.vertices)
        raw = np.concatenate([m.vertices for m in meshes])
        source, inverse = np.unique(raw, axis=0, return_inverse=True)
        faces = inverse[np.concatenate(faces)]
        desired = np.zeros_like(source)
        np.add.at(desired, inverse, np.concatenate(proposal))
        desired /= np.bincount(inverse)[:, None]
        edges = np.unique(np.sort(np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]],
                                                  faces[:, [2, 0]]]), axis=1), axis=0)
        fixed = np.flatnonzero(source[:, 1] >= 1.07)
        result = fit_profile_network(source, desired, edges, fixed=fixed,
                                     target_weights=np.full(len(source), .03),
                                     rotation_regularization=0, iterations=24, wall_seconds=6)
        offset = 0
        for row, m in zip(selected, meshes):
            ids = inverse[offset:offset+len(m.vertices)]
            offset += len(m.vertices)
            attrs = {k: v.copy() for k, v in m.attributes.items()}
            normal = attrs['normals'].astype(float)*attrs.get('normalSigns', np.ones(len(ids)))[:, None]
            normal = np.einsum('nij,nj->ni', result['rotations'][ids], normal)
            normal /= np.linalg.norm(normal, axis=1)[:, None]
            attrs.update(normals=np.rint(normal*32767).astype('<i2'), normalSigns=np.ones(len(ids), 'i1'))
            mesh = Mesh(result['centres'][ids], m.faces, m.frame, attrs, m.source)
            path = out/f'part-{len(rows):02}.npz'
            mesh.save(path)
            rows.append(dict(row, path=str(path), sha256=sha(path), geometrySha256=mesh.geometry_sha,
                             sourcePath=row['path'], sourceSha256=row['sha256']))
        checks.append(dict(side=side, fixedControls=len(fixed), fixedControlsLiteral=result['fixed_nodes_literal'],
                           sourceNodes=len(source), finalEnergy=result['energies'][-1],
                           maximumTargetResidual_m=float(result['target_residual_m'].max()),
                           nativeElapsedSeconds=result['elapsed_seconds']))
        if time.monotonic()-started > 18:
            raise TimeoutError('18-second whole hand ARAP proposal')
    report = copy.deepcopy(data['source'])
    report.update(rows=rows, constructedJoins=[], sharedGraphChecks=checks, complete=True,
                  inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  sourceOnly=False, installed=False, installationAllowed=False,
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(parts=len(rows), checks=checks, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
