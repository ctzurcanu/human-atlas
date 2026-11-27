"""Restore source sections on one joined hand field; fit only branch cores."""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.algorithms.source_tube_topology import recover_regular_source_annuli
from geometry_tools.algorithms.arap_core import arap_core


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate', required=True)
    p.add_argument('--output', required=True)
    args = p.parse_args()
    path, out = Path(args.candidate).resolve(), Path(args.output).resolve()
    data = json.loads(path.read_text())
    if out.exists():
        raise ValueError('Use a new saved profile stage')
    out.mkdir()
    started = time.monotonic()
    rows, checks = [], []
    for side in ['left', 'right']:
        selected = [r for r in data['rows'] if r['side'] == side]
        sources, posed, faces, offset = [], [], [], 0
        for row in selected:
            if sha(row['path']) != row['sha256'] or sha(row['sourcePath']) != row['sourceSha256']:
                raise ValueError('Source/initial pose changed')
            a, b = Mesh.load(row['sourcePath']), Mesh.load(row['path'])
            if not np.array_equal(a.faces, b.faces):
                raise ValueError('Initial field changed the source face bank')
            sources.append(a)
            posed.append(b)
            faces.append(a.faces+offset)
            offset += len(a.vertices)
        raw = np.concatenate([m.vertices for m in sources])
        source, representative, inverse = np.unique(raw, axis=0, return_index=True, return_inverse=True)
        initial = np.concatenate([m.vertices for m in posed])[representative]
        f = inverse[np.concatenate(faces)]
        rec = recover_regular_source_annuli(source, f, np.arange(len(f)//2*2).reshape(-1, 2))
        profiles = [rec['representative_raw_ids'][ring] for ring in rec['rings']]
        target = initial.copy()
        fixed, used, rotations = [], set(), {}
        for ids in profiles:
            if set(map(int, ids)) & used:
                raise ValueError('Recovered full-network profiles overlap')
            used.update(map(int, ids))
            a, b = source[ids], initial[ids]
            u, _, vt = np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)))
            vt[-1] *= np.linalg.det(vt.T@u.T)
            R = vt.T@u.T
            target[ids] = (a-a.mean(0))@R.T+b.mean(0)
            fixed.extend(ids.tolist())
            for i in ids:
                rotations[int(i)] = R
        proximal = np.flatnonzero(source[:, 1] >= 1.12)
        target[proximal] = source[proximal]
        fixed = np.unique(np.r_[fixed, proximal]).astype(int)
        edges = np.unique(np.sort(np.concatenate([f[:, [0, 1]], f[:, [1, 2]], f[:, [2, 0]]]), axis=1), axis=0)
        q, proof = arap_core(source, edges, fixed, target[fixed], initial,
                             maximum_steps=24, wall_seconds=4)
        dp, dq = source[edges[:, 0]]-source[edges[:, 1]], q[edges[:, 0]]-q[edges[:, 1]]
        covariance = np.zeros((len(source), 3, 3))
        outer = np.einsum('ni,nj->nij', dq, dp)
        np.add.at(covariance, edges[:, 0], outer)
        np.add.at(covariance, edges[:, 1], outer)
        u, _, vt = np.linalg.svd(covariance)
        u[:, :, -1] *= np.linalg.det(u@vt)[:, None]
        R = u@vt
        for i, rotation in rotations.items():
            R[i] = rotation
        R[proximal] = np.eye(3)
        gauge = max((float(np.max(abs(np.linalg.norm(q[ids][:, None]-q[ids][None, :], axis=2)-
                                     np.linalg.norm(source[ids][:, None]-source[ids][None, :], axis=2))))
                     for ids in profiles), default=0)
        if gauge > 1e-12:
            raise ValueError('Rigid source profile gauge failed')
        offset = 0
        for row, m in zip(selected, sources):
            ids = inverse[offset:offset+len(m.vertices)]
            offset += len(m.vertices)
            attrs = {k: v.copy() for k, v in m.attributes.items()}
            n = attrs['normals'].astype(float)*attrs.get('normalSigns', np.ones(len(ids)))[:, None]
            n = np.einsum('nij,nj->ni', R[ids], n)
            n /= np.linalg.norm(n, axis=1)[:, None]
            attrs.update(normals=np.rint(n*32767).astype('<i2'), normalSigns=np.ones(len(ids), 'i1'))
            result = Mesh(q[ids], m.faces, m.frame, attrs, m.source)
            file = out/f'part-{len(rows):02}.npz'
            result.save(file)
            rows.append(dict(row, path=str(file), sha256=sha(file), geometrySha256=result.geometry_sha,
                             sourcePath=row['sourcePath'], sourceSha256=row['sourceSha256']))
        checks.append(dict(side=side, rigidSourceProfiles=len(profiles), gaugeError_m=gauge,
                           proximalControlsLiteral=bool(np.array_equal(q[proximal], source[proximal])),
                           sharedSourceNodes=len(source), coreSolve=proof))
        if time.monotonic()-started > 18:
            raise TimeoutError('18-second shared profile/core proposal')
    report = copy.deepcopy(data)
    report.update(rows=rows, rigidProfileChecks=checks, sourceCodeSha256=sha(__file__),
                  initialCandidateSha256=sha(path), complete=True, installed=False,
                  installationAllowed=False, elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(profiles=sum(r['rigidSourceProfiles'] for r in checks),
                          gaugeError_m=max(r['gaugeError_m'] for r in checks),
                          elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
