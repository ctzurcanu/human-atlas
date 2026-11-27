"""Bounded actual-triangle crossing screen for a file-backed mesh group.

The caller supplies an explicit context report/key. This is not a containment
test or an anatomical decision; cropped/partition interfaces remain qualified.
"""
import argparse
from collections import Counter
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.course_io import sha, write_json
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate', required=True)
    p.add_argument('--context', required=True)
    p.add_argument('--context-key', default='fixedBonePieces')
    p.add_argument('--output', required=True)
    args = p.parse_args()
    start = time.monotonic()
    candidate, context = Path(args.candidate), Path(args.context)
    data, targets = json.loads(candidate.read_text()), json.loads(context.read_text())[args.context_key]
    output = Path(args.output)
    if output.exists():
        raise ValueError('Use a new screen report')
    moving = {r['partId'] for r in data['rows']}
    triangles, owners, source_rows = [], [], []
    for row in targets:
        if row['partId'] in moving:
            continue
        if sha(row['path']) != row['sha256']:
            raise ValueError('Declared context changed')
        mesh = Mesh.load(row['path'])
        triangles.append(mesh.triangles)
        owners.extend((row.get('sourceObject', row['partId']), row['partId'], i) for i in range(len(mesh.faces)))
        source_rows.append(dict(partId=row['partId'], geometrySha256=mesh.geometry_sha))
        if time.monotonic()-start > 20:
            raise TimeoutError('20-second native context screen')
    t = np.concatenate(triangles)
    fixed = Mesh(t.reshape(-1, 3), np.arange(t.size//3).reshape(-1, 3), mesh.frame)
    tree = BVH(fixed)
    rows = []
    for row in data['rows']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Candidate changed')
        mesh = Mesh.load(row['path'])
        if mesh.frame != fixed.frame:
            raise ValueError('Candidate/context frame mismatch')
        first = BVH(mesh); tests = 0; counts = Counter(); witnesses = []
        for batch in first.candidate_batches(tree, 0., batch_size=512):
            tests += len(batch)
            if tests > 100000:
                raise ValueError('Finite 100000-pair budget reached')
            result = triangle_pairs(first.triangles[batch[:, 0]], tree.triangles[batch[:, 1]])
            for k in np.flatnonzero(np.isin(result['relation'], ['crossing', 'coplanar-overlap'])):
                a, z = map(int, batch[k]); owner, pid, face = owners[z]; counts[owner] += 1
                if len(witnesses) < 100:
                    witnesses.append(dict(owner=owner, contextPart=pid, contextFace=int(face),
                                          candidateFace=a, kind=str(result['relation'][k])))
            if time.monotonic()-start > 20:
                raise TimeoutError('20-second native context screen')
        rows.append(dict(partId=row['partId'], geometrySha256=mesh.geometry_sha,
                         triangleTests=tests, crossingFacePairsByOwner=dict(counts), firstWitnesses=witnesses))
    report = dict(schema='declared-context-group-crossing-screen/v1', candidateSha256=sha(candidate),
                  contextSha256=sha(context), contextKey=args.context_key, contextSources=source_rows,
                  rows=rows, wholeDeclaredSurfaceScreenComplete=True,
                  containmentAndAnatomicalAcceptanceEstablished=False,
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-start)
    write_json(output, report)
    print(json.dumps(dict(parts=len(rows), crossingFacePairs=sum(sum(r['crossingFacePairsByOwner'].values()) for r in rows),
                          elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
