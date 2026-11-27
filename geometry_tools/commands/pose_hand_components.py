"""Pose reviewed original hand components with anatomical group similarities.

Digit transforms use paired phalangeal length and a stable palm-plane roll;
palm connectors share one palm transform. No mesh is bent, shrunk against an
obstacle, split, capped or installed. Scale is reported explicitly.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.course_io import sha, write_json
from geometry_tools.fit import landmark_transform, transform_mesh
from geometry_tools.model import Mesh


def frame(axis, normal):
    z = axis/np.linalg.norm(axis)
    x = normal-z*(normal@z); x /= np.linalg.norm(x)
    return np.column_stack([x, np.cross(z, x), z])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['binding', 'context', 'guides', 'reference', 'bone-envelopes', 'output']:
        p.add_argument('--'+key, required=True)
    p.add_argument('--side', choices=['left', 'right'], required=True)
    args = p.parse_args(); started = time.monotonic()
    paths = {key: Path(getattr(args, key.replace('-', '_'))).resolve()
             for key in ['binding', 'context', 'guides', 'reference', 'bone-envelopes']}
    binding, context, guides, reference, envelope = [json.loads(q.read_text()) for q in paths.values()]
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new candidate stage')
    if reference['contextSha256'] != sha(paths['context']) or envelope['contextSha256'] != sha(paths['context']):
        raise ValueError('Source/reference or fixed-skeleton context mismatch')
    rows = [r for r in binding['rows'] if r['side'] == args.side]
    if not rows or not all(r['anatomicallyReviewed'] for r in rows):
        raise ValueError('Only explicitly reviewed source component owners may move')
    chosen = [r for r in guides['rows'] if r['side'] == args.side]
    source_palm = np.array([r['namedBoneCorrespondences'][0]['sourceCentre_m'] for r in chosen])
    target_palm = np.array([r['namedBoneCorrespondences'][0]['targetCentre_m'] for r in chosen])
    palm, palm_metrics = landmark_transform(source_palm, target_palm, similarity=True, scale_bounds=(.9, 1.35))
    def normal(x):
        _, _, vh = np.linalg.svd(x-x.mean(0)); n = vh[-1]
        return n if n[2] >= 0 else -n
    old_normal, new_normal = normal(source_palm), normal(target_palm)
    matrices = {'palm': palm}; metrics = {'palm': palm_metrics}
    for group in chosen:
        pairs = group['namedBoneCorrespondences'][1:]
        s = np.array([r['sourceCentre_m'] for r in pairs])
        t = np.array([r['targetCentre_m'] for r in pairs])
        old_frame, new_frame = frame(s[-1]-s[0], old_normal), frame(t[-1]-t[0], new_normal)
        R = new_frame@old_frame.T
        cs, ct = s.mean(0), t.mean(0)
        old, new = s-cs, t-ct
        scale = float(np.sum(new*(old@R.T))/np.sum(old*old))
        if not .9 <= scale <= 1.35:
            raise ValueError('Measured digit scale outside declared biological-guide range')
        M = np.eye(4); M[:3, :3] = scale*R; M[:3, 3] = ct-scale*R@cs
        key = str(group['digit']); matrices[key] = M
        metrics[key] = dict(scale=scale, rotationDeterminant=float(np.linalg.det(R)),
                            phalangeGuideRms_m=float(np.sqrt(np.mean(np.sum((s@M[:3,:3].T+M[:3,3]-t)**2,axis=1)))),
                            sourcePoints=s.tolist(), targetPoints=t.tolist(),
                            axialRollSource='Corresponding source/recipient palm-plane normals')
    sources = {r['partId']: r for r in context['blood']}; output = []
    out.mkdir()
    for pid in dict.fromkeys(r['partId'] for r in rows):
        original = sources[pid]
        if sha(original['path']) != original['sha256']:
            raise ValueError('Original vascular mesh changed')
        mesh = Mesh.load(original['path']); v = mesh.vertices.copy()
        attrs = {k: a.copy() for k, a in mesh.attributes.items()}
        owned = np.zeros(len(v), bool); faces = []; movements = []
        for row in [r for r in rows if r['partId'] == pid]:
            if row['originalSha256'] != original['sha256']:
                raise ValueError('Component source binding changed')
            ids = np.array(row['rawVertexIds'], int)
            if owned[ids].any():
                raise ValueError('Distinct source components share a movable raw control')
            key = str(row['proposedCourses'][0]) if len(row['proposedCourses']) == 1 else 'palm'
            M = matrices[key]; candidate = transform_mesh(mesh, M)
            v[ids] = candidate.vertices[ids]
            for name in attrs:
                attrs[name][ids] = candidate.attributes[name][ids]
            owned[ids] = True; faces.extend(row['sourceFaceIds'])
            fid = np.array(row['sourceFaceIds'], int)
            tri = mesh.vertices[mesh.faces[fid]]; posed = v[mesh.faces[fid]]
            old_length = np.linalg.norm(tri[:,1]-tri[:,0],axis=1)
            new_length = np.linalg.norm(posed[:,1]-posed[:,0],axis=1)
            error = float(np.max(abs(new_length-old_length*metrics[key]['scale'])))
            if error > 1e-12:
                raise ValueError('Source component similarity did not replay exactly')
            movements.append(dict(sourceComponent=row['physicalComponent'], anatomicalGroup=key,
                                   matrix=M.tolist(), scale=metrics[key]['scale'],
                                   scaledEdgeReplayError_m=error))
        if sorted(faces) != list(range(len(mesh.faces))):
            raise ValueError('Reviewed components must own every actual source triangle exactly once')
        q = Mesh(v, mesh.faces, mesh.frame, attrs, dict(sourcePartId=pid))
        path = out/f'part-{len(output):02}.npz'; q.save(path)
        output.append(dict(partId=pid, name=original['name'], system=original['sourceRecord']['system'],
                           path=str(path), sha256=sha(path), sourcePath=original['path'],
                           sourceSha256=original['sha256'], geometrySha256=q.geometry_sha,
                           sourceFacesLiteral=True, sourceComponentShapesPreserved=True, componentTransforms=movements))
        if time.monotonic()-started > 18:
            raise TimeoutError('18-second native component similarity stage')
    report = dict(schema='anatomical-component-hand-pose/v1', side=args.side, rows=output,
                  inputs={key: dict(path=str(q),sha256=sha(q)) for key,q in paths.items()},
                  sourceCodeSha256=sha(__file__), anatomicalGroupMetrics=metrics,
                  explicitUniformScale=True, calibreChangeOnlyByReportedUniformGroupScale=True,
                  skeletonChanged=False, sourceFacesAndUVsRetained=True,
                  installed=False, installationAllowed=False,
                  pending=['Recipient surfaces/skin and intrinsic geometry.',
                           'Shared palmar and wrist interfaces between independently posed anatomical groups.',
                           'Actual rendered pose and Float32 storage.'],
                  elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(side=args.side, sourceObjects=len(output), groups=metrics,
                          elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
