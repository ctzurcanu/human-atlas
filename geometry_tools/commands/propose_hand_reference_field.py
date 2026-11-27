"""Bounded diagnostic shared hand registration from named bone correspondences.

One smooth field acts on all connected arterial parts. No skeleton moves, no
meshes are installed, and a landmark residual is not anatomical acceptance.
"""
import argparse
import copy
import json
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RBFInterpolator

from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True)
    parser.add_argument('--guides', required=True)
    parser.add_argument('--reference-source', required=True,
                        help='Original source mesh report belonging to the landmark guides; prevents applying a pose twice')
    parser.add_argument('--output', required=True)
    parser.add_argument('--mode', choices=['rbf', 'rigid', 'similarity'], default='rbf')
    parser.add_argument('--arterial-courses', help='Hash-bound matched source/recipient artery course report')
    parser.add_argument('--side', choices=['left', 'right'])
    args = parser.parse_args()
    started = time.monotonic()
    source, guides, out = [Path(x).resolve() for x in [args.source, args.guides, args.output]]
    if out.exists():
        raise ValueError('Use a new diagnostic stage')
    data = json.loads(source.read_text())
    guide = json.loads(guides.read_text())
    arterial = json.loads(Path(args.arterial_courses).read_text()) if args.arterial_courses else None
    if arterial and (args.mode != 'rbf' or not args.side):
        raise ValueError('Arterial course field requires explicit side and RBF mode')
    meshes = {r['partId']: Mesh.load(r['path']) for r in data['rows']}
    for row in data['rows']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Source array changed')
    reference_path = Path(args.reference_source).resolve()
    reference_data = json.loads(reference_path.read_text())
    reference = {}
    for row in reference_data['rows']:
        if row['partId'] not in meshes:
            continue
        if sha(row['path']) != row['sha256']:
            raise ValueError('Landmark source reference changed')
        reference[row['partId']] = Mesh.load(row['path'])
    from geometry_tools.algorithms.landmark_source_binding import verify_landmark_source
    source_binding = verify_landmark_source(meshes, reference, coordinate_tolerance_m=1e-6)
    out.mkdir()
    rows, checks = [], []
    for side in ['left', 'right']:
        selected = [r for r in data['rows'] if r['side'] == side]
        references = [x for r in guide['rows'] if r['side'] == side
                      for x in r['namedBoneCorrespondences']]
        p = np.array([r['sourceCentre_m'] for r in references])
        q = np.array([r['targetCentre_m'] for r in references])
        if arterial and side == args.side:
            from .pose_hand_vessel_courses import smooth_samples
            p = np.concatenate([smooth_samples(np.array(r['sourceAxis']), 12)[0]
                                for r in arterial['courseMatches']])
            q = np.concatenate([smooth_samples(np.array(r['targetAxis']), 12)[0]
                                for r in arterial['courseMatches']])
        # Full proximal context stays fixed. These support controls have no
        # anatomical landmark claim; they stop extrapolation into the arm.
        mean = p.mean(0)
        anchors = np.array([[mean[0]+x, y, mean[2]+z]
                            for x in [-.04, .04] for y in [1.02, 1.12]
                            for z in [-.04, .04]])
        centres = np.r_[p, anchors]
        delta = np.r_[q-p, np.zeros_like(anchors)]
        field = RBFInterpolator(centres, delta, kernel='thin_plate_spline',
                                degree=1, smoothing=1e-8 if arterial else 1e-10)
        rotation, scale, translation = np.eye(3), 1., np.zeros(3)
        if args.mode != 'rbf':
            a, b = p-p.mean(0), q-q.mean(0)
            u, singular, vt = np.linalg.svd(a.T@b)
            sign = np.linalg.det(vt.T@u.T)
            vt[-1] *= sign
            rotation = vt.T@u.T
            if args.mode == 'similarity':
                scale = float((singular*np.array([1, 1, sign])).sum()/(a*a).sum())
                if not .9 <= scale <= 1.3:
                    raise ValueError('Whole-hand reference scale outside declared bounds')
            translation = q.mean(0)-scale*rotation@p.mean(0)

        def mapping(points):
            if args.side and side != args.side:
                return points.copy()
            t = np.clip((1.12-points[:, 1])/.16, 0, 1)
            fade = t*t*t*(10-15*t+6*t*t)
            displacement = field(points) if args.mode == 'rbf' else points@(scale*rotation).T+translation-points
            return points + fade[:, None]*displacement

        lengths = [len(meshes[r['partId']].vertices) for r in selected]
        unique, inverse = np.unique(np.concatenate([meshes[r['partId']].vertices
                                                    for r in selected]), axis=0, return_inverse=True)
        mapped = mapping(unique)
        h = 1e-6
        jac = np.stack([(mapping(unique + np.eye(3)[i]*h) -
                         mapping(unique - np.eye(3)[i]*h))/(2*h)
                        for i in range(3)], axis=2)
        determinants = np.linalg.det(jac)
        if determinants.min() <= 0:
            raise ValueError('Reference field inverts sampled vessel controls')
        normal_map = np.linalg.inv(jac).transpose(0, 2, 1)
        offset = 0
        for i, row in enumerate(selected):
            m = meshes[row['partId']]
            ids = inverse[offset:offset+lengths[i]]
            offset += lengths[i]
            attrs = {k: v.copy() for k, v in m.attributes.items()}
            n = attrs['normals'].astype(float)*attrs.get('normalSigns', np.ones(len(ids)))[:, None]
            n = np.einsum('nij,nj->ni', normal_map[ids], n)
            n /= np.linalg.norm(n, axis=1)[:, None]
            attrs['normals'] = np.rint(n*32767).astype('<i2')
            attrs['normalSigns'] = np.ones(len(ids), 'i1')
            candidate = Mesh(mapped[ids], m.faces, m.frame, attrs, m.source)
            path = out/f'part-{len(rows):02}.npz'
            candidate.save(path)
            rows.append(dict(row, path=str(path), sha256=sha(path), geometrySha256=candidate.geometry_sha,
                             sourcePath=row['path'], sourceSha256=row['sha256']))
        singular = np.linalg.svd(jac, compute_uv=False)
        checks.append(dict(side=side, namedBoneCentres=len(p),
                           controlsAreArterialCourseSamples=bool(arterial and side == args.side),
                           landmarkResidualMax_m=float(np.linalg.norm(mapping(p)-q, axis=1).max()),
                           sampledMinimumJacobian=float(determinants.min()),
                           sampledPrincipalStretchRange=[float(singular.min()), float(singular.max())],
                           continuousInjectivityNotProven=True, mode=args.mode,
                           uniformScale=scale, rotation=rotation.tolist(), translation_m=translation.tolist()))
        if time.monotonic()-started > 18:
            raise TimeoutError('18-second shared registration proposal')
    report = copy.deepcopy(data)
    report.update(rows=rows, sourceReportSha256=sha(source), guidesSha256=sha(guides),
                  landmarkSourceBinding=dict(source_binding, referenceReport=str(reference_path),
                                             referenceReportSha256=sha(reference_path)),
                  sharedReferenceFieldChecks=checks, sourceOnly=False, complete=True,
                  installed=False, installationAllowed=False, skeletonChanged=False,
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started,
                  qualification='Diagnostic shared registration; emitted gauge, contacts and self intersections pending.')
    if arterial:
        report['arterialCourseReport'] = dict(path=str(Path(args.arterial_courses).resolve()),
                                            sha256=sha(args.arterial_courses))
    for key in ['catalogueSha256', 'bankSha256', 'sourcePackedCandidateSha256']:
        if key in report:
            report['input'+key[0].upper()+key[1:]] = report.pop(key)
    for pair in report.get('pairs', []):
        pair.update(literalCoordinatesAndEffectiveNormals=False,
                    packedCoordinatesAndEffectiveNormalWordsLiteral=False,
                    coordinatesUseSharedField=True, effectiveNormalsRequireRecheck=True)
    # Face ancestry is now relative to the immediately preceding packed bank.
    # Existing annuli remain literal; no source faces were replaced in this step.
    report['constructedJoins'] = []
    write_json(out/'report.json', report)
    print(json.dumps(dict(parts=len(rows), checks=checks, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
