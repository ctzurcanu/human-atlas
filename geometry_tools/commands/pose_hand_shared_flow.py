"""Pose a whole imported hand together, then restore original section gauges.

Named donor/recipient bone centres guide the common incremental map. They are
not articular centres or proof of arterial placement. Only actual stored source
sections are rigidly restored; no faces, branches or skeleton parts are removed.
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
from geometry_tools.algorithms.landmark_flow import landmark_flow
from geometry_tools.algorithms.source_profile_axis import source_profile_graph


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'profiles', 'guides', 'output']:
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--side',choices=['left','right'],help='Retain the other side literally')
    args = parser.parse_args(); started = time.monotonic()
    paths = {k: Path(getattr(args, k)).resolve() for k in ['source', 'profiles', 'guides']}
    data = {k: json.loads(p.read_text()) for k, p in paths.items()}
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new saved hand stage')
    meshes = {}
    for row in data['source']['rows']:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Imported source changed')
        meshes[row['partId']] = Mesh.load(row['path'])
    result = dict(meshes); evidence = []
    for side in ['left', 'right']:
        if args.side and side!=args.side:continue
        rows = [r for r in data['source']['rows'] if r['side'] == side]
        lengths = [len(meshes[r['partId']].vertices) for r in rows]
        raw = np.concatenate([meshes[r['partId']].vertices for r in rows])
        points, inv = np.unique(raw, axis=0, return_inverse=True)
        sections = set(); offset = 0
        for row in rows:
            m = meshes[row['partId']]
            bank = next((b for b in data['profiles']['rows'] if b['partId'] == row['partId']), None)
            pairs=load_source_pairing(bank,m.faces)
            graph = source_profile_graph(m.vertices, m.faces, triangle_pairing=pairs, extend_annuli=True)
            for cycle in graph['profiles']+graph['boundary_raw']:
                sections.add(tuple(sorted(set(inv[offset+cycle].tolist()))))
            offset += len(m.vertices)
        # A recovered section may share source controls with an adjoining
        # mouth chart. Move the entire overlap closure rigidly, rather than
        # independently prescribing contradictory positions to those controls.
        sections = [set(s) for s in sorted(sections)]
        merged = []
        for section in sections:
            touching = [i for i, group in enumerate(merged) if group.intersection(section)]
            for i in reversed(touching):
                section.update(merged.pop(i))
            merged.append(section)
        source_section_count = len(sections)
        sections = [np.array(sorted(s), int) for s in merged]
        guides = [n for r in data['guides']['rows'] if r['side'] == side
                  for n in r['namedBoneCorrespondences']]
        old = np.array([r['sourceCentre_m'] for r in guides])
        new = np.array([r['targetCentre_m'] for r in guides])
        mean = old.mean(0)
        anchors = np.array([[mean[0]+x, y, mean[2]+z] for x in [-.04, .04]
                            for y in [1.02, 1.12] for z in [-.04, .04]])
        centres = np.array([points[s].mean(0) for s in sections])
        warped, jac, proof = landmark_flow(np.r_[points, centres], np.r_[old, anchors],
                                          np.r_[new, anchors], wall_seconds=7.)
        q = warped[:len(points)].copy()
        # Proper polar rotations retain source section shape and area exactly.
        U, _, Vt = np.linalg.svd(jac[len(points):]); R = U@Vt
        if np.linalg.det(R).min(initial=1) <= 0:
            raise ValueError('A section rotation would reflect its source')
        for i, section in enumerate(sections):
            q[section] = (points[section]-centres[i])@R[i].T+warped[len(points)+i]
        gauge_error = max((float(np.max(np.abs(
            np.linalg.norm(q[s][:, None]-q[s][None], axis=2)-
            np.linalg.norm(points[s][:, None]-points[s][None], axis=2)))) for s in sections), default=0.)
        offset = 0
        for row, length in zip(rows, lengths):
            m = meshes[row['partId']]; ids = inv[offset:offset+length]; offset += length
            # Source UVs/material attributes stay literal. Packed normal words
            # are rebuilt from actual final triangles and synchronized later.
            normals = np.zeros_like(points)
            t = q[ids][m.faces]
            fn = np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0])
            np.add.at(normals, ids[m.faces.ravel()], np.repeat(fn, 3, axis=0))
            normal = normals[ids]; norm = np.linalg.norm(normal, axis=1)
            if (norm <= 1e-20).any():
                raise ValueError('Final source normal is undefined')
            attrs = dict(m.attributes, normals=np.rint(normal/norm[:, None]*32767).astype('<i2'),
                         normalSigns=np.ones(length, 'i1'))
            result[row['partId']] = Mesh(q[ids], m.faces, m.frame, attrs, m.source)
        evidence.append(dict(side=side, rigidOverlapGroups=len(sections),
                             sourceSections=source_section_count,
                             maximumRigidGroupControls=max(map(len, sections), default=0), gaugeError_m=gauge_error,
                             commonFlow=proof, allSourceFacesRetained=True))
        if time.monotonic()-started > 20:
            raise TimeoutError('20-second whole-hand source flow deadline')
    out.mkdir(); emitted = []
    for row in data['source']['rows']:
        m = result[row['partId']]; path = out/f'part-{len(emitted):02}.npz'; m.save(path)
        emitted.append(dict(row, path=str(path), sha256=sha(path), geometrySha256=m.geometry_sha,
                            sourcePath=row['path'], sourceSha256=row['sha256']))
    report = copy.deepcopy(data['source'])
    for key in ['wholeAssembly', 'trials', 'catalogueSha256', 'bankSha256', 'sourcePackedCandidateSha256']:
        report.pop(key, None)
    report.update(rows=emitted, sharedFlow=evidence, complete=False, installed=False,
                  installationAllowed=False, skeletonChanged=False, sourceAnnularPropagation=True,
                  pending=['Complete triangle, body packing, joins, storage and render review'],
                  inputs={k: dict(path=str(p), sha256=sha(p)) for k, p in paths.items()},
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    for pair in report['pairs']:
        pair.update(coordinatesLiteral=True, effectiveNormalWordsRequireSync=True,
                    literalCoordinatesAndEffectiveNormals=False)
    write_json(out/'report.json', report)
    print(json.dumps(dict(checks=evidence, elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
