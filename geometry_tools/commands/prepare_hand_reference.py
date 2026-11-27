"""Export paired source/recipient hand bones and replay the vascular carrier.

This is a read-only landmark preparation stage. It neither moves the recipient
skeleton nor modifies or installs vessels. Source meshes retain original faces.
Run under geometry_tools.runtime; the native stage also has a 20-second limit.
"""
from pathlib import Path
import argparse
import importlib.util
import json
import time
from collections import Counter
import numpy as np
from scipy.spatial import cKDTree
from geometry_tools.course_io import sha, write_json
from geometry_tools.model import Mesh


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--context', required=True)
    parser.add_argument('--workbench', required=True)
    parser.add_argument('--source-catalogue', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.monotonic()
    root = Path(__file__).resolve().parents[2]
    b = Path(args.workbench).resolve()
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Keep earlier evidence; use a new output directory')
    context_path = Path(args.context).resolve()
    context = json.loads(context_path.read_text())
    source_path = Path(args.source_catalogue).resolve()
    source = json.loads(source_path.read_text())
    decoder_path = root/'server/biodigital-registration-landmarks.py'
    decoder = module(decoder_path, 'hand_reference_decoder')
    field_path = b/'global-envelope-v1/body-envelope-map-proposal-v2/body_field.py'
    body = module(field_path, 'hand_reference_carrier')
    contract_path = field_path.with_name('contract.json')
    transform_path = b/'global-envelope-v1/global-skin-similarity-v1/transform.json'
    contract = json.loads(contract_path.read_text())
    transform = np.array(json.loads(transform_path.read_text())['combined_Gnew'])
    entries = {n['objectId']: (n, matrix) for n, matrix, _ in
               decoder.flatten(source['sceneGraph'], source['records'])}

    def decode(sid):
        node, matrix = entries[sid]
        ref = source['geometry'][decoder.resource_id(node['geometryId'], 'geometry')]
        with (source_path.parent/ref['file']).open('rb') as stream:
            stream.seek(ref['offset'])
            raw = stream.read(ref['length'])
        v, f, version = decoder.decode_bda(raw)
        p = v.astype(float)@matrix[:3, :3].T+matrix[:3, 3]
        p = p@transform[:3, :3].T+transform[:3, 3]
        p, _, _ = body.apply(p, contract)
        if np.linalg.det(matrix[:3, :3])*np.linalg.det(transform[:3, :3]) < 0:
            f = f[:, ::-1].copy()
        return p, f, dict(sourceId=sid, name=node['displayName'],
                         originalSlice=ref, originalSliceSha256=__import__('hashlib').sha256(raw).hexdigest(),
                         bdaVersion=version, hierarchyMatrix=matrix.tolist())

    # Verify that the source bone guide and the installed vessels really share
    # this carrier. A plausible transform name is not sufficient evidence.
    replay = []
    for row in context['blood']:
        current = Mesh.load(row['path'])
        if sha(row['path']) != row['sha256']:
            raise ValueError('Changed exported vascular source')
        p, f, meta = decode(row['sourceRecord']['sourceId'])
        distance, source_ids = cKDTree(p).query(current.vertices)
        error = float(distance.max())
        _, aliases = np.unique(p, axis=0, return_inverse=True)
        faces_equal = bool(np.array_equal(f, current.faces))
        source_faces = Counter(map(tuple, aliases[f].tolist()))
        installed_faces = Counter(map(tuple, aliases[source_ids[current.faces]].tolist()))
        subset = all(count <= source_faces[face] for face, count in installed_faces.items())
        replay.append(dict(partId=row['partId'], maximumCoordinateError_m=error,
                           decodedVertexCount=len(p), installedVertexCount=len(current.vertices),
                           decodedFaceCount=len(f), installedFaceCount=len(current.faces),
                           originalFacesLiteral=faces_equal,
                           installedFacesAreOriginalSubset=subset,
                           installedControlToOriginalSource=source_ids.tolist(),
                           historicalMissingSourceFaceCount=len(f)-len(current.faces),
                           sourceGeometrySha256=current.geometry_sha))
        # Historical removal of source faces does not alter the registered
        # controls. Record it explicitly; do not restore or discard any face.
        if not subset or error > 1e-7:
            out.mkdir()
            write_json(out/'failed-replay.json', dict(rows=replay, complete=False,
                       sourceFrameVerified=False, elapsedSeconds=time.monotonic()-started))
            raise ValueError('Hand carrier does not replay installed vessel: '+json.dumps(replay[-1]))
        if time.monotonic()-started > 20:
            raise TimeoutError('20-second source reference deadline')

    ordinal = ['First', 'Second', 'Third', 'Fourth', 'Fifth']
    digit_names = ['thumb', 'index_finger', 'middle_finger', 'ring_finger', 'little_finger']
    rows = []
    out.mkdir()
    for side, suffix in [('left', '.l'), ('right', '.r')]:
        for digit, (number, label) in enumerate(zip(ordinal, digit_names), 1):
            segments = [('metacarpal', f'{side}_{number.lower()}_metacarpal',
                         f'{number} metacarpal bone{suffix}')]
            for segment in ['proximal', 'middle', 'distal']:
                if digit == 1 and segment == 'middle':
                    continue
                segments.append((segment, f'{segment}_phalanx_of_{side}_{label}',
                                 f'{segment.capitalize()} phalanx of {number.lower()} finger of hand{suffix}'))
            for segment, sid_suffix, target_owner in segments:
                sid = 'human_20_male_skeletal_system-'+sid_suffix+'_ID'
                p, f, meta = decode(sid)
                targets = [r for r in context['fixedBonePieces'] if r['sourceObject'] == target_owner]
                if not targets:
                    raise ValueError('Missing recipient source owner: '+target_owner)
                for target in targets:
                    if sha(target['path']) != target['sha256']:
                        raise ValueError('Changed fixed bone context')
                path = out/f'{side}-{digit}-{segment}.npz'
                mesh = Mesh(p, f, 'model:local-ta98:assembled-metres', source=meta)
                mesh.save(path)
                rows.append(dict(side=side, digit=digit, segment=segment,
                                 sourcePath=str(path), sha256=sha(path),
                                 sourceGeometrySha256=mesh.geometry_sha, source=meta,
                                 targetOwner=target_owner, targetPieces=targets))
                if time.monotonic()-started > 20:
                    raise TimeoutError('20-second source reference deadline')
    assert len(rows) == 38
    report = dict(schema='paired-hand-reference/v1', sourceCatalogueSha256=sha(source_path),
                  contextSha256=sha(context_path), decoderSha256=sha(decoder_path),
                  carrierSha256=sha(field_path), contractSha256=sha(contract_path),
                  transformSha256=sha(transform_path), vascularReplay=replay, rows=rows,
                  sourceFrameVerified=True, skeletonChanged=False, vesselsChanged=False,
                  anatomicalAcceptance=False, elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(rows=len(rows), vascularReplay=len(replay),
                          maximumReplayError_m=max(r['maximumCoordinateError_m'] for r in replay),
                          elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__':
    main()
