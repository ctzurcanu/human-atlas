"""Bind hand section proposals to current directed source face ancestry.

Original BDA triangle order is independently replayed for previously missing
banks. Existing banks must replay every retained face, including raw control
compaction at constructed joins. No geometry is fitted or installed here.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.course_io import sha, write_json
from geometry_tools.commands.prepare_hand_reference import module
from geometry_tools.algorithms.source_face_pairing import remap_source_pairs
from geometry_tools.algorithms.source_profile_axis import source_profile_graph, complete_source_sections


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['source', 'profiles', 'source-catalogue', 'output']:
        p.add_argument('--'+name, required=True)
    args = p.parse_args(); started = time.monotonic()
    paths = {k: Path(getattr(args, k.replace('-', '_'))).resolve()
             for k in ['source', 'profiles', 'source-catalogue']}
    data = {k: json.loads(v.read_text()) for k, v in paths.items()}
    out = Path(args.output).resolve()
    if out.exists(): raise ValueError('Use a new source evidence stage')
    out.mkdir()
    decoder = module(Path(__file__).resolve().parents[2]/'server/biodigital-registration-landmarks.py',
                     'hand_pairing_decoder')
    original = data['source-catalogue']
    entries = {n['objectId']: (n, matrix) for n, matrix, _ in
               decoder.flatten(original['sceneGraph'], original['records'])}
    banks = {r['partId']: r for r in data['profiles']['rows']}
    joins = {r['firstPart']: r for r in data['source'].get('constructedJoins', [])}
    rows = []
    for row in data['source']['rows']:
        if sha(row['path']) != row['sha256']: raise ValueError('Actual source mesh changed')
        mesh = Mesh.load(row['path']); bank = banks.get(row['partId'])
        if bank:
            if sha(bank['path']) != bank['sha256']: raise ValueError('Recorded section bank changed')
            with np.load(bank['path'], allow_pickle=False) as z:
                old_faces = z['source_faces'].copy()
                old_pairs = z['retained_source_triangle_pairs'].copy()
            evidence = dict(kind='recorded-source-bank', path=bank['path'], sha256=bank['sha256'])
        else:
            node, matrix = entries[row['sourceRecord']['sourceId']]
            ref = original['geometry'][decoder.resource_id(node['geometryId'], 'geometry')]
            with (paths['source-catalogue'].parent/ref['file']).open('rb') as stream:
                stream.seek(ref['offset']); raw = stream.read(ref['length'])
            _, old_faces, version = decoder.decode_bda(raw)
            if np.linalg.det(matrix[:3, :3]) < 0: old_faces = old_faces[:, ::-1].copy()
            # Proposals replay ORIGINAL exporter order, never current row order.
            old_pairs = np.arange(len(old_faces)//2*2).reshape(-1, 2)
            evidence = dict(kind='original-bda-triangle-order-proposals', slice=ref,
                            sliceSha256=__import__('hashlib').sha256(raw).hexdigest(), bdaVersion=version)
        join = joins.get(row['partId'])
        kept = np.array(join['keptSourceFaceIds'] if join else np.arange(len(old_faces)))
        pairing, proof = remap_source_pairs(old_faces, mesh.faces, kept, old_pairs)
        g = source_profile_graph(mesh.vertices, mesh.faces, triangle_pairing=pairing, extend_annuli=True)
        _, physical = np.unique(mesh.vertices, axis=0, return_inverse=True)
        sections = complete_source_sections(g, physical)
        file = out/f'profiles-{len(rows):02}.npz'
        np.savez_compressed(file, retained_source_triangle_pairs=pairing,
                            source_faces=mesh.faces, source_vertices=mesh.vertices)
        rows.append(dict(partId=row['partId'], sourcePath=row['path'], sourceSha256=row['sha256'],
                         path=str(file), sha256=sha(file), sourceFaceAncestryVerified=True,
                         ancestry=proof, pairingOrigin=evidence, profiles=len(g['profiles']),
                         completeSourceControlClosures=len(sections),
                         physicalMouths=len(g['boundary_raw']), anatomicalAcceptance=False))
        if time.monotonic()-started > 18: raise TimeoutError('18-second source ancestry stage')
    report = dict(schema='source-vascular-profile-bank/v2', rows=rows, complete=True,
                  installed=False, geometryChanged=False, installationAllowed=False,
                  pairingPolicy='Hash-bound directed source ancestry; strict source annulus checks still required.',
                  inputs={k: dict(path=str(v), sha256=sha(v)) for k, v in paths.items()},
                  sourceCodeSha256=sha(__file__), elapsedSeconds=time.monotonic()-started)
    write_json(out/'report.json', report)
    print(json.dumps(dict(parts=len(rows), profiles=sum(r['profiles'] for r in rows),
                          remappedJoins=len(joins), elapsedSeconds=report['elapsedSeconds'])))


if __name__ == '__main__': main()
