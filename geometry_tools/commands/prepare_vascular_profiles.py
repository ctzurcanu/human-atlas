"""Recover validated annular source profiles and retain explicit branch cores.

Uses original indexed BioDigital triangles, not a guessed cylindrical surface.
Consecutive triangle pairs are proposals; only strict annulus topology may
certify a profile. This stage does not move anatomy or establish lumen patency.
"""
import argparse
import json
import time
from collections import defaultdict, deque
from pathlib import Path
import numpy as np
from geometry_tools.commands.prepare_hand_reference import module
from geometry_tools.course_io import sha, write_json
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_tube_topology import recover_regular_source_annuli
from geometry_tools.algorithms.tube_continuity import boundary_loops


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--context', required=True)
    parser.add_argument('--replay', required=True)
    parser.add_argument('--source-catalogue', required=True)
    parser.add_argument('--contains', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    start = time.monotonic()
    context_path, replay_path = Path(args.context), Path(args.replay)
    context = json.loads(context_path.read_text())
    replay = json.loads(replay_path.read_text())
    if replay['contextSha256'] != sha(context_path) or not replay['sourceFrameVerified']:
        raise ValueError('Source profile extraction needs the actual verified source frame')
    source_path = Path(args.source_catalogue)
    if sha(source_path) != replay['sourceCatalogueSha256']:
        raise ValueError('Original source catalogue changed')
    source = json.loads(source_path.read_text())
    root = Path(__file__).resolve().parents[2]
    decoder = module(root/'server/biodigital-registration-landmarks.py', 'profile_source_decoder')
    entries = {n['objectId']: (n, m) for n, m, _ in decoder.flatten(source['sceneGraph'], source['records'])}
    by_id = {r['partId']: r for r in replay['vascularReplay']}
    selected = [r for r in context['blood'] if args.contains in r['partId']]
    if not selected:
        raise ValueError('No source records match the explicit selection')
    out = Path(args.output).resolve()
    if out.exists():
        raise ValueError('Use a new stage directory')
    out.mkdir()
    rows = []
    for row in selected:
        if sha(row['path']) != row['sha256']:
            raise ValueError('Exported actual mesh changed')
        mesh = Mesh.load(row['path'])
        n, matrix = entries[row['sourceRecord']['sourceId']]
        ref = source['geometry'][decoder.resource_id(n['geometryId'], 'geometry')]
        with (source_path.parent/ref['file']).open('rb') as stream:
            stream.seek(ref['offset'])
            raw = stream.read(ref['length'])
        v, f, _ = decoder.decode_bda(raw)
        if np.linalg.det(matrix[:3, :3]) < 0:
            f = f[:, ::-1].copy()
        source_controls = np.array(by_id[row['partId']]['installedControlToOriginalSource'])
        actual_faces = source_controls[mesh.faces]
        # Preserve the installed face subset: some original exports contain
        # collapsed triangles already removed by earlier source repairs. Map
        # original pair proposals to retained face IDs without adding faces.
        _, aliases = np.unique(v, axis=0, return_inverse=True)
        available = defaultdict(deque)
        for i, face in enumerate(aliases[actual_faces]):
            available[tuple(face)].append(i)
        original_to_actual = np.full(len(f), -1, int)
        for i, face in enumerate(aliases[f]):
            ids = available[tuple(face)]
            if ids:
                original_to_actual[i] = ids.popleft()
        if any(available.values()):
            raise ValueError('Installed triangles do not replay the original source')
        original_pairs = np.arange(len(f)//2*2).reshape(-1, 2)
        pairing = original_to_actual[original_pairs]
        pairing = pairing[(pairing >= 0).all(axis=1)]
        recovered = recover_regular_source_annuli(v, actual_faces, pairing, ring_size=6)
        physical = recovered['raw_to_physical'][source_controls]
        representatives = {int(k): int(np.flatnonzero(physical == k)[0]) for k in np.unique(physical)}
        profiles, excluded, ring_ids = [], [], {}
        for ring_id, ring in enumerate(recovered['rings']):
            if not set(ring) <= set(representatives):
                excluded.append(ring)
                continue
            ring_ids[ring_id] = len(profiles)
            profiles.append([representatives[k] for k in ring])
        links, link_faces = [], []
        for strip in recovered['strips']:
            a, z = strip['rings']
            if a in ring_ids and z in ring_ids:
                links.append([ring_ids[a], ring_ids[z]])
                link_faces.append(strip['source_face_ids'])
        # Boundary loops are original surface boundaries, not invented closures.
        _, rep, inv = np.unique(physical, return_index=True, return_inverse=True)
        mouths = [rep[loop].tolist() for loop in boundary_loops(inv[mesh.faces])]
        path = out/f'profiles-{len(rows):02}.npz'
        np.savez_compressed(path, physical_ids=physical,
                            profiles=np.array(profiles, int).reshape(-1, 6),
                            profile_links=np.array(links, int).reshape(-1, 2),
                            profile_link_faces=np.array(link_faces, int).reshape(-1, 12),
                            retained_source_triangle_pairs=pairing,
                            source_vertices=mesh.vertices, source_faces=mesh.faces)
        rows.append(dict(partId=row['partId'], sourcePath=row['path'], sourceSha256=row['sha256'],
                         path=str(path), sha256=sha(path), profiles=len(profiles),
                         actualProfileLinks=len(links),
                         strictSourceAnnularStrips=len(recovered['strips']),
                         unrecoveredSourceFaceCount=len(recovered['unassigned_source_face_ids']),
                         sourceFaces=len(f), currentFaces=len(mesh.faces),
                         topologyUsesInstalledFaceSubset=True,
                         originalPhysicalMouths=mouths, excludedProfiles=excluded,
                         sourceProfileCoordinatesChanged=False))
        write_json(out/'report.json', dict(schema='source-vascular-profile-bank/v1',
                   contextSha256=sha(context_path), replaySha256=sha(replay_path), rows=rows,
                   complete=len(rows)==len(selected), installed=False,
                   pairingPolicy='Consecutive original triangle-pair proposal; strict source annulus validation required.',
                   sourceFacesAndBranchCoresRetained=True,
                   elapsedSeconds=time.monotonic()-start))
        if time.monotonic()-start > 20:
            raise TimeoutError('20-second native source-profile stage')
    print(json.dumps(dict(parts=len(rows), profiles=sum(r['profiles'] for r in rows),
                          elapsedSeconds=time.monotonic()-start)))


if __name__ == '__main__':
    main()
