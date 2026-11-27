"""Check the emitted bladder wall against the unchanged current pelvic meshes.

The database supplies a broad phase only. Each retained neighbour is tested on
its actual triangles. Crossings are evidence for review, not inferred bonds.
"""
from pathlib import Path
import json

import numpy as np

from .algorithms.constrained_boundary_disk import cancel_opposed_faces
from .audit import audit
from .bladder_workbench import crop
from .common import atomic_json, digest_file
from .contact import contacts
from .db import Database
from .model import Mesh
from .registry import load_registered
from .review import prepare
from .spatial import neighbors


def emitted_boundary(folder):
    folder = Path(folder)
    manifest = json.loads((folder / 'tissue.json').read_text())
    if digest_file(manifest['bank']) != manifest['bank_sha256']:
        raise ValueError('Emitted tissue bank changed')
    bank = Path(manifest['bank']).read_bytes()
    vertices, faces, count = [], [], 0
    for row in manifest['parts'].values():
        n = row['vertexCount']
        p = np.ndarray((n, 3), '<f4', buffer=bank, offset=row['positions']).astype(float)
        p += np.asarray(row['geometryOrigin'])
        f = np.ndarray((row['indexCount'] // 3, 3), '<u4', buffer=bank, offset=row['indices'])
        vertices.append(p); faces.append(f.astype(int) + count); count += n
    p, mapping = np.unique(np.concatenate(vertices), axis=0, return_inverse=True)
    f, cancelled = cancel_opposed_faces(mapping[np.concatenate(faces)])
    return Mesh(p, f, manifest['frame']), dict(
        tissue_manifest_sha256=digest_file(folder / 'tissue.json'),
        bank_sha256=manifest['bank_sha256'],
        internal_triangle_pairs_cancelled=cancelled['opposedPairsCancelled'],
        identical_coordinates_only=True)


def check(folder, output):
    folder, output = Path(folder).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    wall, evidence = emitted_boundary(folder)
    wall.save(output / 'emitted-wall.npz')
    whole_check = audit(wall, require_closed=True, expected_components=1,
                        expected_genus=[2], max_pairs=1000000)
    if not whole_check['passed']:
        raise ValueError('Emitted whole bladder fails geometry checks')
    db = Database(); db.assert_current_catalogues('local-ta98')
    mid = db.register_artifact_mesh(wall, output / 'emitted-wall.npz',
        part_id='BLADDER:male:pelvic-check:' + wall.geometry_sha,
        name='Emitted male bladder candidate', metadata=evidence)
    broad = neighbors(db, mid, distance_m=.006, limit=1000)
    if broad['truncated']:
        raise ValueError('The pelvic neighbour inventory is truncated')
    rows, entries = [], [dict(mesh=str(output / 'emitted-wall.npz'),
                            name='New bladder wall', color=[.76, .43, .38, 1])]
    bounds = wall.bounds + [[-.018, -.025, -.018], [.018, .025, .018]]
    for candidate in broad['neighbors']:
        pid = candidate['part_id']
        if pid.startswith('TA98PART:ZA:Urinary bladder:'):
            continue  # These are the replaced original, not pelvic neighbours.
        record, m = load_registered(db, candidate['id'])
        _, result = contacts(db, wall, m, 0., max_pairs=1000000)
        rows.append(dict(part_id=pid, name=candidate['name'], system=candidate['system'],
                         geometry_sha256=m.geometry_sha, result=result))
        n = candidate['name'].lower()
        show = (pid == 'ZA:Prostate' or 'seminal gland' in n or
                'prostatic urethra' in n or pid.endswith('ureter-wall:adventitia') or
                pid.endswith('rectum:serosa') or 'pubic symphysis' in n or
                'pubococcygeus' in n or 'hip bone' in n)
        if show:
            clipped = crop(m, bounds)
            if clipped is not None:
                p = output / ('context-' + str(len(entries)) + '.npz'); clipped.save(p)
                color = ([.72, .72, .59, 1] if candidate['system'] == 'skeletal' else
                         [.42, .56, .78, 1] if 'ureter' in n else
                         [.54, .42, .64, 1] if 'rectum' in n else [.72, .62, .35, 1])
                entries.append(dict(mesh=str(p), name=candidate['name'] + ' (uncapped crop)', color=color))
        atomic_json(output / 'contacts-partial.json', rows, overwrite=True)
    crossing = [q['part_id'] for q in rows if q['result']['crossing_face_pairs'] or
                q['result']['coplanar_overlap_face_pairs']]
    result = dict(schema='bladder-pelvic-crossings/v1', source=evidence,
                  whole_check=whole_check, neighbors_tested=len(rows),
                  broad_phase_margin_m=.006, rows=rows, crossing_items=crossing,
                  complete_triangle_scan=all(not q['result']['truncated'] for q in rows),
                  containment_not_tested=True, skeleton_modified=False,
                  anatomical_acceptance=False, installed=False)
    center = wall.bounds.mean(0)
    result['review'] = prepare(db, dict(frame=wall.frame, meshes=entries, resolution=1000,
        views=[dict(name='front-pelvis', direction=[.2, .1, 1]),
               dict(name='right-pelvis', direction=[1, .1, 0]),
               dict(name='posterior-ports', direction=[-.2, .1, -1]),
               dict(name='sagittal-pelvis', direction=[1, .05, 0],
                    section=dict(normal=[1, 0, 0], point=center.tolist(), clear_outer=True))]))
    atomic_json(output / 'pelvis.json', result)
    return result


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder'); p.add_argument('--output', required=True)
    args = p.parse_args(); r = check(args.folder, args.output)
    print(json.dumps({k: r[k] for k in ['neighbors_tested', 'crossing_items', 'complete_triangle_scan', 'review']}, indent=2))
