"""Conform one checked native vessel wall to the installed cardiac cell bank.

All existing tissue interfaces stay reciprocal. Ownership is explicit; source
overlap is not evidence of an anatomically correct or physiologically open join.
"""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.conforming_assembly_ownership import subtract_assembly
from geometry_tools.algorithms.conforming_fields import replay_corner_fields
from geometry_tools.algorithms.surface_role_replay import split_by_source_role

ap = argparse.ArgumentParser()
ap.add_argument('--bank', required=True)
addition = ap.add_mutually_exclusive_group(required=True)
addition.add_argument('--root')
addition.add_argument('--additions')
ap.add_argument('--stage', required=True)
ap.add_argument('--native-intersections',action='store_true')
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
emission = json.loads((h / args.bank / 'emission.json').read_text())
accepted = json.loads((h / args.bank / 'checked-geometry.json').read_text())
assert accepted['passed'] and accepted['complete_bank']
assert digest_file(h / args.bank / 'cells.bin') == emission['bank_sha256']
source = h / emission['source_stage']
meta = json.loads((source / 'ownership.json').read_text())
z = np.load(source / 'arrangement.npz')
positions = z['positions'].astype('<f4').astype(float)
surfaces = {r['id']: z[r['faces']] for r in meta['surfaces']}
roles = {r['id']: r.get('role','') for r in meta['surfaces']}
normals = {r['id']: z[r['normals']] for r in meta['surfaces']}
uvs = {r['id']: z[r['uvs']] for r in meta['surfaces']}
cells = meta['cells']
path = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas = json.loads(path.read_text())
assert atlas['finiteCardiacIntegration']['coupledBankSha256'] == emission['bank_sha256']
origin = np.array(meta['origin'])
count = len(positions)
hosts = list(cells)
if args.additions:
    extra_emission = json.loads((h / args.additions / 'emission.json').read_text())
    extra_checks = json.loads((h / args.additions / 'checked-geometry.json').read_text())
    assert extra_checks['passed'] and extra_checks['complete_bank']
    assert digest_file(h / args.additions / 'cells.bin') == extra_emission['bank_sha256']
    extra_source = h / extra_emission['source_stage']
    extra = json.loads((extra_source / 'ownership.json').read_text())
    assert extra['catalogue_sha256'] == digest_file(path) and extra['origin'] == origin.tolist()
    zz = np.load(extra_source / 'arrangement.npz')
    positions = np.r_[positions, zz['positions'].astype('<f4').astype(float)]
    cutters = list(extra['cells'])
    assert not set(cutters) & set(hosts)
    for r in extra['surfaces']:
        sid = r['id']
        assert sid not in surfaces
        surfaces[sid] = zz[r['faces']] + count
        roles[sid] = r.get('role','')
        normals[sid], uvs[sid] = zz[r['normals']], zz[r['uvs']]
    cells.update(extra['cells'])
else:
    root_meta = json.loads((h / args.root / 'construction.json').read_text())
    assert root_meta['passed'] and root_meta['source_catalogue_sha256'] == digest_file(path)
    root = Mesh.load(h / args.root / 'finite-root.npz')
    pid = root_meta['source_part_id']
    assert pid not in cells
    positions = np.r_[positions, root.vertices - origin]
    sid = pid + '/native-finite-wall'
    surfaces[sid] = root.faces + count
    tri = root.triangles
    face_normal = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    face_normal /= np.linalg.norm(face_normal, axis=1)[:, None]
    normals[sid] = np.repeat(np.rint(face_normal * 32767).astype('<i2')[:, None], 3, axis=1)
    uvs[sid] = np.zeros((len(tri), 3, 2), '<f4')
    cells[pid] = [dict(surface=sid, side=1)]
    roles[sid] = ''
    cutters = [pid]
# Exact-coordinate welding only; duplicate source chart controls must share
# one topological position while retaining their independent corner fields.
positions, inverse = np.unique(positions, axis=0, return_inverse=True)
surfaces = {sid: inverse[fs] for sid, fs in surfaces.items()}
atomic_npz(out / 'source-arrangement.npz', positions=positions)
atomic_json(out / 'source.json', dict(
    bank_sha256=emission['bank_sha256'], root_stage=args.root, additions_bank=args.additions,
    source_part_ids=cutters, catalogue_sha256=digest_file(path)))
print('Joining', len(cutters), 'checked outflow/root cells to', len(hosts), 'checked cardiac cells', flush=True)
result = subtract_assembly(positions, surfaces, cells, hosts, cutters,exact_source_replay=not args.native_intersections)
result,roles=split_by_source_role(result,roles)
q = result['positions']
arrays = dict(positions=q)
definitions = []
words, tex = replay_corner_fields(positions, surfaces, normals, uvs, result)
for i, (sid, fs) in enumerate(result['surfaces'].items()):
    if time.monotonic() - start > 110:
        raise TimeoutError('Vessel-root source-field deadline')
    arrays[f'surface_{i}_faces'] = fs
    arrays[f'surface_{i}_normals'] = words[sid]
    arrays[f'surface_{i}_uvs'] = tex[sid]
    definitions.append(dict(id=sid, faces=f'surface_{i}_faces',
                            normals=f'surface_{i}_normals', uvs=f'surface_{i}_uvs',role=roles[sid]))
contacts = [dict(host=host, vessel=u['neighbor'], surface=u['surface'],
                 triangles=len(result['surfaces'][u['surface']]))
            for host in hosts for u in result['cells'][host] if u.get('neighbor') in cutters]
assert contacts, 'No actual shared host attachment; do not emit a disconnected root'
atomic_npz(out / 'arrangement.npz', **arrays)
atomic_json(out / 'ownership.json', dict(
    constructed=True, numerical_acceptance_pending=True, source_stage=args.stage, exact_source_replay=not args.native_intersections,
    catalogue_sha256=digest_file(path), origin=origin.tolist(), hosts=hosts, cutters=cutters,
    surfaces=definitions, cells=result['cells'], contacts=contacts,
    evidence=result['evidence'], source_parts={p['id']: p for p in atlas['parts'] if p['id'] in cells},
    elapsed_seconds=time.monotonic() - start, installed=False))
print(json.dumps(dict(constructed=True, cells=len(cells), contacts=contacts,
                      elapsed_seconds=time.monotonic() - start)), flush=True)
