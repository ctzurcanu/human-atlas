"""Replace selected myocardial support tissue by an owned annular region."""
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
from geometry_tools.algorithms.coupled_mask_partition import partition_mask

ap = argparse.ArgumentParser()
ap.add_argument('--bank', required=True)
ap.add_argument('--masks', required=True)
ap.add_argument('--side', choices=['right', 'left'], required=True)
ap.add_argument('--stage', required=True)
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
emission = json.loads((h / args.bank / 'emission.json').read_text())
assert emission['passed_storage'] and digest_file(h / args.bank / 'cells.bin') == emission['bank_sha256']
source = h / emission['source_stage']
meta = json.loads((source / 'ownership.json').read_text())
z = np.load(source / 'arrangement.npz')
p = z['positions'].astype('<f4').astype(float)
surfaces = {r['id']: z[r['faces']] for r in meta['surfaces']}
normals = {r['id']: z[r['normals']] for r in meta['surfaces']}
uvs = {r['id']: z[r['uvs']] for r in meta['surfaces']}
spec = json.loads((h / args.masks / 'construction.json').read_text())
assert spec['source_bank_sha256'] in [emission['bank_sha256'], meta.get('annular_base_bank_sha256')]
region = next(r for r in spec['records'] if r['side'] == args.side)
mask = Mesh.load(h / args.masks / f'{args.side}-mask.npz')
assert mask.geometry_sha == region['mask_geometry_sha256']
atomic_npz(out / 'source-arrangement.npz', positions=p)
def save(value, name='paired-checkpoint'):
    arrays, rows = dict(positions=value['positions']), []
    for i, (sid, faces) in enumerate(value['surfaces'].items()):
        arrays[f'f{i}'] = faces
        arrays[f'n{i}'] = value['normals'][sid]
        arrays[f'u{i}'] = value['uvs'][sid]
        rows.append(dict(id=sid, faces=f'f{i}', normals=f'n{i}', uvs=f'u{i}'))
    atomic_npz(out / (name + '.npz'), **arrays)
    atomic_json(out / (name + '.json'), dict(cells=value['cells'], surfaces=rows,
        validation_pending=True, source_bank_sha256=emission['bank_sha256']))
result = partition_mask(p, surfaces, meta['cells'], normals, uvs,
    spec['selected_myocardial_parts'], mask.vertices, mask.faces,
    region=region['region_id'], wall_seconds=105., output_checkpoint=save,
    checkpoint=lambda pid, proof: print('Annular partition', args.side, pid, flush=True))
arrays, rows = dict(positions=result['positions']), []
for i, (sid, faces) in enumerate(result['surfaces'].items()):
    arrays[f'surface_{i}_faces'] = faces
    arrays[f'surface_{i}_normals'] = result['normals'][sid]
    arrays[f'surface_{i}_uvs'] = result['uvs'][sid]
    rows.append(dict(id=sid, faces=f'surface_{i}_faces',
        normals=f'surface_{i}_normals', uvs=f'surface_{i}_uvs'))
atomic_npz(out / 'arrangement.npz', **arrays)
meta.update(source_stage=args.stage, cells=result['cells'], surfaces=rows,
    annular_base_bank_sha256=spec['source_bank_sha256'],
    annular_regions=meta.get('annular_regions', []) + [region],
    annular_partition_evidence=meta.get('annular_partition_evidence', []) + [result['evidence']],
    numerical_acceptance_pending=True, installed=False, elapsed_seconds=time.monotonic()-start)
atomic_json(out / 'ownership.json', meta)
print(json.dumps(dict(constructed=True, cells=len(result['cells']),
    region=region['region_id'], elapsed_seconds=time.monotonic()-start)), flush=True)
