"""Insert a complementary serosal split into all reciprocal cardiac charts."""
from pathlib import Path
import argparse, json, sys, time
import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.algorithms.coupled_mask_partition import partition_mask

ap = argparse.ArgumentParser()
ap.add_argument('--split', required=True)
ap.add_argument('--stage', required=True)
ap.add_argument('--direct-plc', action='store_true')
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
split = h / args.split
construction = json.loads((split / 'construction.json').read_text())
src = h / construction['source_stage']
meta = json.loads((src / 'ownership.json').read_text())
z = np.load(src / 'arrangement.npz')
p = z['positions'].astype('<f4').astype(float)
surfaces = {r['id']: z[r['faces']] for r in meta['surfaces']}
normals = {r['id']: z[r['normals']] for r in meta['surfaces']}
uvs = {r['id']: z[r['uvs']] for r in meta['surfaces']}
epi = 'TA98HEART:male:epicardium'
region = 'TA98HEART:male:pericardium:subserosa'
uses = meta['cells'][epi]
host = np.concatenate([surfaces[u['surface']] if u['side'] == 1 else
                       surfaces[u['surface']][:, [0, 2, 1]] for u in uses])
parents = [(u['surface'], i) for u in uses for i in range(len(surfaces[u['surface']]))]
mask = np.load(split / 'cells.npz')
v, f = mask['positions'], mask['subserosa_faces']
tri = v[f]
d, which, _ = igl.point_mesh_squared_distance(tri.mean(1), p, host)
parent_tri = p[host[which]]
fn = np.cross(parent_tri[:, 1]-parent_tri[:, 0], parent_tri[:, 2]-parent_tri[:, 0])
fn /= np.linalg.norm(fn, axis=1)[:, None]
residual = abs(np.einsum('fki,fi->fk', tri-parent_tri[:, :1], fn)).max(1)
match = (d < 4e-24) & (residual < 2e-12)
refs = [dict(surface=parents[i][0], source_face=parents[i][1]) if keep else None
        for i, keep in zip(which, match)]
print('Serosal mask:', len(f), 'faces;', int(match.sum()), 'literal carrier-plane fragments', flush=True)
def save(value, name='paired-checkpoint'):
    arrays, rows = dict(positions=value['positions']), []
    for i, (sid, faces) in enumerate(value['surfaces'].items()):
        arrays[f'f{i}'] = faces
        arrays[f'n{i}'] = value['normals'][sid]
        arrays[f'u{i}'] = value['uvs'][sid]
        rows.append(dict(id=sid, faces=f'f{i}', normals=f'n{i}', uvs=f'u{i}'))
    atomic_npz(out / (name+'.npz'), **arrays)
    atomic_json(out / (name+'.json'), dict(cells=value['cells'], surfaces=rows,
                validation_pending=True, source_bank_sha256=construction['source_bank_sha256']))
if args.direct_plc:
    from geometry_tools.algorithms.adopt_tissue_partition import adopt_partition
    calculation = h / 'serosal-calculation-v620'
    original_ids = np.load(calculation / 'tetrahedra.npz')['source_controls']
    result = adopt_partition(p, surfaces, meta['cells'], normals, uvs, host=epi,
                 partition_positions=v, partition_cells={epi: mask['serosa_faces'], region: f},
                 original_control_ids=original_ids, wall_seconds=105.)
else:
    result = partition_mask(p, surfaces, meta['cells'], normals, uvs, [epi], v, f,
                            region=region, mask_face_sources=refs, wall_seconds=105.,
                            output_checkpoint=save,
                            checkpoint=lambda pid, proof: print('Serosal complementary volume', proof, flush=True))
save(result, 'arrangement')
new = json.loads((out / 'arrangement.json').read_text())
meta.update(source_stage=args.stage, cells=new['cells'], surfaces=new['surfaces'],
            serosal_partition=dict(source_bank_sha256=construction['source_bank_sha256'],
                construction=str(split / 'construction.json'),
                construction_sha256=digest_file(split / 'construction.json'),
                original_epicardium=epi, serosa=epi, subserosa=region,
                terms=['A12.1.08.008', 'A12.1.08.009'], evidence=result['evidence'],
                illustrative_histological_ratio=True, no_exterior_expansion=True),
            numerical_acceptance_pending=True, installed=False,
            elapsed_seconds=time.monotonic()-start)
atomic_json(out / 'ownership.json', meta)
print('Shared serosal tissue arrangement saved:', len(result['cells']), 'cells;',
      time.monotonic()-start, 's', flush=True)
