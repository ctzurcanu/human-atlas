"""Read the supplied local model as landmark evidence, never a runtime donor."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'scripts')]
from ta98_glb_geometry import load_roots
from geometry_tools.common import atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.fit import landmark_transform

ap = argparse.ArgumentParser()
ap.add_argument('--stage', required=True)
args = ap.parse_args()
start = time.monotonic()
out = ROOT / '.local-models/geometry-workbench/heart' / args.stage
out.mkdir(exist_ok=False)
directory = ROOT.parent / '3d/ext/meshes'
metadata_path = directory / '964db2dd4f98052f03baa9ca5f2dbcae/assets.json'
metadata = json.loads(metadata_path.read_text())
names = {r['cid']: r['name'] for r in json.loads((directory / 'en.json').read_text())['items']}
associations = {
    24297: 'TA98HEART:male:mitral-anterior-cusp',
    19627: 'ZA:Posterior leaflet of left atrioventricular valve',
    18198: 'TA98HEART:male:tricuspid-anterior-cusp',
    20944: 'ZA:Inferior leaflet of right atrioventricular valve',
    20999: 'ZA:Septal leaflet of right atrioventricular valve',
    18980: 'ZA:Anterior papillary muscle of right ventricle',
    20592: 'ZA:Inferior papillary muscle of right ventricle',
    18331: 'ZA:Septal papillary muscle of right ventricle'}
additional = {18598, 22809, 23777, 20644, 23809, 20454, 19988, 18054,
              22621, 21663, 24354, 23271, 22800, 18687, 16459}
rows = {r['sni'][0]['id']: r for r in metadata['structures']
        if r.get('cid') in set(associations) | additional}
source = ROOT / '.local-models/models/Male_2024.glb'
meshes = load_roots(source, set(rows), allow_zero_weight_targets=True)
by_cid = {rows[m['node']]['cid']: m for m in meshes}
catalogue = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
atlas = json.loads(catalogue.read_text())
parts = {p['id']: p for p in atlas['parts']}
reader = AtlasReader()
def centroid(v, f):
    tri = v[f]
    area = np.linalg.norm(np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0]), axis=1)
    return area @ tri.mean(1) / area.sum()
x, y = [], []
for cid, target in associations.items():
    m = by_cid[cid]
    x.append(centroid(m['unique'], m['welded']))
    native = reader.read_part(atlas, parts[target], frame='body-metres', catalogue_path=catalogue)
    y.append(centroid(native.vertices, native.faces))
matrix, proof = landmark_transform(np.array(x), np.array(y), similarity=True, scale_bounds=(.001, .1))
records = []
for m in meshes:
    cid = rows[m['node']]['cid']
    vertices = m['unique'] @ matrix[:3, :3].T + matrix[:3, 3]
    Mesh(vertices, m['welded'], 'model:local-ta98:assembled-metres').save(out / f'reference-{cid}.npz')
    records.append(dict(cid=cid, name=names[cid], node=m['node'],
        source_structure_id=rows[m['node']]['id'], triangles=len(m['welded']),
        morph_target_count=m['morphTargetCount'], default_morph_weights=m['defaultMorphWeights'],
        reference_only=True, runtime_copy_permitted=False))
atomic_npz(out / 'landmarks.npz', source=np.array(x), target=np.array(y), matrix=matrix)
atomic_json(out / 'registration.json', dict(source_glb=str(source), source_sha256=digest_file(source),
    metadata_sha256=digest_file(metadata_path), source_catalogue_sha256=digest_file(catalogue),
    association_order=list(associations.items()), matrix=matrix.tolist(), registration=proof,
    reference_meshes=records, anatomical_acceptance=False, installed=False,
    elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(reference_count=len(records), registration=proof,
    elapsed_seconds=time.monotonic()-start)), flush=True)
