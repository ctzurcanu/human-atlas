"""Bounded masks for finite annular tissue within the existing myocardial host.

Masks follow measured peripheral leaflet attachments. They are construction
data, not extra torus solids, lumen caps or a histological thickness claim.
"""
from pathlib import Path
from collections import defaultdict
import argparse
import json
import sys
import time
import numpy as np
import manifold3d as md
from scipy.spatial import ConvexHull

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.finite_tissue_boolean import make_solid

ap = argparse.ArgumentParser()
ap.add_argument('--bank', required=True)
ap.add_argument('--stage', required=True)
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
out = h / args.stage
out.mkdir(exist_ok=False)
e = json.loads((h / args.bank / 'emission.json').read_text())
assert digest_file(h / args.bank / 'cells.bin') == e['bank_sha256']
source = h / e['source_stage']
meta = json.loads((source / 'ownership.json').read_text())
z = np.load(source / 'arrangement.npz')
p = z['positions'].astype('<f4').astype(float)
surfaces = {r['id']: z[r['faces']] for r in meta['surfaces']}
a = json.loads((ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text())
marker = a['finiteCardiacIntegration']
myocardium = [pid for pid in marker['wallPartIds'] if ':myocardium:' in pid]
directions = np.array([v for v in np.ndindex(3, 3, 3) if v != (1, 1, 1)], float) - 1
directions /= np.linalg.norm(directions, axis=1)[:, None]
records = []
for side, key, term in [('right', 'tricuspid-valve', 'A12.1.00.027'), ('left', 'mitral-valve', 'A12.1.00.032')]:
    leaflets = next(g['partIds'] for g in marker['termGroups'] if g['key'] == key)
    vertices = np.unique(np.concatenate([p[surfaces[u['surface']]].reshape(-1, 3)
        for pid in leaflets for u in meta['cells'][pid]]), axis=0)
    center = vertices.mean(0)
    _, _, basis = np.linalg.svd(vertices - center)
    xy = (vertices - center) @ basis[:2].T
    hull = ConvexHull(xy)
    polygon = xy[hull.vertices]
    groups = defaultdict(list)
    contacts, selected_count = [], 0
    for pid in myocardium:
        for use in meta['cells'][pid]:
            if use.get('neighbor') not in leaflets:
                continue
            sid = use['surface']
            tri = p[surfaces[sid]]
            points = (tri.mean(1) - center) @ basis[:2].T
            distance = np.full(len(points), np.inf)
            for first, second in zip(polygon, np.roll(polygon, -1, axis=0)):
                edge = second - first
                t = np.clip((points-first) @ edge / (edge @ edge), 0, 1)
                distance = np.minimum(distance, np.linalg.norm(points-first-t[:, None]*edge, axis=1))
            chosen = np.flatnonzero(distance <= .002)
            selected_count += len(chosen)
            contacts.append(dict(host=pid, leaflet=use['neighbor'], source_surface=sid,
                peripheral_faces=chosen.tolist(), peripheral_band_m=.002))
            for triangle in tri[chosen]:
                tile = tuple(np.floor(triangle.mean(0) / .002).astype(int))
                groups[tile].append(triangle)
    assert groups, 'No measured peripheral myocardial attachment for ' + side
    pieces = []
    for tile, triangles in groups.items():
        if time.monotonic() - start > 30:
            raise TimeoutError('Annular construction mask deadline')
        points = (np.unique(np.concatenate(triangles), axis=0)[:, None] + .0008*directions).reshape(-1, 3)
        hull = ConvexHull(points)
        f = hull.simplices.copy()
        t = points[f]
        flip = np.einsum('ti,ti->t', np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), hull.equations[:, :3]) < 0
        f[flip] = f[flip][:, [0, 2, 1]]
        pieces.append(make_solid(points, f)[0])
    mask = md.Manifold.batch_boolean(pieces, md.OpType.Add)
    assert mask.status() == md.Error.NoError
    mm = mask.to_mesh64()
    mesh = Mesh(mm.vert_properties[:, :3], np.asarray(mm.tri_verts), 'heart-native-local-metres')
    mesh.save(out / f'{side}-mask.npz')
    records.append(dict(side=side, term=term, region_id=f'TA98HEART:male:fibrous:{side}-annulus',
        measured_contacts=contacts, peripheral_triangle_count=selected_count,
        mask_geometry_sha256=mesh.geometry_sha, convex_tiles=len(groups),
        illustrative_tissue_support_m=.0008, finite_region_partition_pending=True,
        complete_histological_annulus_established=False))
atomic_json(out / 'construction.json', dict(source_bank_sha256=e['bank_sha256'], records=records,
    selected_myocardial_parts=myocardium, computational_only=True, anatomical_acceptance=False,
    elapsed_seconds=time.monotonic()-start))
print(json.dumps(dict(masks=len(records), triangles=[r['peripheral_triangle_count'] for r in records],
    elapsed_seconds=time.monotonic()-start)), flush=True)
