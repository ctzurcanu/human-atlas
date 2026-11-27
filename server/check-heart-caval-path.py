"""Check actual exported tissue against the declared caval entry course."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
import igl
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import Mesh, AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.tube_sections import path_tangents, rotation_minimizing_frames
from geometry_tools.algorithms.exact_triangle_separation import exact_segment_triangle_intersection

ap = argparse.ArgumentParser()
ap.add_argument('--bank', required=True)
ap.add_argument('--root', required=True)
args = ap.parse_args()
start = time.monotonic()
h = ROOT / '.local-models/geometry-workbench/heart'
stage = h / args.bank
emission = json.loads((stage / 'emission.json').read_text())
assert digest_file(stage / 'cells.bin') == emission['bank_sha256']
reader = AtlasReader()
mini = dict(chunks=[dict(url='/local-models/' + str((stage / 'cells.bin').relative_to(ROOT / '.local-models')))])
meshes = [reader.read_part(mini, dict(row, chunk=0, geometryOrigin=emission['geometry_origin']),
                          frame='model:local-ta98:assembled-metres') for row in emission['parts'].values()]
offsets = np.r_[0, np.cumsum([len(m.vertices) for m in meshes])]
all_points = np.concatenate([m.vertices for m in meshes])
v, inverse = np.unique(all_points, axis=0, return_inverse=True)
all_faces = np.concatenate([inverse[m.faces + offsets[i]] for i, m in enumerate(meshes)])
_, first, counts = np.unique(np.sort(all_faces, axis=1), axis=0, return_index=True, return_counts=True)
assert counts.max() == 2
faces = all_faces[first[counts == 1]]
boundary = Mesh(v, faces, meshes[0].frame)
tree = BVH(boundary)
root = np.load(h / args.root / 'source-walls.npz')
centers = root['centers']
tangent = path_tangents(centers)
# A short explicit inward continuation must reach the known right blood space.
curve = np.r_[centers[:1] - .006 * tangent[:1], centers]
frames = rotation_minimizing_frames(curve)
fluid = Mesh.load(h / 'blood-spaces-v130/fluid-boundary-3.npz')
reaches_atrium = bool(abs(igl.winding_number(fluid.vertices, fluid.faces, curve[:1])[0]) > .5)
rows = []
tests = 0
for radius in [0., .003]:
    for angle in ([0.] if radius == 0 else np.arange(8) * np.pi / 4):
        track = curve + radius * (np.cos(angle) * frames[:, :, 0] + np.sin(angle) * frames[:, :, 1])
        winding = abs(igl.winding_number(v, faces, track))
        hits = []
        for k, (a, b) in enumerate(zip(track, track[1:])):
            near = tree.query_bounds(np.minimum(a, b), np.maximum(a, b), distance_m=0.)
            tests += len(near)
            if tests > 200000 or time.monotonic() - start > 25:
                raise TimeoutError('Caval path complete-segment query budget')
            for fi in near:
                if exact_segment_triangle_intersection(a, b, v[faces[fi]])['intersects']:
                    hits.append(dict(segment=k, face=int(fi)))
        rows.append(dict(radius_m=radius, angle_rad=float(angle), stations_outside_tissue=bool(np.all(winding < .5)),
                         tissue_crossings=hits, curve_m=track.tolist()))
result = dict(passed=reaches_atrium and all(r['stations_outside_tissue'] and not r['tissue_crossings'] for r in rows),
              bank_sha256=emission['bank_sha256'], boundary_geometry_sha256=boundary.geometry_sha,
              reaches_right_atrium=reaches_atrium, complete_segment_tests=True, tracks=rows,
              tests=tests, elapsed_seconds=time.monotonic() - start,
              scope='Centre course and eight 3 mm offset tracks; not a complete valve flow simulation')
atomic_json(stage / 'caval-path.json', result)
print(json.dumps({k: result[k] for k in ['passed', 'reaches_right_atrium', 'tests', 'elapsed_seconds']}), flush=True)
raise SystemExit(0 if result['passed'] else 2)
