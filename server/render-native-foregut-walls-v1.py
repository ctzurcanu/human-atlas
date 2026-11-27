"""Rasterize actual source/packed triangles and true triangle-plane sections."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/foregut-arterial-batch-v4/native-walls-v1'
parser = argparse.ArgumentParser()
parser.add_argument('stage', choices=['source', 'packed-wall', 'packed-exterior'])
parser.add_argument('--directory', type=Path, default=OUT)
parser.add_argument('--parts', nargs='+', default=['gda', 'coeliac'])
args = parser.parse_args()
OUT = args.directory.resolve()
started = time.monotonic()
image = Image.new('RGB', (1440, 500 * len(args.parts)), '#141d24')
draw = ImageDraw.Draw(image)
bindings = []
for row, short in enumerate(args.parts):
    path = OUT / (short + '-' + args.stage + '.npz')
    mesh = np.load(path)
    v, f = mesh['vertices'], mesh['faces']
    t = v[f]
    source = np.load(OUT / (short + '-source.npz'))
    rings = np.load(OUT / (short + '-rings.npz'))['rings']
    centres = source['vertices'][rings].mean(1)
    centre = v.mean(0)
    radius = np.ptp(v, axis=0).max()
    normals = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
    normals /= np.linalg.norm(normals, axis=1)[:, None]
    shade = .28 + .72 * np.abs(normals @ np.array([.25, .5, .829156]))
    for col, direction in enumerate([[0, 0, 1], [1, .1, 0]]):
        d = np.array(direction, float)
        d /= np.linalg.norm(d)
        r = np.cross([0, 1, 0], d)
        r /= np.linalg.norm(r)
        up = np.cross(d, r)
        p = t - centre
        scale = 420 / radius
        xy = np.stack([col * 450 + 225 + (p @ r) * scale,
                       row * 500 + 265 - (p @ up) * scale], axis=-1)
        for i in np.argsort((p @ d).mean(1)):
            draw.polygon([tuple(x) for x in xy[i]],
                         fill=tuple(int(c * shade[i]) for c in (229, 91, 95)))
        draw.text((col * 450 + 20, row * 500 + 16), short + ' ' + args.stage + (' anterior' if col == 0 else ' lateral'), fill='white')
    for number, station in enumerate([len(rings) // 3, len(rings) * 2 // 3]):
        point = centres[station]
        normal = centres[station + 1] - centres[station - 1]
        normal /= np.linalg.norm(normal)
        axis = np.eye(3)[np.argmin(abs(normal))]
        r = np.cross(axis, normal)
        r /= np.linalg.norm(r)
        up = np.cross(normal, r)
        dist = (t - point) @ normal
        selected = (dist.min(1) < 0) & (dist.max(1) > 0)
        segments = []
        for tri, dd in zip(t[selected], dist[selected]):
            hits = []
            for a, b in [(0, 1), (1, 2), (2, 0)]:
                if dd[a] * dd[b] < 0:
                    hits.append(tri[a] + (tri[b] - tri[a]) * dd[a] / (dd[a] - dd[b]))
            if len(hits) == 2:
                q = np.array(hits) - point
                segments.append(np.stack([q @ r, q @ up], axis=1))
        assert segments
        allp = np.concatenate(segments)
        size = np.ptp(allp, axis=0).max()
        mid = (allp.min(0) + allp.max(0)) / 2
        for segment in segments:
            xy = (segment - mid) * (180 / size) + [1130, row * 500 + number * 230 + 140]
            draw.line([tuple(x) for x in xy], fill='#ffb3aa', width=3)
        draw.text((985, row * 500 + number * 230 + 25), 'Actual transverse section ' + str(station), fill='white')
    bindings.append(dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
path = OUT / (args.stage + '-surface-sections.png')
image.save(path)
(OUT / (args.stage + '-render-binding.json')).write_text(json.dumps(dict(
    meshes=bindings, image=str(path), image_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    actual_triangles_and_sections=True, elapsed_seconds=time.monotonic() - started), indent=2) + '\n')
print(str(path))
