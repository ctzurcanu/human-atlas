"""Actual packed cell surfaces and annular-interface section, no mesh fitting."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / '.local-models/geometry-workbench/circulatory-completion/splenic-owned-branches-v1/candidate-v1'
start = time.monotonic()
meshes = [np.load(OUT / (role + '-packed.npz')) for role in ['trunk', 'branches']]
triangles = [m['vertices'][m['faces']] for m in meshes]
image = Image.new('RGB', (1500, 1000), '#141d24')
draw = ImageDraw.Draw(image)
allv = np.concatenate([m['vertices'] for m in meshes])
centre = (allv.min(0) + allv.max(0)) / 2
size = np.ptp(allv, axis=0).max()
plan = json.loads((OUT / 'construction-plan.json').read_text())
splitplan = json.loads((OUT.parents[1] / 'native-foregut-source-partition-check-v1/two-cell-hilar-split-plan-v1/split-plan-and-section-proof.json').read_text())
normal = np.array(splitplan['selected_plane']['distal_plane_normal'])
capcount = json.loads((OUT / 'packing-proof.json').read_text())['split_evidence']['interfaceTriangles']
for col, direction in enumerate([[0, .2, 1], [1, .1, -.7]]):
    d = np.array(direction, float)
    d /= np.linalg.norm(d)
    r = np.cross([0, 1, 0], d)
    r /= np.linalg.norm(r)
    up = np.cross(d, r)
    rows = []
    for role, t in enumerate(triangles):
        # Both owners displayed: the reciprocal internal cut is hidden.
        t = t[:-capcount]
        p = t - centre
        xy = np.stack([col * 750 + 375 + (p @ r) * (680 / size),
                       250 - (p @ up) * (680 / size)], axis=-1)
        n = np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0])
        n /= np.linalg.norm(n, axis=1)[:, None]
        light = .28 + .72 * abs(n @ np.array([.3, .4, .8660254]))
        for i in range(len(t)):
            color = (229, 96, 99) if role == 0 else (211, 130, 151)
            rows.append(((p[i] @ d).mean(), xy[i], tuple(int(c * light[i]) for c in color)))
    for _, points, color in sorted(rows, key=lambda q: q[0]):
        draw.polygon([tuple(p) for p in points], fill=color)
    draw.text((col * 750 + 18, 18), 'Actual assembled exterior: internal shared chart hidden', fill='white')
# Paint the literal common interface face-on, including its inner contour.
mesh = meshes[0]
t = mesh['vertices'][mesh['faces'][-capcount:]]
origin = np.array(splitplan['selected_plane']['origin_m'])
r = np.cross(np.eye(3)[np.argmin(abs(normal))], normal)
r /= np.linalg.norm(r)
up = np.cross(normal, r)
q = t - origin
xy = np.stack([q @ r, q @ up], axis=-1)
mid = (xy.min(axis=(0, 1)) + xy.max(axis=(0, 1))) / 2
scale = 360 / np.ptp(xy.reshape(-1, 2), axis=0).max()
screen = (xy - mid) * scale + [400, 750]
for tri in screen:
    draw.polygon([tuple(p) for p in tri], fill='#e19b9a')
draw.text((180, 525), 'Literal shared annulus: lumen remains open', fill='white')
# True transverse intersections through each side of the ownership boundary.
for j, shift in enumerate([-.0007, .0007]):
    point = origin + shift * normal
    lines = []
    for tt in triangles:
        dd = (tt - point) @ normal
        active = (dd.min(1) < 0) & (dd.max(1) > 0)
        for tri, dist in zip(tt[active], dd[active]):
            hits = []
            for a, b in [(0, 1), (1, 2), (2, 0)]:
                if dist[a] * dist[b] < 0:
                    hits.append(tri[a] + (tri[b] - tri[a]) * dist[a] / (dist[a] - dist[b]))
            if len(hits) == 2:
                p = np.array(hits) - point
                lines.append(np.stack([p @ r, p @ up], axis=1))
    points = np.concatenate(lines)
    mid = (points.min(0) + points.max(0)) / 2
    s = 220 / np.ptp(points, axis=0).max()
    for line in lines:
        p = (line - mid) * s + [850 + j * 360, 750]
        draw.line([tuple(q) for q in p], fill='#ffb0b4', width=3)
    draw.text((755 + j * 360, 565), 'True cut at ' + str(shift * 1000) + ' mm', fill='white')
draw.text((25, 960), 'Diagnostic ownership colors only. Source course unchanged; anatomical boundary, collars and organ interfaces remain provisional.', fill='white')
path = OUT / 'candidate-surface-sections.png'
image.save(path)
(OUT / 'render-binding.json').write_text(json.dumps(dict(
    image=str(path), image_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
    bank_sha256=hashlib.sha256((OUT / 'cells.bin').read_bytes()).hexdigest(),
    actual_packed_meshes_and_triangle_sections=True, elapsed_seconds=time.monotonic() - start), indent=2) + '\n')
print(str(path))
