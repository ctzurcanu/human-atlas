"""Original representative tongue fascicles; never installs or copies a donor.

These small finite bundles illustrate the four intrinsic fibre directions. They
are not full-volume donor segmentation or verified muscular attachments.
Construction and independent emitted-geometry checking are separate bounded
stages, invoked through geometry_tools.runtime.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np
import igl

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.open3d_backend import diagnostics

OUT = ROOT / '.local-models/geometry-workbench/tongue-completion/intrinsics-v1'
SOURCE = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
FRAME = 'model:local-ta98:assembled-metres'
PAPER = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/'
LIMITATION = (
    'Representative original finite fascicle bundles, not donor-specific whole '
    'muscle segmentation. Bundle diameter, count, terminal caps and local spacing '
    'are illustrative, not measured histology. Native external-muscle envelopes '
    'mark interdigitating territories, not literal solid exclusion volumes. '
    'Complete muscular insertions and fibre continuity remain unverified.'
)


def source():
    a = json.loads(SOURCE.read_text())
    p = next(p for p in a['parts'] if p['id'] == 'ZA:Tongue')
    reader = AtlasReader()
    m = reader.read_part(a, p, frame=FRAME, catalogue_path=SOURCE,
                         catalogue_sha=digest_file(SOURCE))
    return a, m, reader


class Envelope:
    def __init__(self, mesh):
        self.mesh = mesh
        self.tri = mesh.triangles
        q = self.tri[:, :, [0, 2]]
        self.base = q[:, 0]
        mat = np.stack((q[:, 1] - q[:, 0], q[:, 2] - q[:, 0]), axis=2)
        det = np.linalg.det(mat)
        self.valid = np.flatnonzero(abs(det) > 1e-16)
        self.inverse = np.linalg.inv(mat[self.valid])

    def width(self, z):
        points = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            a, b = self.tri[:, i], self.tri[:, j]
            hit = (a[:, 2] <= z) & (b[:, 2] >= z) | (a[:, 2] >= z) & (b[:, 2] <= z)
            hit &= abs(b[:, 2] - a[:, 2]) > 1e-15
            aa, bb = a[hit], b[hit]
            points.extend((aa[:, 0] + (z - aa[:, 2]) * (bb[:, 0] - aa[:, 0]) / (bb[:, 2] - aa[:, 2])).tolist())
        if not points:
            raise ValueError('No native tongue section at this AP coordinate')
        return min(points), max(points)

    def interval(self, x, z):
        uv = np.einsum('nij,nj->ni', self.inverse, np.array([x, z]) - self.base[self.valid])
        bary = np.c_[1 - uv.sum(1), uv]
        hit = np.all(bary >= -1e-10, axis=1)
        yy = np.sort(np.einsum('ni,ni->n', bary[hit], self.tri[self.valid[hit], :, 1]))
        # Deduplicate ray/edge ties only; no mesh controls are welded or moved.
        yy = yy[np.r_[True, np.diff(yy) > 1e-10]]
        if len(yy) != 2:
            raise ValueError(f'Expected one native interior interval at {(x,z)}, got {yy}')
        return yy[0], yy[1]

    def point(self, fraction, z, depth=None, height_fraction=None):
        lo, hi = self.width(z)
        x = (hi + lo) / 2 + fraction * (hi - lo) / 2
        bottom, top = self.interval(x, z)
        if height_fraction is not None:
            y = bottom + height_fraction * (top - bottom)
        elif depth > 0:
            y = top - min(depth, .18 * (top - bottom))
        else:
            y = bottom + min(-depth, .25 * (top - bottom))
        return np.array([x, y, z])


def rounded_bundle(path, radius, env, radial=6):
    path = np.asarray(path)
    # Use actual source distances to bound the illustrative bundle radius.
    d2 = igl.point_mesh_squared_distance(path, env.mesh.vertices, env.mesh.faces)[0]
    rr = np.minimum(radius, .32 * np.sqrt(d2))
    if rr.min() < 2e-5:
        raise ValueError('A fascicle centre is too close to the source boundary')
    rr *= .38 + .62 * np.sin(np.linspace(0, np.pi, len(path))) ** .55
    tangent = np.gradient(path, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    normal = np.cross(tangent[0], [0., 1., 0.])
    if np.linalg.norm(normal) < .1:
        normal = np.cross(tangent[0], [0., 0., 1.])
    normal /= np.linalg.norm(normal)
    vertices, normals = [], []
    for center, t, r in zip(path, tangent, rr):
        normal -= t * np.dot(normal, t)
        normal /= np.linalg.norm(normal)
        binormal = np.cross(t, normal)
        for theta in np.arange(radial) * 2 * np.pi / radial:
            n = normal * np.cos(theta) + binormal * np.sin(theta)
            vertices.append(center + r * n)
            normals.append(n)
    vertices.extend([path[0], path[-1]])
    normals.extend([-tangent[0], tangent[-1]])
    faces = []
    for j in range(len(path) - 1):
        for k in range(radial):
            a = j * radial + k
            b = j * radial + (k + 1) % radial
            c = (j + 1) * radial + k
            d = (j + 1) * radial + (k + 1) % radial
            faces.extend([[a, b, c], [b, d, c]])
    cap0, cap1 = len(vertices) - 2, len(vertices) - 1
    end = (len(path) - 1) * radial
    for k in range(radial):
        faces.extend([[cap0, (k + 1) % radial, k], [cap1, end + k, end + (k + 1) % radial]])
    v, f, n = np.array(vertices), np.array(faces), np.array(normals)
    q = v[f] - v.mean(0)
    volume = np.einsum('ij,ij->i', q[:, 0], np.cross(q[:, 1], q[:, 2])).sum() / 6
    if volume < 0:
        f = f[:, ::-1]
    return v, f, n, dict(centerline_m=path.tolist(), nominal_radius_m=radius,
                         actual_radius_min_m=float(rr.min()), actual_radius_max_m=float(rr.max()))


def construct():
    began = time.monotonic()
    assert (ROOT / 'reports/ta98-tongue-inventory.json').exists(), 'Inventory must precede construction'
    a, native, reader = source()
    assert topology(native)['closed_oriented_manifold'], 'A surface cannot establish a tissue compartment'
    env = Envelope(native)
    paths = {}
    z = np.linspace(.030, .0802, 24)
    paths['superior-longitudinal'] = [np.array([env.point(x, q, depth=.0015) for q in z]) for x in [-.60, -.40, -.20, 0., .20, .40, .60]]
    for side, sign in [('left', 1.), ('right', -1.)]:
        paths['inferior-longitudinal-' + side] = [np.array([env.point(sign * x, q, depth=-.003) for q in z]) for x in [.49, .67]]
        transverse = []
        vertical = []
        for q in np.linspace(.035, .075, 6):
            transverse.append(np.array([env.point(sign * x, q, height_fraction=.54) for x in np.linspace(.04, .78, 12)]))
        for q in np.linspace(.0385, .0785, 6):
            for x in [.31, .59]:
                base = env.point(sign * x, q, height_fraction=.5)
                bottom, top = env.interval(base[0], q)
                vertical.append(np.c_[np.full(8, base[0]), np.linspace(bottom + .002, top - .0034, 8), np.full(8, q)])
        paths['transverse-' + side] = transverse
        paths['vertical-' + side] = vertical
    terms = {'superior-longitudinal': 'A05.1.04.106', 'inferior-longitudinal': 'A05.1.04.107',
             'transverse': 'A05.1.04.108', 'vertical': 'A05.1.04.109'}
    rows = []
    OUT.mkdir(parents=True, exist_ok=True)
    for key, bundles in paths.items():
        vv, ff, nn, bundle_rows = [], [], [], []
        count = 0
        for path in bundles:
            v, f, n, record = rounded_bundle(path, .00046, env)
            vv.append(v); ff.append(f + count); nn.append(n)
            record.update(vertex_start=count, vertex_count=len(v), triangle_start=sum(len(q) for q in ff[:-1]), triangle_count=len(f))
            bundle_rows.append(record)
            count += len(v)
        vertices, faces, normals = np.concatenate(vv), np.concatenate(ff), np.concatenate(nn)
        # Construction is written in Float32-relative coordinates like the viewer.
        origin = np.array([0., 1.625, .055])
        emitted = (vertices - origin).astype(np.float32).astype(float) + origin
        base = key.replace('-left', '').replace('-right', '')
        mesh = Mesh(emitted, faces, FRAME, {'normals': normals}, dict(
            kind='original-representative-intrinsic-fascicles', source_catalogue_sha256=digest_file(SOURCE),
            source_tongue_geometry_sha256=native.geometry_sha, source_inventory_sha256=digest_file(ROOT / 'reports/ta98-tongue-inventory.json'),
            primary_reference=PAPER, anatomical_registration_verified=False,
            source_mesh_vertices_copied=False, proprietary_geometry_copied=False,
            representation='representative-fibre-bundles', limitation=LIMITATION,
            construction_origin_m=origin.tolist(), bundles=bundle_rows))
        target = OUT / (key + '.npz')
        mesh.save(target)
        rows.append(dict(key=key, term=terms[base], mesh=str(target), sha256=digest_file(target), geometry_sha256=mesh.geometry_sha,
                         partId='TA98TONGUE:male:intrinsic:' + key, bundles=len(bundles), triangles=len(faces), vertices=len(vertices),
                         color={'superior-longitudinal':[.72,.30,.29,1.], 'inferior-longitudinal':[.68,.34,.31,1.], 'transverse':[.78,.43,.36,1.], 'vertical':[.62,.28,.32,1.]}[base],
                         source_corners_copied=False, anatomical_acceptance=False, representation='representative-fibre-bundles', coverageLimitation=LIMITATION))
    report = dict(schema='tongue-original-intrinsics/v1', source_catalogue_sha256=digest_file(SOURCE),
                  source_tongue_geometry_sha256=native.geometry_sha, source_bank_sha256=native.source['bank_sha256'], frame=FRAME,
                  primary_reference=PAPER, source_geometry_modified=False, installed=False, anatomical_acceptance=False,
                  meshes=rows, total_meshes=len(rows), total_triangles=sum(r['triangles'] for r in rows),
                  total_npz_bytes=sum(Path(r['mesh']).stat().st_size for r in rows), total_raw_geometry_bytes=sum(r['vertices']*18+r['triangles']*12 for r in rows),
                  all_canonical_muscle_terms=list(terms.values()), complete_muscle_volume_segmentation=False,
                  construction_script_sha256=digest_file(__file__),
                  state='constructed-awaiting-independent-audit-and-rendered-review', remaining_limits=[LIMITATION], elapsed_seconds=time.monotonic()-began)
    assert report['total_meshes'] <= 8 and report['total_triangles'] < 10000
    atomic_json(OUT / 'construction.json', report)
    atomic_json(OUT / 'review-spec.json', dict(frame=FRAME, meshes=[dict(mesh=r['mesh'], name=r['key'], color=r['color']) for r in rows],
                resolution=640, wall_seconds=60, views=[dict(name='front',direction=[0,0,1]),dict(name='right',direction=[1,0,0]),dict(name='top',direction=[0,1,0]),dict(name='right-section',direction=[1,0,0],section=dict(normal=[1,0,0],point=[0,1.625,.055],clear_outer=True))]))
    print(json.dumps({k:report[k] for k in ['total_meshes','total_triangles','total_raw_geometry_bytes','elapsed_seconds','state']}))


def check():
    began = time.monotonic()
    record = json.loads((OUT / 'construction.json').read_text())
    a, native, reader = source()
    assert digest_file(SOURCE) == record['source_catalogue_sha256']
    assert native.geometry_sha == record['source_tongue_geometry_sha256']
    rows, meshes = [], []
    for row in record['meshes']:
        path = Path(row['mesh']); assert digest_file(path) == row['sha256']
        m = Mesh.load(path); topo = topology(m)
        assert topo['closed_oriented_manifold'] and len(topo['components']) == row['bundles']
        assert all(c['closed_orientable_genus'] == 0 for c in topo['components'])
        q = np.concatenate([m.vertices, m.triangles.mean(1), (m.triangles[:,0]+m.triangles[:,1])/2,
                            (m.triangles[:,1]+m.triangles[:,2])/2, (m.triangles[:,2]+m.triangles[:,0])/2])
        winding = igl.winding_number(native.vertices, native.faces, q)
        assert np.all(abs(winding) > .99), (row['key'], winding.min())
        d2 = igl.point_mesh_squared_distance(q, native.vertices, native.faces)[0]
        crossing, outer = pair_contacts(native, m, 0., max_pairs=100000, max_tests=1000000)
        assert outer['complete_within_pair_budget'] and outer['crossing_face_pairs'] == 0 and outer['coplanar_overlap_face_pairs'] == 0
        native_checks = diagnostics(m)
        assert native_checks['self_intersection_pairs'] == 0
        rows.append(dict(key=row['key'], sha256=digest_file(path), geometry_sha256=m.geometry_sha, topology=topo,
                         independent_winding_queries=len(q), minimum_absolute_winding=float(abs(winding).min()),
                         minimum_sampled_envelope_clearance_m=float(np.sqrt(d2).min()), exterior_crossings=outer, cpp_diagnostics=native_checks,
                         numerical_checks_passed=True, anatomical_acceptance=False))
        meshes.append(m)
    cross_rows = []
    for i, first in enumerate(meshes):
        for j, second in enumerate(meshes[:i]):
            _, summary = pair_contacts(first, second, 0., max_pairs=50000, max_tests=1000000)
            cross_rows.append(dict(first=record['meshes'][i]['key'],second=record['meshes'][j]['key'],**summary))
    # Cross-muscle contacts remain evidence, not silently excused as tissue fusion.
    actual_crossings = sum(r['crossing_face_pairs'] + r['coplanar_overlap_face_pairs'] for r in cross_rows)
    report = dict(schema='tongue-intrinsic-numerical-audit/v1', passed=actual_crossings==0 and all(r['complete_within_pair_budget'] for r in cross_rows),
                  construction_sha256=digest_file(OUT/'construction.json'), source_catalogue_sha256=digest_file(SOURCE),
                  source_geometry_modified=False, rows=rows, cross_muscle_contacts=cross_rows, cross_muscle_defect_pairs=actual_crossings,
                  anatomical_acceptance=False, visual_review='pending', complete_muscle_volume_segmentation=False, installed=False,
                  elapsed_seconds=time.monotonic()-began)
    atomic_json(OUT / 'checked-geometry.json', report)
    print(json.dumps(dict(passed=report['passed'],cross_muscle_defect_pairs=actual_crossings,elapsed_seconds=report['elapsed_seconds'],visual_review='pending')))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('stage', choices=['construct','check'])
    ap.add_argument('--output', type=Path, default=OUT)
    args = ap.parse_args()
    OUT = args.output.resolve()
    if not OUT.is_relative_to(ROOT / '.local-models/geometry-workbench/tongue-completion'):
        raise ValueError('Candidates must remain in the local tongue workbench')
    construct() if args.stage == 'construct' else check()
