"""One source-bound coronal-contour loft proposal; never install or accept it.

The selected v3 carrier cross sections are the construction controls. A tensor
product cubic B-spline rounds their polygonal interpolation without inventing
halfspace anatomy. No normal-only repair is claimed. All geometry and pair
checks are saved in separate finite stages. A failed stage remains diagnostic.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path
import shutil
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'server')]
import numpy as np
from scipy.interpolate import BSpline
from geometry_tools.audit import audit, topology
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.constraints import check_constraint
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH

WORK = ROOT / '.local-models/geometry-workbench/tongue-completion'
FRAME = 'model:local-ta98:assembled-metres'
SOURCE = WORK / 'coupled-male-v3-emitted-review'
WHOLE = WORK / 'basal-whole-protected-v1/whole-tongue.npz'
ORIGIN = np.array([0., 1.625, 0.])
PARTS = {'superior-longitudinal': 'part-02.npz',
         'inferior-longitudinal-left': 'part-03.npz',
         'inferior-longitudinal-right': 'part-04.npz',
         'transverse-vertical-shared': 'part-05.npz'}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def section(mesh, z, allow_multiple=False):
    """Return the actual single oriented polygon, retaining triangle ancestry."""
    t = mesh.triangles
    d = t[:, :, 2] - z
    ids = np.flatnonzero((d.min(1) < 0) & (d.max(1) > 0))
    points, edges, rows = {}, [], []
    for face in ids:
        hits = []
        for k in range(3):
            j = (k + 1) % 3
            if d[face, k] * d[face, j] >= 0:
                continue
            # Canonical edge order gives adjacent faces literally identical
            # intersections without quantizing away tiny genuine source faces.
            pa, pb = t[face, k], t[face, j]
            if tuple(pa) > tuple(pb):
                pa, pb = pb, pa
            u = (z-pa[2])/(pb[2]-pa[2])
            xyz = pa+u*(pb-pa)
            key = tuple(xyz[:2])
            points.setdefault(key, xyz[:2])
            hits.append(key)
        require(len(hits) == 2 and hits[0] != hits[1], 'Ambiguous source section triangle')
        edges.append(tuple(hits))
        rows.append(dict(source_face=int(face), edge_keys=[list(x) for x in hits]))
    adjacent = {key: [] for key in points}
    for a, b in edges:
        adjacent[a].append(b)
        adjacent[b].append(a)
    require(adjacent and all(len(v) == 2 for v in adjacent.values()),
            'Source section is not a collection of simple closed loops')
    unused, loops = set(points), []
    while unused:
        start = min(unused)
        loop, previous, current = [], None, start
        while current not in loop:
            loop.append(current)
            unused.discard(current)
            nxt = [x for x in adjacent[current] if x != previous]
            previous, current = current, nxt[0]
        require(current == start, 'Source section has inconsistent connectivity')
        loops.append(np.asarray([points[x] for x in loop]))
    if allow_multiple:
        return loops, rows
    require(len(loops) == 1,
            'Multiple source section loops at Z='+str(z)+' need anatomical correspondence, not automatic bridging')
    polygon = loops[0]
    following = np.roll(polygon, -1, axis=0)
    area = (polygon[:, 0]*following[:, 1]-polygon[:, 1]*following[:, 0]).sum() / 2
    if area < 0:
        polygon = polygon[::-1]
        area = -area
    return polygon, float(area), rows


def resample(polygon, count=32):
    # A rightmost source point fixes phase; subsequent stations use the least
    # squared cyclic correspondence. No rotation of the physical cross section.
    polygon = np.roll(polygon, -int(np.argmax(polygon[:, 0])), axis=0)
    closed = np.r_[polygon, polygon[:1]]
    length = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    cumulative = np.r_[0., np.cumsum(length)]
    s = np.arange(count) * cumulative[-1] / count
    return np.c_[np.interp(s, cumulative, closed[:, 0]),
                 np.interp(s, cumulative, closed[:, 1])]


def periodic_cubic(controls, count=32):
    # Uniform periodic cubic B-spline: each point is a convex combination of
    # four actual resampled section controls, rather than an unconstrained fit.
    count_controls = len(controls)
    phase = np.arange(count) * count_controls / count
    i = np.floor(phase).astype(int)
    u = phase - i
    weights = np.c_[(1-u)**3, 3*u**3-6*u*u+4,
                    -3*u**3+3*u*u+3*u+1, u**3] / 6
    return sum(weights[:, k, None] * controls[(i+k-1) % count_controls]
               for k in range(4))


def quality(mesh):
    tri = mesh.triangles
    raw = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    area = np.linalg.norm(raw, axis=1) / 2
    normal = raw / np.maximum(np.linalg.norm(raw, axis=1)[:, None], 1e-30)
    edges = np.concatenate([np.linalg.norm(tri[:, (i+1) % 3]-tri[:, i], axis=1)
                            for i in range(3)])
    planes = np.c_[np.round(normal, 5),
                   np.round(np.einsum('ij,ij->i', normal, tri[:, 0]), 7)]
    _, inverse = np.unique(planes, axis=0, return_inverse=True)
    return dict(triangles=len(mesh.faces), indexed_vertices=len(mesh.vertices),
                edge_median_m=float(np.median(edges)), edge_95_m=float(np.quantile(edges, .95)),
                largest_triangle_area_m2=float(area.max()),
                largest_near_coplanar_area_fraction=float(np.bincount(inverse, weights=area).max()/area.sum()),
                volume_m3=topology(mesh)['oriented_triangle_integral_m3'])


def loft(source, key):
    station_count = 21
    count_ring = 32
    lo, hi = source.bounds[:, 2]
    zs = np.linspace(lo + .02*(hi-lo), hi-.02*(hi-lo), station_count)
    control, evidence, previous = [], [], None
    for z in zs:
        polygon, area, ancestry = section(source, float(z))
        row = resample(polygon, count_ring)
        if previous is not None:
            phase = int(np.argmin([np.square(np.roll(row, k, axis=0)-previous).sum()
                                  for k in range(count_ring)]))
            row = np.roll(row, phase, axis=0)
        rounded = periodic_cubic(row, count_ring)
        control.append(np.c_[rounded, np.full(count_ring, z)])
        evidence.append(dict(z_m=float(z), source_area_m2=area,
                             actual_source_polygon=polygon.tolist(), ancestry=ancestry))
        previous = row
    # Poles use only the measured source extreme-face coordinates. The body
    # remains entirely controlled by the twenty-one actual section contours.
    start = source.vertices[np.isclose(source.vertices[:, 2], lo, rtol=0, atol=1e-11)].mean(0)
    end = source.vertices[np.isclose(source.vertices[:, 2], hi, rtol=0, atol=1e-11)].mean(0)
    controls = np.asarray([np.repeat(start[None], count_ring, axis=0), *control,
                           np.repeat(end[None], count_ring, axis=0)])
    n = len(controls)
    knots = np.r_[np.zeros(4), np.linspace(0, 1, n-2)[1:-1], np.ones(4)]
    spline = BSpline(knots, controls, 3)
    sample_count = 45 if key == 'transverse-vertical-shared' else 37
    parameter = np.linspace(0, 1, sample_count)
    rings = spline(parameter)[1:-1]
    vertices = np.r_[start[None], rings.reshape(-1, 3), end[None]]
    faces = []
    for j in range(count_ring):
        faces.append([0, 1+(j+1) % count_ring, 1+j])
    for i in range(len(rings)-1):
        a, b = 1+i*count_ring, 1+(i+1)*count_ring
        for j in range(count_ring):
            k = (j+1) % count_ring
            faces.extend([[a+j, a+k, b+k], [a+j, b+k, b+j]])
    tip = len(vertices)-1
    base = 1+(len(rings)-1)*count_ring
    for j in range(count_ring):
        faces.append([base+j, base+(j+1) % count_ring, tip])
    faces = np.asarray(faces, dtype=np.int64)
    # The exact proposed runtime Float32 storage is the object being audited.
    vertices = (vertices-ORIGIN).astype('<f4').astype(float)+ORIGIN
    mesh = Mesh(vertices, faces, FRAME)
    if topology(mesh)['oriented_triangle_integral_m3'] < 0:
        mesh = Mesh(vertices, faces[:, [0, 2, 1]], FRAME)
    normals = np.zeros_like(mesh.vertices)
    tri = mesh.triangles
    raw = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    for corner in range(3):
        np.add.at(normals, mesh.faces[:, corner], raw)
    normals /= np.maximum(np.linalg.norm(normals, axis=1)[:, None], 1e-30)
    mesh = Mesh(mesh.vertices, mesh.faces, FRAME,
                dict(normals=np.rint(normals*32767).astype('<i2'),
                     normalSigns=np.ones(len(normals), dtype='i1')),
                dict(derivation='actual-coronal-source-contour-bicubic-loft',
                     exclusive_tissue_ownership=False))
    nearest, _ = BVH(source).closest_points(mesh.vertices, backend='bvh')
    shape = dict(maximum_source_surface_distance_m=float(np.max(nearest['distance_m'])),
                 median_source_surface_distance_m=float(np.median(nearest['distance_m'])),
                 source=quality(source), candidate=quality(mesh),
                 control_stations=evidence, ring_samples=count_ring,
                 long_axis_samples=sample_count, physical_geometry_changed=True)
    return mesh, shape


def construct(output):
    require(output.is_relative_to(WORK) and not output.exists(), 'Use one new local output stage')
    source_audit = read(WORK/'source-v1/source-audit.json')
    bones = [r for r in source_audit['contextParts'] if r['sourcePart'].get('system') == 'skeletal']
    require(len(bones) == 8, 'Eight fixed bone sources required')
    bindings = {str(SOURCE/p): digest_file(SOURCE/p) for p in PARTS.values()}
    bindings.update({str(WHOLE): digest_file(WHOLE),
                     str(ROOT/'reports/ta98-tongue-inventory.json'): digest_file(ROOT/'reports/ta98-tongue-inventory.json')})
    output.mkdir(parents=True)
    atomic_json(output/'plan.json', dict(schema='tongue-coronal-contour-loft-plan/v1',
        source_bindings=bindings, primary_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
        whole_assembly='Male tongue, intrinsic/extrinsic muscles, exposed finite mucosa, aponeurosis and fixed skeleton',
        planned_parts=list(PARTS), unchanged=['superior-longitudinal', 'mucosa', 'aponeurosis', 'all skeleton', 'all extrinsics'],
        interfaces=['All intrinsic pairs remain disjoint macroscopically as in the selected source proposal',
                    'All changed meshes remain inside corrected whole',
                    'No changed mesh crosses exposed mucosa or aponeurosis',
                    'SL aponeurotic 122-facet chart and mucosa/apo 62-facet chart remain literal',
                    'All eight fixed bones remain unchanged and uncrossed'],
        anatomy_scope='Source-supported approximate macroscopic envelopes; no new histological boundary or completion claim',
        interpolation='Periodic convex cubic section controls; positive tensor-product AP B-spline, no halfspace cuts',
        installed=False, anatomical_acceptance=False))
    rows = []
    for key, filename in PARTS.items():
        began = time.monotonic()
        source_path = SOURCE/filename
        source = Mesh.load(source_path)
        require(source.frame == FRAME, 'Source coordinate frame changed')
        if key == 'superior-longitudinal':
            target = output/(key+'.npz')
            shutil.copyfile(source_path, target)
            rows.append(dict(key=key, mesh=str(target), sha256=digest_file(target),
                             geometry_sha256=source.geometry_sha, unchanged_literal=True,
                             source_mesh=str(source_path), source_sha256=digest_file(source_path)))
            continue
        mesh, shape = loft(source, key)
        target = output/(key+'.npz')
        mesh.save(target)
        shape_path = output/(key+'-shape.json')
        atomic_json(shape_path, shape)
        rows.append(dict(key=key, mesh=str(target), sha256=digest_file(target),
                         geometry_sha256=mesh.geometry_sha, source_mesh=str(source_path),
                         source_sha256=digest_file(source_path), shape_evidence=str(shape_path),
                         shape_evidence_sha256=digest_file(shape_path), elapsed_seconds=time.monotonic()-began))
        print(json.dumps(dict(key=key, triangles=len(mesh.faces), elapsed_seconds=time.monotonic()-began)), flush=True)
    atomic_json(output/'construction.json', dict(schema='tongue-coronal-contour-loft/v1',
        meshes=rows, source_bindings=bindings, source_files_modified=False,
        skeleton_modified=False, representation_mode='macroscopic-muscle-reference-envelope',
        tissue_ownership_claim=False, anatomical_acceptance=False, installed=False,
        checks_required=['self', 'intrinsic-pairs', 'coverings', 'bones', 'whole-containment', 'actual-render-review']))


def check_stage(output, stage):
    construction = read(output/'construction.json')
    rows = construction['meshes']
    meshes = {}
    for row in rows:
        require(digest_file(row['mesh']) == row['sha256'], 'Candidate changed before audit')
        meshes[row['key']] = Mesh.load(row['mesh'])
    checks = []
    if stage == 'self':
        for row in rows:
            if row.get('unchanged_literal'):
                continue
            began = time.monotonic()
            result = audit(meshes[row['key']], require_closed=True, expected_components=1, expected_genus=[0])
            checks.append(dict(key=row['key'], seconds=time.monotonic()-began, result=result))
    elif stage == 'intrinsic-pairs':
        for a, b in itertools.combinations(meshes, 2):
            began = time.monotonic()
            result = check_constraint(dict(kind='disjoint', parameters=dict(occupancy='solid')), meshes[a], meshes[b])
            checks.append(dict(first=a, second=b, seconds=time.monotonic()-began, result=result))
    elif stage == 'coverings':
        covering_paths = {'mucosa': SOURCE/'part-00.npz', 'aponeurosis': SOURCE/'part-01.npz'}
        for key, mesh in meshes.items():
            if key == 'superior-longitudinal':
                continue
            for cover, path in covering_paths.items():
                began = time.monotonic()
                result = check_constraint(dict(kind='disjoint', parameters=dict(occupancy='solid')), mesh, Mesh.load(path))
                checks.append(dict(first=key, second=cover, seconds=time.monotonic()-began, result=result))
    elif stage.startswith('bones-'):
        row = next(r for r in rows if r['key'] == stage[len('bones-'):])
        source = read(WORK/'source-v1/source-audit.json')
        for bone in source['contextParts']:
            if bone['sourcePart'].get('system') != 'skeletal':
                continue
            require(digest_file(bone['meshSavePath']) == bone['meshSaveSha256'], 'A fixed bone changed')
            began = time.monotonic()
            mesh = Mesh.load(bone['meshSavePath'])
            params = dict(occupancy='solid') if topology(mesh)['closed_oriented_manifold'] else {}
            result = check_constraint(dict(kind='disjoint', parameters=params), meshes[row['key']], mesh)
            checks.append(dict(first=row['key'], second=bone['partId'], seconds=time.monotonic()-began, result=result))
    elif stage == 'whole-containment':
        whole = Mesh.load(WHOLE)
        for key, mesh in meshes.items():
            began = time.monotonic()
            result = check_constraint(dict(kind='nested', parameters={}), mesh, whole)
            checks.append(dict(first=key, second=str(WHOLE), seconds=time.monotonic()-began, result=result))
    else:
        raise ValueError('Unknown finite check stage')
    passed = all(row['result'].get('passed') for row in checks)
    atomic_json(output/('check-'+stage+'.json'), dict(schema='tongue-coronal-loft-check/v1',
        stage=stage, checks=checks, passed=passed, anatomical_acceptance=False, installed=False,
        construction_sha256=digest_file(output/'construction.json')))
    print(json.dumps(dict(stage=stage, passed=passed, checks=len(checks),
                         elapsed_seconds=sum(row.get('seconds', 0) for row in checks))), flush=True)


def diagnose_existing(output):
    """Document stopped TV reconstruction; retain its literal source as context.

    This is explicitly a partial IL-only diagnostic assembly, not a substitute
    complete TV construction and not a candidate that the installer can use.
    No geometry is regenerated, copied or modified by this stage.
    """
    require(output.is_relative_to(WORK) and (output/'plan.json').is_file()
            and not (output/'construction.json').exists(), 'Use an existing unfinished stage exactly once')
    sections = []
    source = Mesh.load(SOURCE/PARTS['transverse-vertical-shared'])
    lo, hi = source.bounds[:, 2]
    for z in np.linspace(lo+.02*(hi-lo), hi-.02*(hi-lo), 21):
        loops, ancestry = section(source, float(z), allow_multiple=True)
        areas = []
        for p in loops:
            nxt = np.roll(p, -1, axis=0)
            areas.append(abs(float((p[:, 0]*nxt[:, 1]-p[:, 1]*nxt[:, 0]).sum()/2)))
        sections.append(dict(z_m=float(z), loop_count=len(loops), loop_areas_m2=areas,
                             polygons=[p.tolist() for p in loops], ancestry=ancestry))
    atomic_json(output/'tv-source-topology-blocker.json', dict(
        schema='tongue-loft-source-section-topology/v1', source_mesh=str(SOURCE/PARTS['transverse-vertical-shared']),
        source_sha256=digest_file(SOURCE/PARTS['transverse-vertical-shared']), sections=sections,
        multiple_loop_stations=[r['z_m'] for r in sections if r['loop_count'] != 1],
        reason='A single source-supported loft cannot bridge or discard separate source section regions without a new anatomical correspondence decision',
        geometry_reconstructed=False, installed=False))
    rows = []
    for key, filename in PARTS.items():
        source_path = SOURCE/filename
        target = output/(key+'.npz')
        literal = key in ('superior-longitudinal', 'transverse-vertical-shared')
        if key == 'transverse-vertical-shared':
            target = source_path
        require(target.is_file(), 'A previously constructed diagnostic mesh is absent')
        mesh = Mesh.load(target)
        row = dict(key=key, mesh=str(target), sha256=digest_file(target), geometry_sha256=mesh.geometry_sha,
                   source_mesh=str(source_path), source_sha256=digest_file(source_path), unchanged_literal=literal)
        if literal:
            require(digest_file(target) == digest_file(source_path), 'Retained diagnostic source is not literal')
        else:
            shape = output/(key+'-shape.json')
            require(shape.is_file(), 'Actual IL shape evidence is absent')
            row.update(shape_evidence=str(shape), shape_evidence_sha256=digest_file(shape))
        rows.append(row)
    atomic_json(output/'construction.json', dict(schema='tongue-coronal-contour-loft/v1', meshes=rows,
        source_bindings=read(output/'plan.json')['source_bindings'],
        partial_diagnostic_only=True, tv_reconstruction_blocked=True, tv_source_retained_for_context_only=True,
        no_candidate_packing_authorized=True, source_files_modified=False, skeleton_modified=False,
        anatomical_acceptance=False, installed=False, actual_completed_lofts=2,
        representation_mode='macroscopic-muscle-reference-envelope', tissue_ownership_claim=False))
    print(json.dumps(dict(partial_diagnostic_only=True, actual_completed_lofts=2,
                         multiple_loop_stations=[r['z_m'] for r in sections if r['loop_count'] != 1])), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--stage', default='construct')
    args = parser.parse_args()
    if args.stage == 'construct':
        construct(args.output.resolve())
    elif args.stage == 'diagnose-existing':
        diagnose_existing(args.output.resolve())
    else:
        check_stage(args.output.resolve(), args.stage)
