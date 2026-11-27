"""Bounded common bladder-neck field, fitted against unchanged pelvic sources.

Native signed distances only guide a candidate. Original Float64 triangles,
finite-wall topology, sections and conduit interfaces remain the final gates.
"""
from pathlib import Path
import json

import numpy as np
from scipy.optimize import least_squares

from .audit import audit, topology
from .common import atomic_json, digest_file
from .db import Database
from .harmonic_wall import tetrahedralize
from .model import Mesh
from .open3d_backend import SurfaceIndex, library
from .registry import load_registered


def smooth(t):
    t = np.clip(t, 0., 1.)
    return t * t * (3. - 2. * t)


def field(points, anchor, parameters):
    """One position field carries the complete wall, never separate coats."""
    p = np.asarray(points, float); a = np.asarray(anchor, float)
    lower_scale, height_mm, lift_mm, forward_mm = parameters
    weight = 1. - smooth((p[:, 1] - a[1]) / (height_mm / 1000.))
    radius = np.linalg.norm(p[:, [0, 2]] - a[[0, 2]], axis=1)
    broad = smooth((radius - .003) / .011)
    q = p.copy()
    q[:, 1] += lift_mm / 1000. * weight * broad
    # Compose monotone maps instead of combining their displacements. The
    # axial derivative is >= 1 - 1.5*lift/height; then radial scaling depends
    # only on the new height. The last z derivative is >= 1-1.5*forward/11.
    weight = 1. - smooth((q[:, 1] - a[1]) / (height_mm / 1000.))
    scale = 1. - (1. - lower_scale) * weight
    q[:, [0, 2]] = a[[0, 2]] + (p[:, [0, 2]] - a[[0, 2]]) * scale[:, None]
    radius = np.linalg.norm(q[:, [0, 2]] - a[[0, 2]], axis=1)
    q[:, 2] += forward_mm / 1000. * weight * smooth((radius - .003) / .011)
    return q


def fit(source, scan, output):
    source, scan, output = Path(source).resolve(), Path(scan).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    output.mkdir(parents=True)
    mesh = Mesh.load(source / 'whole-wall.npz')
    initial = json.loads((source / 'fit.json').read_text())
    anchor = np.asarray(initial['target_anchor_m'])
    db = Database(); db.assert_current_catalogues('local-ta98')
    scan_rows = json.loads(scan.read_text())['rows']
    guides = []; excluded = []; bindings = []
    for row in scan_rows:
        s = row['result']
        if not s['crossing_face_pairs'] and not s['coplanar_overlap_face_pairs']:
            continue
        if row['part_id'].startswith('RENALDETAIL:'):
            excluded.append(dict(part_id=row['part_id'], reason='Explicit adjoining ureter; distal course and ostial join require coupled construction'))
            continue
        record, other = load_registered(db, row['part_id'], 'local-ta98')
        top = topology(other)
        if not top['closed_oriented_manifold']:
            excluded.append(dict(part_id=row['part_id'], reason='Open/nonmanifold guide: native occupancy is not meaningful; retained for final triangle checks'))
            continue
        guides.append((row['part_id'], SurfaceIndex(other)))
        bindings.append(dict(part_id=row['part_id'], geometry_sha256=other.geometry_sha))
    if not guides:
        raise ValueError('No verified closed pelvic fitting guides')
    p = np.r_[mesh.vertices, mesh.triangles.mean(1)]
    p0, tets, _, _, _ = tetrahedralize(mesh.vertices, mesh.faces)
    def determinant(points):
        t = points[tets]
        return np.linalg.det(t[:, 1:] - t[:, :1])
    original_determinants = determinant(p0)
    library()  # Native import once, before timing the bounded solve.
    evaluations = []
    def residual(parameters):
        q = field(p, anchor, parameters); errors = []
        for pid, index in guides:
            tensor = library().core.Tensor(((q - index.origin) / index.scale).astype(np.float32))
            d = index.scene.compute_signed_distance(tensor, nthreads=2, nsamples=3).numpy().astype(float) * index.scale
            # The immediate urethral/prostatic recipient needs its own shared
            # lumen interface; do not force an empty gap there to hide it.
            if pid == 'ZA:Prostate':
                port = np.linalg.norm(q[:, [0, 2]] - anchor[[0, 2]], axis=1) < .003
                d[port & (q[:, 1] < anchor[1] + .003)] = .001
            errors.extend(np.maximum(0., .00025 - d) * 1000.)
        # Retain the measured source envelope wherever packing permits it.
        errors.extend((np.asarray(parameters) - [1., 26., 0., 0.]) * [.025, .0005, .002, .002])
        ratios = determinant(field(p0, anchor, parameters)) / original_determinants
        errors.extend(10. * np.maximum(0., .15 - ratios))
        result = np.asarray(errors)
        evaluations.append(dict(parameters=np.asarray(parameters).tolist(), squared_residual=float(result @ result)))
        return result
    result = least_squares(residual, [.5, 30., 10., 3.],
                           bounds=([.27, 27., 0., 0.], [.8, 40., 16., 7.]),
                           diff_step=.002, max_nfev=55, ftol=1e-5, xtol=1e-6, gtol=1e-5)
    q = field(mesh.vertices, anchor, result.x)
    fitted = Mesh(q, mesh.faces, mesh.frame)
    check = audit(fitted, require_closed=True, expected_components=1,
                  expected_genus=[2], max_pairs=1000000)
    # Check the same piecewise-linear field throughout a constrained wall,
    # rather than trusting a smooth-looking boundary alone.
    p1 = field(p0, anchor, result.x)
    ratios = determinant(p1) / original_determinants
    fitted.save(output / 'whole-wall.npz')
    evidence = dict(schema='bladder-pelvic-neck-field/v1', source=str(source),
        source_sha256=mesh.geometry_sha, scan_sha256=digest_file(scan),
        fixed_pelvic_sources=bindings, excluded_from_native_occupancy_loss=excluded,
        parameters=result.x.tolist(), parameter_names=['inferior_radial_scale', 'transition_height_mm', 'broad_base_lift_mm', 'broad_base_forward_mm'],
        bounds=[[.27, 27., 0., 0.], [.8, 40., 16., 7.]],
        anchor_m=anchor.tolist(), evaluations=evaluations, solver_converged=bool(result.success),
        sample_points=len(p), exact_final_triangles_checked=check,
        wall_tetrahedra=len(tets), minimum_tetrahedron_volume_ratio=float(ratios.min()),
        maximum_displacement_m=float(np.linalg.norm(q - mesh.vertices, axis=1).max()),
        passed=check['passed'] and bool(np.all(ratios > 0)),
        source_faces_literal=True, skeleton_modified=False, installed=False,
        all_adjoining_interfaces_verified=False, anatomical_acceptance=False)
    atomic_json(output / 'fit.json', evidence)
    return evidence


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('source'); p.add_argument('--scan', required=True); p.add_argument('--output', required=True)
    args = p.parse_args(); r = fit(args.source, args.scan, args.output)
    print(json.dumps({k: r[k] for k in ['passed', 'parameters', 'minimum_tetrahedron_volume_ratio', 'maximum_displacement_m']}, indent=2))
