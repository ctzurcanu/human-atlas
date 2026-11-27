"""Bounded trust-region ACA fit with separate signed-tissue residuals.

This writes a construction input only after a completed solve. It never edits
an atlas, relaxes vessel calibre, or accepts its own geometry for installation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np
import open3d as o3d
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_sections import (
    arclength, path_tangents, redirect_local_sections, bend_diagnostics,
)

WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1'
CAT = ROOT / '.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
ORIGIN = np.array([0., 1.725, .028])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, default=lambda x:
        x.tolist() if hasattr(x, 'tolist') else str(x)) + '\n')


def scene(mesh):
    result = o3d.t.geometry.RaycastingScene(nthreads=1)
    result.add_triangles(
        o3d.core.Tensor(((mesh.vertices - ORIGIN) * 1000).astype(np.float32)),
        o3d.core.Tensor(mesh.faces.astype(np.uint32)),
    )
    return result


def run(side, output, seed=None, joint_root=False, analytic=False):
    started = time.monotonic()
    output.mkdir(exist_ok=False, parents=True)
    info = json.loads((WORK / 'source-inputs.json').read_text())
    assert sha(CAT) == info['source_sha256']
    if joint_root:
        inventory = json.loads((WORK / 'native-carotid-terminal-cohort-v7b/complete-authority-inventory.json').read_text())
        assert inventory['term_count'] == 129 and inventory['source_catalogue_sha256'] == sha(CAT)
        assert seed is not None, 'The earlier junction needs the measured nonvascular graph seed'
    sign = 1 if side == 'left' else -1
    source_id = f'ZA:Anterior cerebral artery.{"l" if sign == 1 else "r"}'
    row = next(r for r in info['sources'] if r['part_id'] == source_id)
    source = Mesh.load(WORK / row['mesh_file'])
    original = np.load(WORK / 'superior-course-candidate-v2-clamped/source-and-candidate-fields.npz')
    fields = original[side + '_source_fields']
    rings = original[side + '_rings']
    centres = fields[0].mean(1)
    tangents = path_tangents(centres)
    target = np.array([sign * .0035, 1.723, .0285])
    arc = arclength(centres)
    # Four free spline knots give a smooth centre course with twelve degrees
    # of freedom. Fixed endpoints and communication point are intrinsic.
    start = 2 if joint_root else 6
    knots = np.array([2, 5, 8, 10, 13, 16, 20, 23]) if joint_root else np.array([6, 9, 11, 13, 16, 20, 23])
    free = np.array([1, 2, 3, 5, 6]) if joint_root else np.array([1, 2, 4, 5])
    dimensions = len(free) * 3
    positions = np.array([centres[start], target, centres[23]])
    lengths = arclength(positions)
    initial_curve = CubicSpline(lengths, positions, axis=0,
        bc_type=((1, tangents[start]), (1, tangents[23])))
    parameter = np.r_[
        lengths[1] * (arc[start:14] - arc[start]) / (arc[13] - arc[start]),
        lengths[1] + (lengths[2] - lengths[1]) * (arc[14:24] - arc[13]) / (arc[23] - arc[13]),
    ]
    initial = centres.copy()
    initial[start:24] = initial_curve(parameter)
    initial[[start, 13, 23]] = positions
    if seed is not None:
        proof = json.loads((seed.parent / 'route-stage.json').read_text())
        assert proof['complete'] and proof['passed_all_four_routes']
        initial = np.load(seed)[side + '_centres'].copy()
        assert np.array_equal(initial[[start, 23]], centres[[start, 23]])
        assert np.allclose(initial[13], target, atol=1e-14, rtol=0)
        parameter = arclength(initial[start:24])
    base = initial[knots].copy()
    frozen = np.r_[np.arange(start + 1), np.arange(23, 37)]
    dense_parameter = np.linspace(parameter[0], parameter[-1], 257)
    source_radius = np.linalg.norm(fields[0] - centres[:, None], axis=2).max(1)
    dense_radius = np.interp(dense_parameter, parameter, source_radius[start:24])

    contract = json.loads((WORK / 'active-fixed-context-v1/context-contract.json').read_text())
    assert contract['complete'] and contract['master_sha256'] == info['source_sha256']
    manifest = json.loads((WORK / 'closed-sign-tissues-v4/closed-tissue-manifest.json').read_text())
    assert manifest['complete']
    signed = []
    for entry in manifest['rows']:
        if not entry['closed']:
            continue
        m = Mesh.load(entry['path'])
        assert m.geometry_sha == entry['geometry_sha256'] and entry['Float32_rebased_closed_oriented']
        signed.append((entry['id'], scene(m), m.bounds))
    # Do not include those same closed tissues again as unsigned barriers.
    # Inside a nerve the unsigned gradient points farther inward; each signed
    # tissue gets its own residual rows, without hard-minimum tie switching.
    vv, ff, offset = [], [], 0
    for entry in contract['neighbors']:
        if any(word in entry['id'] for word in ['Optic nerve', 'Optic chiasm', 'Olfactory nerve', 'Straight gyrus']):
            continue
        m = Mesh.load(entry['path'])
        assert sha(entry['path']) == entry['sha256']
        if joint_root and entry['id'] == f'ZA:Middle cerebral artery (M1-segment).{"l" if sign == 1 else "r"}':
            # A1 enters the intended terminal arterial junction here. Only
            # this declared root patch is omitted from seed/fit penalties;
            # acceptance requires its actual union, never an ignored clash.
            table = np.load(WORK / f'native-carotid-terminal-cohort-v7b/mca-{"l" if sign == 1 else "r"}-stations.npz')
            all_mca = Mesh.load(WORK / f'native-carotid-terminal-cohort-v7b/mca-{"l" if sign == 1 else "r"}.npz')
            mca_rings = table['rings']
            mca_station = np.zeros(len(all_mca.vertices), dtype=int)
            mca_station[mca_rings] = np.arange(len(mca_rings))[:, None]
            whole_stations = mca_station[all_mca.faces].min(1)
            keep = whole_stations[np.asarray(entry['source_face_ids'])] >= 18
            m = Mesh(m.vertices, m.faces[keep], m.frame)
            if not len(m.faces):
                continue
        vv.append(m.vertices)
        ff.append(m.faces + offset)
        offset += len(m.vertices)
    other = scene(Mesh(np.concatenate(vv), np.concatenate(ff), source.frame))
    station = np.zeros(len(source.vertices), dtype=int)
    station[rings] = np.arange(37)[:, None]
    faces = source.faces[(station[source.faces].min(1) >= start) & (station[source.faces].max(1) <= 23)]
    edges = np.unique(np.sort(np.r_[faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]], axis=1), axis=0)
    used = np.unique(faces)
    cache, history = {}, []
    evaluations = 0

    def evaluate(x):
        nonlocal evaluations
        if time.monotonic() - started > 32:
            raise TimeoutError('32-second stage deadline; no completed construction candidate')
        if 'x' in cache and np.array_equal(cache['x'], x):
            return cache
        p = base.copy()
        p[free] += x.reshape(-1, 3) / 1000
        curve = CubicSpline(parameter[knots - start], p, axis=0,
            bc_type=((1, tangents[start]), (1, tangents[23])))
        c = centres.copy()
        c[start:24] = curve(parameter)
        c[[start, 13, 23]] = positions
        f, rot = redirect_local_sections(fields, c, dtype=np.float64)
        if analytic:
            # A clamped spline's exact tangent must drive the section frame.
            # Central differences missed tight sub-station bends and produced
            # a plane jump when the fixed endpoint profile was restored.
            target_tangent = curve(parameter, 1)
            target_tangent /= np.linalg.norm(target_tangent, axis=1)[:, None]
            rot = np.broadcast_to(np.eye(3), (37, 3, 3)).copy()
            for station_index, target_direction in zip(range(start, 24), target_tangent):
                a, b = tangents[station_index], target_direction
                v = np.cross(a, b)
                cosine = float(np.dot(a, b))
                assert cosine > -.999999, 'Antiparallel section redirection is not a valid local fit'
                skew = np.array([[0., -v[2], v[1]], [v[2], 0., -v[0]], [-v[1], v[0], 0.]])
                rot[station_index] += skew + skew @ skew / (1 + cosine)
            f = np.einsum('sij,ksvj->ksvi', rot, fields - centres[None, :, None]) + c[None, :, None]
        f[:, frozen] = fields[:, frozen]
        rot[frozen] = np.eye(3)
        vertices = source.vertices.copy()
        vertices[rings] = f[0]
        samples = np.r_[vertices[used], vertices[edges].mean(1), vertices[faces].mean(1)]
        tensor = o3d.core.Tensor(((samples - ORIGIN) * 1000).astype(np.float32))
        distances = [other.compute_distance(tensor, nthreads=1).numpy().astype(float)]
        for _, sc, box in signed:
            d = np.full(len(samples), 1., dtype=float)
            near = np.all(samples >= box[0] - .0003, axis=1) & np.all(samples <= box[1] + .0003, axis=1)
            if near.any():
                d[near] = sc.compute_signed_distance(
                    o3d.core.Tensor(((samples[near] - ORIGIN) * 1000).astype(np.float32)),
                    nthreads=1, nsamples=3).numpy().astype(float)
            distances.append(d)
        distances = np.array(distances)
        if analytic:
            first, second = curve(dense_parameter, 1), curve(dense_parameter, 2)
            speed = np.linalg.norm(first, axis=1)
            support = np.linalg.norm(np.cross(first, second), axis=1) / np.maximum(speed, 1e-12) ** 3 * dense_radius
        else:
            support = np.asarray(bend_diagnostics(c, f)['curvature_support'])[start:24]
        t = path_tangents(c)
        residual = np.r_[
            (40 * np.maximum(.25 - distances, 0)).ravel(),
            (200 if analytic else 20) * np.maximum(support - (.70 if analytic else .80), 0),
            .12 * np.diff(t[start - 1:25], axis=0).ravel(),
            .01 * x,
        ]
        evaluations += 1
        cache.clear()
        cache.update(x=x.copy(), residual=residual, centres=c, fields=f,
            rotations=rot, distances=distances, support=support, vertices=vertices)
        return cache

    def fun(x):
        return evaluate(x)['residual']

    def jac(x):
        # Absolute 1-micrometre differences are above local Float32 query
        # precision. Relative steps at zero-valued displacement knots are not.
        first = fun(x).copy()
        result = np.empty((len(first), len(x)))
        for k in range(len(x)):
            shifted = x.copy()
            shifted[k] += .001
            result[:, k] = (fun(shifted) - first) / .001
        return result

    initial_check = evaluate(np.zeros(dimensions))
    initial_min = float(initial_check['distances'].min())
    try:
        result = least_squares(fun, np.zeros(dimensions), jac=jac,
            bounds=(-5 * np.ones(dimensions), 5 * np.ones(dimensions)), method='trf',
            max_nfev=100, ftol=1e-7, xtol=1e-7, gtol=1e-7)
        value = evaluate(result.x)
        completed = True
        message = result.message
        physical = bool(value['distances'].min() >= .2 and value['support'].max() <= .85)
        # The construction-stage objective is feasibility. A finite completed
        # evaluation budget can supply a fully specified feasible course even
        # when the optional energy minimum has not converged. Its actual mesh
        # and context must still pass independent acceptance gates.
        passed = bool(physical and (analytic or result.success))
    except TimeoutError as exc:
        result, value, completed, passed, message = None, cache, False, False, str(exc)
    np.savez_compressed(output / 'diagnostic-checkpoint.npz',
        **{k: value[k] for k in ['centres', 'fields', 'rotations', 'distances', 'support', 'x']})
    report = dict(schema='aca-signed-soft-trust-region/v9' if analytic else 'aca-signed-soft-trust-region/v6' if joint_root else 'aca-signed-soft-trust-region/v5', source_catalogue_sha256=sha(CAT),
        side=side, completed=completed, passed=passed, installed=False,
        evaluations=evaluations, elapsed_seconds=time.monotonic() - started,
        solver_message=str(message), initial_signed_clearance_mm=initial_min,
        optimizer_energy_converged=bool(result.success) if result is not None else False,
        completed_geometry_constraint_stage=bool(completed and passed),
        curvature_method='257 analytic cubic derivative samples with interpolated held section support; independent dense and actual mesh proof required' if analytic else 'source-station central differences',
        minimum_sampled_clearance_mm=float(value['distances'].min()),
        maximum_middle_curvature_support=float(value['support'].max()),
        free_spline_knots=knots[free], fixed_source_profiles=frozen,
        measured_seed_sha256=sha(seed) if seed is not None else None,
        intended_terminal_root_union_required=joint_root,
        solver_script_sha256=sha(Path(__file__)),
        target_profile13_m=target, declared_clearance_mm=.2,
        method='Separate signed closed-tissue squared-hinge residuals, unsigned other surfaces, trust-region least-squares; source-calibre sections retained intrinsically. Actual Float64 triangle and rendered acceptance remains required.')
    save(output / 'solve-proof.json', report)
    if passed:
        attrs = {k: a.copy() for k, a in source.attributes.items()}
        words = source.attributes['normals'][rings].astype(float)
        words = np.einsum('sij,svj->svi', value['rotations'], words)
        words /= np.linalg.norm(words, axis=2)[..., None]
        words = np.rint(words * 32767).astype(np.int16)
        words[frozen] = source.attributes['normals'][rings[frozen]]
        attrs['normals'][rings] = words
        mesh = Mesh(value['vertices'], source.faces, source.frame, attrs)
        mesh.save(output / 'complete-course-input.npz')
        np.savez_compressed(output / 'complete-fields.npz',
            centres=value['centres'], fields=value['fields'], rotations=value['rotations'],
            rings=rings, source_fields=fields, fixed_profiles=frozen)
    print(json.dumps(report, default=lambda x: x.tolist() if hasattr(x, 'tolist') else str(x)))
    return passed


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--side', choices=['left', 'right'], required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--seed', type=Path)
    parser.add_argument('--joint-root', action='store_true')
    parser.add_argument('--analytic', action='store_true')
    args = parser.parse_args()
    sys.exit(0 if run(args.side, args.output.resolve(), args.seed, args.joint_root, args.analytic) else 2)
