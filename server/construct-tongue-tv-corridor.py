"""One explicitly bounded original central TV reference-core construction.

This deliberately does not reconstruct the rejected Boolean whole-minus-IL
residual. Coronal whole, exposed mucosal, SL and paired IL profiles determine
the admissible central corridor. Maximal inscribed section ellipses describe
an approximate interwoven TV territory, not an exclusive histological volume.
No source, covering, SL, skeleton or installed catalogue is modified.
"""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT/'server')]
import numpy as np
from scipy.interpolate import BSpline, PchipInterpolator
from scipy.optimize import brentq
from geometry_tools.audit import topology
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.model import Mesh

spec = importlib.util.spec_from_file_location('tongue_coronal_loft', ROOT/'server/construct-tongue-coronal-lofts.py')
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)
WORK, FRAME, SOURCE, ORIGIN = helper.WORK, helper.FRAME, helper.SOURCE, helper.ORIGIN
OLD = WORK/'coronal-contour-loft-v2'
CLEARANCE = .0002


def require(ok, message):
    if not ok:
        raise ValueError(message)


def loops_at(mesh, z):
    if not mesh.bounds[0, 2] < z < mesh.bounds[1, 2]:
        return []
    loops, _ = helper.section(mesh, float(z), allow_multiple=True)
    return loops


def intervals(loops, value, axis):
    """Actual polygon intersections, even/odd occupied intervals per loop."""
    out = []
    for polygon in loops:
        values = []
        other = 1-axis
        for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
            if (a[axis] <= value < b[axis]) or (b[axis] <= value < a[axis]):
                u = (value-a[axis])/(b[axis]-a[axis])
                values.append(float(a[other]+u*(b[other]-a[other])))
        values.sort()
        require(len(values) % 2 == 0, 'A source section has ambiguous occupied intervals')
        out.extend(zip(values[::2], values[1::2]))
    return out


def patch_points(mesh, z, faces=348):
    """Cut only the literal exposed outer carrier, excluding finite rim/inner."""
    tri = mesh.triangles[:faces]
    distance = tri[:, :, 2]-z
    selected = np.flatnonzero((distance.min(1) < 0) & (distance.max(1) > 0))
    hits = []
    for f in selected:
        for j in range(3):
            k = (j+1) % 3
            if distance[f, j]*distance[f, k] >= 0:
                continue
            u = distance[f, j]/(distance[f, j]-distance[f, k])
            hits.append(tri[f, j]+u*(tri[f, k]-tri[f, j]))
    require(hits, 'Actual exposed mucosal body profile is absent')
    return np.asarray(hits)


def chart_il(path, output):
    """Separate only incompatible corner fields; every physical triangle literal."""
    source = Mesh.load(path)
    tri = source.triangles
    raw = np.cross(tri[:, 1]-tri[:, 0], tri[:, 2]-tri[:, 0])
    face_normals = raw/np.linalg.norm(raw, axis=1)[:, None]
    words = source.attributes['normals'][source.faces].copy()
    signs = source.attributes['normalSigns'][source.faces]
    effective = words.astype(float)*signs[:, :, None]
    bad = np.einsum('fci,fi->fc', effective, raw) <= 0
    effective_words = np.rint(effective).astype('<i2')
    replacement = np.rint(face_normals*32767).astype('<i2')
    effective_words[bad] = np.broadcast_to(replacement[:, None], effective_words.shape)[bad]
    keys = np.c_[source.faces.reshape(-1), effective_words.reshape(-1, 3)]
    charts, inverse = np.unique(keys, axis=0, return_inverse=True)
    mesh = Mesh(source.vertices[charts[:, 0]], inverse.reshape(-1, 3), FRAME,
                dict(normals=charts[:, 1:].astype('<i2'), normalSigns=np.ones(len(charts), dtype='i1')),
                dict(source_mesh=str(path), source_sha256=digest_file(path),
                     physical_geometry_changed=False, incompatible_corner_charts_only=True))
    require(mesh.geometry_sha == source.geometry_sha and np.array_equal(mesh.triangles, tri),
            'Corner chart changed physical IL geometry')
    require(np.all(np.einsum('fci,fi->fc', mesh.attributes['normals'][mesh.faces].astype(float), raw) > 0),
            'A corrected corner field points through its physical triangle')
    mesh.save(output)
    return dict(key=output.stem, mesh=str(output), sha256=digest_file(output),
                geometry_sha256=mesh.geometry_sha, source_mesh=str(path), source_sha256=digest_file(path),
                geometry_unchanged=True, normal_chart_changed=True, repaired_corners=int(bad.sum()),
                indexed_vertices=len(mesh.vertices), triangles=len(mesh.faces))


def bound_section(z, whole, sl, il, muco):
    whole_loops, sl_loops = loops_at(whole, z), loops_at(sl, z)
    require(whole_loops and sl_loops, 'Whole or dorsal SL profile absent')
    whole_points = np.concatenate(whole_loops)
    cx = float((whole_points[:, 0].min()+whole_points[:, 0].max())/2)
    central_whole = [r for r in intervals(whole_loops, cx, 0) if r[1] > r[0]]
    require(central_whole, 'The whole has no central section interval')
    # Highest connected body interval excludes a separate computational basal
    # stem where one exists. The exposed carrier supplies its actual lower edge.
    body_bottom, body_top = max(central_whole, key=lambda r:r[1])
    exposed = patch_points(muco, z)
    mucosal_bottom = float(exposed[:, 1].min())
    # 0.42 mm is the selected coat's maximum modeled thickness; the 0.20 mm
    # separation is an explicit display-envelope allowance, not real histology.
    lower = max(body_bottom+.00062, mucosal_bottom+.00062)
    il_profiles = []
    for name, mesh in il.items():
        loops = loops_at(mesh, z)
        if loops:
            points = np.concatenate(loops)
            upper = float(points[:, 1].max())
            lower = max(lower, upper+CLEARANCE)
            il_profiles.append(dict(key=name, superior_y_m=upper, loops=[p.tolist() for p in loops]))
    dorsal = intervals(sl_loops, cx, 0)
    require(dorsal, 'SL underside has no central intersection')
    upper = min(body_top-.00062, min(r[0] for r in dorsal)-CLEARANCE)
    row = dict(z_m=float(z), center_x_m=cx, mucosal_lower_y_m=mucosal_bottom,
               whole_body_interval_m=[body_bottom, body_top], lower_y_m=float(lower),
               upper_y_m=float(upper), height_m=float(upper-lower),
               whole_loops=[p.tolist() for p in whole_loops],
               sl_loops=[p.tolist() for p in sl_loops], il_profiles=il_profiles,
               corridor_positive=bool(upper > lower))
    if upper <= lower:
        return None, row
    cy, ry = (upper+lower)/2, (upper-lower)/2
    rx_limits = []
    # The maximal admissible radius is computed once from source intervals at
    # every ring angle. It is not reduced and retried against later audits.
    angles = np.arange(128)*2*np.pi/128
    for angle in angles:
        c, s = float(np.cos(angle)), float(np.sin(angle))
        y = cy+ry*s
        rooms = [r for r in intervals(whole_loops, y, 1) if r[0] < cx < r[1]]
        if not rooms:
            row.update(corridor_positive=False, failure='Whole central horizontal interval absent')
            return None, row
        left, right = max(rooms, key=lambda r:r[1]-r[0])
        left, right = left+.00062, right-.00062
        for a, b in intervals(sl_loops, y, 1):
            if a <= cx <= b:
                row.update(corridor_positive=False, failure='SL occupies central vertical profile')
                return None, row
            if b < cx:
                left = max(left, b+CLEARANCE)
            elif a > cx:
                right = min(right, a-CLEARANCE)
        if abs(c) > 1e-8:
            rx_limits.append((right-cx)/c if c > 0 else (left-cx)/c)
    rx = float(min(rx_limits))
    if not np.isfinite(rx) or rx <= 0:
        row.update(corridor_positive=False, failure='No positive lateral width inside the measured bounds')
        return None, row
    row.update(center_y_m=float(cy), radius_x_m=rx, radius_y_m=float(ry),
               area_m2=float(np.pi*rx*ry), maximal_inscribed_radius_computed_once=True)
    return np.array([cx, cy, z, rx, ry]), row


def main(output, posterior_termination=None, refine_bound=None):
    began = time.monotonic()
    require(output.is_relative_to(WORK) and not output.exists(), 'Use a new explicit local stage')
    output.mkdir(parents=True)
    whole = Mesh.load(helper.WHOLE)
    sl = Mesh.load(SOURCE/'part-02.npz')
    muco = Mesh.load(SOURCE/'part-00.npz')
    old_tv = Mesh.load(SOURCE/'part-05.npz')
    inputs = [helper.WHOLE, SOURCE/'part-00.npz', SOURCE/'part-01.npz', SOURCE/'part-02.npz', SOURCE/'part-05.npz',
              OLD/'inferior-longitudinal-left.npz', OLD/'inferior-longitudinal-right.npz',
              ROOT/'reports/ta98-tongue-inventory.json']
    bindings = {str(p):digest_file(p) for p in inputs}
    if posterior_termination is not None:
        require(posterior_termination.is_relative_to(WORK) and posterior_termination.is_file(),
                'Bind the actual previously stopped corridor report')
        bindings[str(posterior_termination)] = digest_file(posterior_termination)
    refinement = None
    if refine_bound is not None:
        require(refine_bound.is_relative_to(WORK) and refine_bound.is_file(), 'Bind actual predecessor crossing evidence')
        refinement = json.loads(refine_bound.read_text())
        require(refinement.get('schema') == 'tongue-tv-SL-candidate-crossings/v1'
                and refinement['summary'].get('complete_within_pair_budget') is True
                and refinement['summary'].get('crossing_face_pairs') == 34
                and refinement['source_SL_sha256'] == digest_file(SOURCE/'part-02.npz')
                and digest_file(refinement['candidate_TV']) == refinement['candidate_TV_sha256']
                and digest_file(refinement['witnesses']) == refinement['witnesses_sha256'],
                'Numerical refinement requires the actual unchanged predecessor witnesses')
        bindings[str(refine_bound)] = digest_file(refine_bound)
        bindings[str(Path(refinement['candidate_TV']).resolve())] = refinement['candidate_TV_sha256']
        bindings[str(Path(refinement['witnesses']).resolve())] = refinement['witnesses_sha256']
    atomic_json(output/'plan.json', dict(schema='tongue-tv-core-corridor-plan/v1', input_bindings=bindings,
        primary_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
        anatomical_extent='Original representative central interwoven transverse/vertical territory below the dorsal SL sheet and above the inferior longitudinal envelopes',
        rejected_assignment='Prior Boolean whole-minus-IL remainder and its automatically disconnected section wings are not anatomical authority',
        tissue_ownership_claim=False, exact_muscle_extent_claim=False,
        exterior='Corrected whole and actual exposed mucosal lower/reflection boundary',
        interfaces=['Fixed skeleton', 'Literal SL and its 122 shared aponeurotic facets', 'Literal mucosa/aponeurosis 62 facets', 'New paired IL superior envelopes'],
        source_section_clearance_m=CLEARANCE, maximum_modeled_mucosal_thickness_m=.00042,
        no_automatic_narrowing_or_parameter_search=True,
        posterior_termination_decision=(None if posterior_termination is None else
            'Explicitly authorized macroscopic posterior termination at the bracketed zero of the unchanged corridor-height function; not a measured anatomical muscle border'),
        numerical_bound_refinement=(None if refinement is None else
            'Physical-Z shape-preserving interpolation through the same actual-facet section bounds, with exact crossed-facet vertex Z events and predecessor witness ring stations added; no anatomical parameter changes'),
        installed=False, anatomical_acceptance=False))
    il_rows = [chart_il(OLD/(key+'.npz'), output/(key+'.npz'))
               for key in ('inferior-longitudinal-left', 'inferior-longitudinal-right')]
    il = {r['key']:Mesh.load(r['mesh']) for r in il_rows}
    lo, hi = old_tv.bounds[:, 2]
    stations = np.linspace(lo+.02*(hi-lo), hi-.02*(hi-lo), 25)
    controls, sections = [], []
    for z in stations:
        value, evidence = bound_section(float(z), whole, sl, il, muco)
        sections.append(evidence)
        controls.append(value)
    failed = [i for i,r in enumerate(sections) if not r['corridor_positive']]
    posterior_end = None
    if failed and posterior_termination is not None:
        prior = json.loads(posterior_termination.read_text())
        require(prior.get('schema') == 'tongue-tv-actual-coronal-bounds/v1'
                and prior.get('failed_stations') == failed == [0]
                and prior['sections'] == sections,
                'The accepted root decision must bind exactly the same unmodified corridor')
        bracket = [sections[0]['z_m'], sections[1]['z_m']]
        evaluations = []
        def height(z):
            _, row = bound_section(float(z), whole, sl, il, muco)
            evaluations.append(dict(z_m=float(z), height_m=row['height_m']))
            return row['height_m']
        root = brentq(height, *bracket, xtol=1e-11, rtol=1e-12, maxiter=40)
        _, boundary = bound_section(float(root), whole, sl, il, muco)
        require(abs(boundary['height_m']) < 1e-8, 'Bracketed posterior closure did not converge')
        posterior_end = np.array([boundary['center_x_m'],
            (boundary['lower_y_m']+boundary['upper_y_m'])/2, root, 0., 0.])
        atomic_json(output/'posterior-termination.json', dict(
            schema='tongue-tv-posterior-zero-corridor/v1', source_report=str(posterior_termination),
            source_report_sha256=digest_file(posterior_termination), bracket_m=bracket,
            root_z_m=float(root), absolute_height_residual_m=abs(boundary['height_m']), evaluations=evaluations,
            boundary=boundary, endpoint=posterior_end.tolist(),
            all_thickness_clearance_ellipse_parameters_unchanged=True,
            all_24_positive_source_sections_retained=True,
            qualification='Geometric approximation of representative posterior TV macrocore termination, not patient-measured anatomical border',
            primary_reference='https://pmc.ncbi.nlm.nih.gov/articles/PMC3687025/',
            anatomical_acceptance=False))
        controls = controls[1:]
        failed = []
    atomic_json(output/'coronal-bounds.json', dict(schema='tongue-tv-actual-coronal-bounds/v1', sections=sections,
        failed_stations=failed, passed=not failed, original_negative_stations=[0] if posterior_end is not None else [],
        accepted_zero_corridor_boundary=posterior_end.tolist() if posterior_end is not None else None,
        input_bindings=bindings,
        entire_requested_AP_range_m=[float(lo), float(hi)], anatomical_acceptance=False))
    if failed:
        atomic_json(output/'construction-blocked.json', dict(schema='tongue-tv-corridor-construction/v1',
            passed=False, failed_stations=failed, reason='No continuous positive central corridor under the declared actual bounds; no forced block or narrowed retry generated',
            IL_corner_charts=il_rows, TV_constructed=False, installed=False, anatomical_acceptance=False,
            elapsed_seconds=time.monotonic()-began))
        print(json.dumps(dict(TV_constructed=False, failed_stations=failed,
                             IL_corner_charts_corrected=True, elapsed_seconds=time.monotonic()-began)), flush=True)
        return
    if refinement is not None:
        require(posterior_end is not None, 'Use the same explicitly accepted posterior closure')
        by_z = {float(v[2]):np.asarray(v) for v in controls}
        facet_z = sorted(set(sl.triangles[refinement['SL_faces'], :, 2].reshape(-1).tolist()))
        with np.load(refinement['witnesses'], allow_pickle=False) as data:
            witness_z = sorted(set(data['TV_triangles'][:, :, 2].reshape(-1).tolist()))
        records = []
        for z in sorted(set(facet_z+witness_z)):
            if not posterior_end[2]+1e-8 < z < hi-1e-8:
                continue
            if z in facet_z:
                # At a literal source vertex, the slicing plane coincides with
                # edges. The two 1-nm one-sided limits establish the continuous
                # facet bound without perturbing the stored event coordinates.
                left, left_row = bound_section(z-1e-9, whole, sl, il, muco)
                right, right_row = bound_section(z+1e-9, whole, sl, il, muco)
                require(left is not None and right is not None,
                        'A genuine nonpositive corridor appeared at a crossed-facet event')
                value = (left+right)/2
                value[2] = z
                records.append(dict(actual_facet_vertex_z_m=z, source_faces=refinement['SL_faces'],
                    one_sided_numeric_limit_m=1e-9, left=left_row, right=right_row,
                    continuous_bound_control=value.tolist()))
            else:
                value, row = bound_section(z, whole, sl, il, muco)
                require(value is not None, 'A genuine nonpositive corridor appeared at an actual witness ring')
                records.append(dict(actual_witness_ring_z_m=z, bound=row,
                                    continuous_bound_control=value.tolist()))
            by_z[z] = value
        controls = [by_z[z] for z in sorted(by_z)]
        atomic_json(output/'numerical-refinement.json', dict(schema='tongue-tv-facet-bound-refinement/v1',
            predecessor_crossing_report=str(refine_bound), predecessor_crossing_report_sha256=digest_file(refine_bound),
            literal_crossed_SL_face_ids=refinement['SL_faces'], actual_source_vertex_z_events_m=facet_z,
            predecessor_witness_ring_z_m=witness_z, records=records,
            original_24_positive_coronal_bounds_unchanged=True,
            posterior_root_and_all_anatomical_parameters_unchanged=True,
            radius_factor_search_performed=False, source_SL_and_shared_charts_modified=False,
            interpolation='PCHIP in physical Z, interpolating actual controls rather than the previous approximating uniform-parameter B-spline',
            installed=False, anatomical_acceptance=False))
    controls = np.asarray(controls)
    # Common smooth low-order loft of measured maximal-inscribed sections.
    # Endpoint radii vanish at explicit source AP extrema, with no planar disks.
    ends = np.array([controls[0], controls[-1]])
    ends[:, 2] = [lo, hi]
    ends[:, 3:] = 0
    if posterior_end is not None:
        ends[0] = posterior_end
    if refinement is not None:
        # Query the unchanged terminal Z. Copying the last interior station's
        # Y while advancing Z was an extrapolation bug, not an anatomical rule.
        actual_tip, tip_row = bound_section(float(hi), whole, sl, il, muco)
        require(actual_tip is not None, 'The unchanged anterior cap Z has no positive actual corridor')
        previous = ends[1].copy()
        ends[1, :3] = actual_tip[:3]
        atomic_json(output/'anterior-cap-correction.json', dict(
            schema='tongue-tv-actual-anterior-cap-centre/v1', previous_extrapolated_endpoint=previous.tolist(),
            corrected_endpoint=ends[1].tolist(), actual_corridor=tip_row,
            AP_extent_and_all_anatomical_parameters_unchanged=True,
            correction_reason='Query actual corridor at the cap Z instead of copying the last interior Y',
            installed=False, anatomical_acceptance=False))
    control = np.r_[ends[:1], controls, ends[1:]]
    count = len(control)
    knots = np.r_[np.zeros(4), np.linspace(0,1,count-2)[1:-1], np.ones(4)]
    if refinement is None:
        spline = BSpline(knots,control,3)
        rings = spline(np.linspace(0,1,49))[1:-1]
    else:
        spline = PchipInterpolator(control[:,2], control, axis=0)
        sample_z = np.unique(np.r_[np.linspace(control[0,2],control[-1,2],49), control[:,2]])
        rings = spline(sample_z)[1:-1]
    angle = np.arange(48)*2*np.pi/48
    xyz = np.stack([rings[:,0,None]+rings[:,3,None]*np.cos(angle),
                    rings[:,1,None]+rings[:,4,None]*np.sin(angle),
                    np.broadcast_to(rings[:,2,None],(len(rings),len(angle)))],axis=-1)
    vertices = np.r_[control[0,None,:3], xyz.reshape(-1,3), control[-1,None,:3]]
    faces = [[0,1+(j+1)%48,1+j] for j in range(48)]
    for i in range(len(rings)-1):
        a,b=1+i*48,1+(i+1)*48
        for j in range(48):
            k=(j+1)%48
            faces.extend([[a+j,a+k,b+k],[a+j,b+k,b+j]])
    base=1+(len(rings)-1)*48
    faces.extend([[base+j,base+(j+1)%48,len(vertices)-1] for j in range(48)])
    vertices=(vertices-ORIGIN).astype('<f4').astype(float)+ORIGIN
    mesh=Mesh(vertices,np.asarray(faces,dtype=np.int64),FRAME)
    if topology(mesh)['oriented_triangle_integral_m3'] < 0:
        mesh=Mesh(vertices,mesh.faces[:,[0,2,1]],FRAME)
    tri=mesh.triangles
    raw=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    normals=np.zeros_like(mesh.vertices)
    for c in range(3):np.add.at(normals,mesh.faces[:,c],raw)
    normals/=np.maximum(np.linalg.norm(normals,axis=1)[:,None],1e-30)
    words=np.rint(normals*32767).astype('<i2')
    mesh=Mesh(mesh.vertices,mesh.faces,FRAME,dict(normals=words,normalSigns=np.ones(len(vertices),dtype='i1')),
        dict(anatomical_extent='Representative central interwoven TV core, not exclusive histological tissue'))
    raw_path=output/'transverse-vertical-shared-raw-smooth-fields.npz'
    mesh.save(raw_path)
    mesh_path=output/'transverse-vertical-shared.npz'
    tv_row=chart_il(raw_path,mesh_path)
    mesh=Mesh.load(mesh_path)
    rows=[dict(key='superior-longitudinal',mesh=str(SOURCE/'part-02.npz'),sha256=digest_file(SOURCE/'part-02.npz'),
               geometry_sha256=sl.geometry_sha,unchanged_literal=True),*il_rows,
          tv_row]
    atomic_json(output/'construction.json',dict(schema='tongue-coronal-contour-loft/v1',meshes=rows,
        input_bindings=bindings,tissue_ownership_claim=False,exact_muscle_extent_claim=False,
        numerical_checks_pending=True,source_files_modified=False,skeleton_modified=False,
        original_tv_geometry=True,anatomical_acceptance=False,installed=False,
        quality=helper.quality(mesh),elapsed_seconds=time.monotonic()-began))
    print(json.dumps(dict(TV_constructed=True,triangles=len(mesh.faces),elapsed_seconds=time.monotonic()-began)),flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--posterior-termination',type=Path)
    parser.add_argument('--refine-bound',type=Path)
    args=parser.parse_args()
    main(args.output.resolve(),args.posterior_termination.resolve() if args.posterior_termination else None,
         args.refine_bound.resolve() if args.refine_bound else None)
