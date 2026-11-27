"""Diagnose the rejected v22 at its unchanged centreline coordinates.

Transport comparisons are diagnostics, never replacement candidates. The
available HRA airway ends are mesh truncations, not arterial terminal targets.
"""
import json
import argparse
import time
from pathlib import Path

import numpy as np

from geometry_tools.common import digest_file
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.contact import self_contact_summary
from geometry_tools.algorithms.tube_sections import (
    arclength, path_tangents, bend_diagnostics, sweep_sections,
)

BASE = Path('.local-models/geometry-workbench/circulatory-completion/pulmonary-circulation-frontier-v1')
OLD = BASE / 'right-S2-terminal-bound-candidate-v22'
OUT = BASE / 'right-S2-station-fold-diagnosis-v25'


def normals(sections):
    q = sections - sections.mean(1)[:, None]
    n = np.cross(q, np.roll(q, -1, axis=1)).sum(1)
    return n / np.linalg.norm(n, axis=1)[:, None]


def metrics(sections, source_sections):
    c = sections.mean(1)
    t = path_tangents(c)
    source_t = path_tangents(source_sections.mean(1))
    source_sign = np.sign(np.sum(normals(source_sections) * source_t, axis=1))
    alignment = source_sign * np.sum(normals(sections) * t, axis=1)
    angle = np.degrees(np.arccos(np.clip(np.abs(alignment), -1, 1)))
    bend = bend_diagnostics(c, sections[None])
    s = arclength(c)
    step = np.diff(s)
    return dict(
        length_m=float(s[-1]), minimum_spacing_m=float(step.min()),
        maximum_spacing_m=float(step.max()),
        maximum_station_turn_deg=float(np.degrees(bend['station_turn_rad']).max()),
        maximum_ring_plane_vs_own_tangent_angle_deg=float(angle.max()),
        ring_plane_tilt_over_30deg=np.flatnonzero(angle > 30).tolist(),
        reversed_ring_orientation_stations=np.flatnonzero(alignment <= 0).tolist(),
        maximum_curvature_support=float(bend['curvature_support'].max()),
        minimum_sampled_sweep_jacobian=float(bend['minimum_sampled_sweep_jacobian']),
        tight_bend_stations=bend['failed_stations'].tolist(),
    ), dict(centres=c, own_tangents=t, signed_ring_alignment=alignment,
            ring_tangent_angle_deg=angle, curvature_support=bend['curvature_support'],
            station_spacing_m=step, station_turn_rad=bend['station_turn_rad'])


def transitions():
    began = time.monotonic()
    a = np.load(OLD / 'literal-proper-ring-transport.npz', allow_pickle=False)
    diagnosis = json.loads((OUT / 'evidence.json').read_text())
    c = a['canonical_source_centres']
    q = a['source_bronchial_projection']
    h = a['current_HRA_course_base']
    m = a['unprotected_offset_candidate']
    w = a['root_protection_weights']
    rows = []
    for row in diagnosis['rows']:
        key = row['key']
        source = a[key + '_source_sections'].mean(1)
        target = a[key + '_candidate_sections'].mean(1)
        ids = np.argmin(np.linalg.norm(source[:, None] - c[None], axis=2), axis=1)
        error = np.linalg.norm(source - c[ids], axis=1)
        assert error.max() < 1e-6
        ds = np.linalg.norm(np.diff(source, axis=0), axis=1)
        dt = np.linalg.norm(np.diff(target, axis=0), axis=1)
        dq = np.linalg.norm(np.diff(q[ids], axis=0), axis=1)
        dh = np.linalg.norm(np.diff(h[ids], axis=0), axis=1)
        dm = np.linalg.norm(np.diff(m[ids], axis=0), axis=1)
        rank = np.argsort(dt / ds)[::-1][:6]
        rows.append(dict(key=key, maximum_source_station_canonical_lookup_error_m=float(error.max()),
                         largest_step_expansions=[dict(stations=[int(i), int(i + 1)],
                              source_step_m=float(ds[i]), candidate_step_m=float(dt[i]),
                              ratio=float(dt[i] / ds[i]), source_guide_projection_step_m=float(dq[i]),
                              target_guide_step_m=float(dh[i]), unprotected_mapping_step_m=float(dm[i]),
                              root_protection_weights=w[ids[[i, i + 1]]].tolist(),
                              source_points_m=source[[i, i + 1]].tolist(),
                              candidate_points_m=target[[i, i + 1]].tolist()) for i in rank]))
    result = dict(schema='right-S2-rejected-guide-projection-transition-diagnosis/v25',
                  candidate_evidence_sha256=digest_file(OLD / 'candidate-evidence.json'),
                  own_tangent_diagnosis_sha256=digest_file(OUT / 'evidence.json'), rows=rows,
                  anatomy_not_inferred_from_projection=True, geometry_changed=False,
                  elapsed_seconds=time.monotonic() - began)
    p = OUT / 'projection-transition-evidence.json'
    p.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(output=str(p), sha256=digest_file(p), rows=rows,
                          elapsed_seconds=result['elapsed_seconds'])))


def main():
    began = time.monotonic()
    OUT.mkdir(exist_ok=True)
    evidence = json.loads((OLD / 'candidate-evidence.json').read_text())
    rejected = json.loads((OLD / 'root-source-versus-candidate-self-contact-v2.json').read_text())
    assert rejected['installation_rejected'] and rejected['new_candidate_defective_components'] == 6
    contract = json.loads((BASE / 'right-S2-ordered-chart-contract-v21/contract.json').read_text())
    family = json.loads((BASE / 'family-contract-v2.json').read_text())
    assert digest_file(Path(family['source_catalogue'])) == evidence['source_catalogue_sha256']
    atlas = json.loads(Path(family['source_catalogue']).read_text())
    parts = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    transport = np.load(OLD / 'literal-proper-ring-transport.npz', allow_pickle=False)
    arrays = {}
    result = []
    for row in contract['selected_source_curves']:
        key = row['key']
        src = transport[key + '_source_sections']
        v22 = transport[key + '_candidate_sections']
        centres = v22.mean(1)
        fragment_row = next(x for x in evidence['selected_fragment_meshes']
                            if (x['part_id'], x['component']) == (row['part_id'], row['component']))
        fragment = Mesh.load(fragment_row['mesh'])
        whole = reader.read_part(atlas, parts[row['part_id']], frame=fragment.frame)
        source_faces = np.asarray(fragment.source['source_face_ids'], int)
        used, inv = np.unique(whole.faces[source_faces], return_inverse=True)
        assert np.array_equal(inv.reshape(-1, 3), fragment.faces)
        lookup = {tuple(p): (si, ri) for si, ring in enumerate(src) for ri, p in enumerate(ring)}
        station_ring = np.array([lookup[tuple(p)] for p in whole.vertices[used]], int)
        assert np.max(np.linalg.norm(v22[station_ring[:, 0], station_ring[:, 1]] - fragment.vertices, axis=1)) < 1e-12
        certificate = next(x for x in rejected['rows']
                           if (x['part_id'], x['component']) == (row['part_id'], row['component']))
        r = dict(part_id=row['part_id'], component=row['component'], key=key, stations=len(src),
                 v22_self_defect_positive_certificate=certificate['candidate_positive_certificate']['defect_pairs'],
                 literal_source_self_audit_complete=certificate['source_self_summary']['complete_within_pair_budget'],
                 literal_source_self_defects=certificate['source_self_summary']['defect_pairs'], phases={})
        for phase, sections in [('source', src), ('v22', v22)]:
            summary, values = metrics(sections, src)
            r['phases'][phase] = summary
            for name, value in values.items():
                arrays[key + '_' + phase + '_' + name] = value
        examples = []
        for first, second in certificate['candidate_positive_certificate']['examples'][:8]:
            aa = sorted(set(station_ring[fragment.faces[first], 0].tolist()))
            bb = sorted(set(station_ring[fragment.faces[second], 0].tolist()))
            examples.append(dict(fragment_faces=[first, second], source_faces=source_faces[[first, second]].tolist(),
                                 first_stations=aa, second_stations=bb))
        r['exact_self_crossing_face_to_station_witnesses'] = examples
        # One comparison at exactly the already-rejected centreline: double-
        # reflection RMF based on this vessel's OWN tangent field. No fit.
        try:
            own, frames, report = sweep_sections(src[None], centres, dtype=np.float64)
            own = own[0]
            summary, values = metrics(own, src)
            arrays[key + '_own_RMF_diagnostic_sections'] = own
            arrays[key + '_own_RMF_frames'] = frames
            for name, value in values.items():
                arrays[key + '_own_RMF_' + name] = value
            vertices = own[station_ring[:, 0], station_ring[:, 1]]
            diagnostic = Mesh(vertices, fragment.faces, fragment.frame,
                              source=dict(diagnostic_only=True, candidate_v22_centres_unchanged=True,
                                          anatomical_acceptance=False, source_part_id=row['part_id'],
                                          source_component=row['component'], source_face_ids=source_faces.tolist()))
            summary['self_audit'], masks = self_contact_summary(diagnostic, max_pairs=100000, max_tests=500000, stop_on_defect=True)
            for side, values in masks.items():
                arrays[key + '_own_RMF_self_' + side + '_face_ids'] = values
            summary['section_transport_report'] = report
            summary['maximum_centre_change_m'] = float(np.linalg.norm(own.mean(1) - centres, axis=1).max())
            summary['root_profile_error_not_freeze_repaired_m'] = float(np.linalg.norm(own[-1] - src[-1], axis=1).max()) if row['component'] == 1 else None
            path = OUT / (key + '-own-tangent-diagnostic.npz')
            if path.exists():
                assert Mesh.load(path).geometry_sha == diagnostic.geometry_sha
            else:
                diagnostic.save(path)
            summary['diagnostic_mesh'] = str(path.resolve())
            r['phases']['own_RMF_at_same_centres'] = summary
        except ValueError as exc:
            r['own_RMF_diagnostic_error'] = str(exc)
        result.append(r)
    np.savez_compressed(OUT / 'literal-station-diagnosis.npz', **arrays)
    data = dict(schema='right-S2-rejected-course-own-tangent-diagnosis/v25',
                candidate_evidence_sha256=digest_file(OLD / 'candidate-evidence.json'),
                independent_self_rejection_sha256=digest_file(OLD / 'root-source-versus-candidate-self-contact-v2.json'),
                source_catalogue_sha256=evidence['source_catalogue_sha256'], rows=result,
                airway_ends_are_available_mesh_truncations_not_physiological_arterial_terminal_targets=True,
                previous_terminal_bijection_not_anatomical_requirement=True,
                no_target_offset_field_or_centreline_modified=True,
                diagnostic_transports_are_not_candidates=True, contacts_held=True,
                catalogue_changed=False, installed=False, anatomical_acceptance=False,
                arrays_sha256=digest_file(OUT / 'literal-station-diagnosis.npz'),
                elapsed_seconds=time.monotonic() - began)
    (OUT / 'evidence.json').write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(dict(output=str(OUT), evidence_sha256=digest_file(OUT / 'evidence.json'),
                          elapsed_seconds=data['elapsed_seconds'],
                          rows=[dict(key=x['key'], phases={k: {a:b for a,b in v.items() if a in [
                              'maximum_ring_plane_vs_own_tangent_angle_deg', 'maximum_curvature_support',
                              'maximum_station_turn_deg', 'self_audit']} for k,v in x['phases'].items()}) for x in result])))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['diagnosis', 'transitions'], default='diagnosis')
    args = parser.parse_args()
    (main if args.stage == 'diagnosis' else transitions)()
