#!/usr/bin/env python3
"""Audit whether a proposed branch leaves a native arterial side wall or mouth.

Read-only anatomy diagnostic. Run under geometry_tools.runtime --seconds 20.
The fan closures only classify the source blood-volume envelope; they are
diagnostic surfaces, never arterial tissue or replacement geometry.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.algorithms.native_tree_exterior import _loops
from geometry_tools.audit import topology, winding_contains
from geometry_tools.contact import pair_contacts
from geometry_tools.model import AtlasReader, Mesh

WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion'
FRAME = 'model:local-ta98:assembled-metres'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def diagnostic(args):
    started = time.monotonic()
    binding_path = WORK / 'clavicular-branch-bilateral-context-audit-v1/active-served-corridor-context.json'
    binding = json.loads(binding_path.read_text())
    body_path = Path(binding['served_body'])
    assert digest(body_path) == binding['served_body_sha256'], 'Served anatomy changed'
    assert digest(binding['source_catalogue']) == binding['source_catalogue_sha256'], 'Master changed'
    atlas = json.loads(body_path.read_text())
    records = {p['id']: p for p in atlas['parts']}
    reader = AtlasReader()
    side = '.l' if args.side == 'left' else '.r'
    load = lambda pid: reader.read_part(atlas, records[pid], frame=FRAME,
                                       catalogue_path=body_path, catalogue_sha=binding['served_body_sha256'])
    parent = load('ZA:Thoraco-acromial artery' + side)
    course_path = args.course_dir / 'main-outer-course.npz'
    main = Mesh.load(course_path, frame=FRAME)
    if args.side == 'right':
        vertices = main.vertices.copy()
        vertices[:, 0] *= -1
        main = Mesh(vertices, main.faces[:, [0, 2, 1]], FRAME)
    loops = _loops(parent.faces)
    vertices = parent.vertices.tolist()
    cap_faces, mouths = [], []
    for index, loop in enumerate(loops):
        centre_id = len(vertices)
        vertices.append(parent.vertices[loop].mean(axis=0).tolist())
        cap_faces.extend([[int(loop[(i + 1) % len(loop)]), int(loop[i]), centre_id]
                          for i in range(len(loop))])
        mouths.append(dict(mouth_index=index, source_control_ids=loop.tolist(),
                           centre_m=vertices[-1], control_count=len(loop)))
    vertices = np.asarray(vertices)
    caps = Mesh(vertices, np.asarray(cap_faces), FRAME)
    classifier = Mesh(vertices, np.vstack([parent.faces, caps.faces]), FRAME)
    assert topology(classifier)['closed_oriented_manifold'], 'Native diagnostic classifier is not closed'
    lateral, lateral_summary = pair_contacts(main, parent, 0., max_pairs=100000, max_tests=1000000)
    cap_contacts, cap_summary = pair_contacts(main, caps, 0., max_pairs=100000, max_tests=1000000)
    centres = main.vertices.reshape(-1, args.profile_controls, 3).mean(axis=1)
    centres_inside = winding_contains(classifier, centres)
    root_inside = winding_contains(classifier, main.vertices[:args.profile_controls])
    cap_strict = np.flatnonzero(np.isin(cap_contacts['relation'], ['crossing', 'coplanar-overlap']))
    targets = []
    target_ids = ['ZA:Axillary vein' + side, 'ZA:Pectoral branches of thoraco-acromial artery' + side,
                  'ZA:Left subclavian vein' if args.side == 'left' else 'ZA:Right subclavian vein']
    for pid in target_ids:
        if pid not in records:
            continue
        other = load(pid)
        _, baseline = pair_contacts(parent, other, 0., max_pairs=100000, max_tests=1000000)
        current, current_summary = pair_contacts(main, other, 0., max_pairs=100000, max_tests=1000000)
        ids = np.flatnonzero(np.isin(current['relation'], ['crossing', 'coplanar-overlap']))
        inside = winding_contains(classifier, current['first'][ids]) if len(ids) else dict(inside=[], boundary_ambiguous=[])
        targets.append(dict(part_id=pid, native_parent_baseline=baseline, proposed_course=current_summary,
                            strict_witness_count=len(ids), strict_witness_inside_native_parent=sum(inside['inside']),
                            strict_witness_ambiguous=sum(inside['boundary_ambiguous']),
                            all_course_strict_witnesses_inside_native_parent=all(inside['inside']),
                            strict_witnesses_m=current['first'][ids].tolist()))
    passed = (lateral_summary['crossing_face_pairs'] > 0
              and cap_summary['returned_triangle_pairs'] == 0
              and all(root_inside['inside'])
              and lateral_summary['complete_within_pair_budget'] and cap_summary['complete_within_pair_budget'])
    report = dict(schema='clavicular-native-side-root-interface-diagnostic/v2', side=args.side,
                  implementation_sha256=digest(__file__), source_catalogue_sha256=binding['source_catalogue_sha256'],
                  served_body_sha256=binding['served_body_sha256'], input_course_sha256=digest(course_path),
                  course_geometry_sha256=main.geometry_sha, native_parent_geometry_sha256=parent.geometry_sha,
                  native_original_mouths=mouths,
                  first_profile_centre_m=centres[0].tolist(),
                  first_to_second_profile_direction=(centres[1] - centres[0]).tolist(),
                  actual_first_profile_points_m=main.vertices[:args.profile_controls].tolist(),
                  first_profile_containment=root_inside,
                  main_centreline_inside_parent_indices=np.flatnonzero(centres_inside['inside']).tolist(),
                  lateral_source_parent_contacts=lateral_summary,
                  original_native_mouth_diagnostic_fan_contacts=cap_summary,
                  mouth_strict_contacts=[dict(main_face=int(cap_contacts['pairs'][i, 0]),
                                               diagnostic_fan_face=int(cap_contacts['pairs'][i, 1]),
                                               point_m=cap_contacts['first'][i].tolist()) for i in cap_strict],
                  forbidden_interface_diagnostics=targets,
                  geometric_side_fork_gate_pass=passed, geometry_modified=False,
                  elapsed_seconds=time.monotonic() - started,
                  qualification='No anatomy acceptance or lumen/perfusion claim. A vein crossing buried within the existing parent may disappear during union; actual final exterior must be rechecked. Original-mouth fans are diagnostic only.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items()
                      if k not in ('actual_first_profile_points_m', 'mouth_strict_contacts', 'forbidden_interface_diagnostics')}, indent=2))
    print('Obstacle diagnostic:', [(t['part_id'], t['native_parent_baseline']['crossing_face_pairs'],
                                    t['proposed_course']['crossing_face_pairs'],
                                    t['strict_witness_inside_native_parent'], t['strict_witness_count']) for t in targets])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--course-dir', type=Path, required=True)
    parser.add_argument('--side', choices=['left', 'right'], default='left')
    parser.add_argument('--profile-controls', type=int, default=16)
    parser.add_argument('--output', type=Path, required=True)
    diagnostic(parser.parse_args())
