#!/usr/bin/env python3
"""Literal-triangle context audit for a proposed local clavicular arterial course.

Run under geometry_tools.runtime with a 40-second wall deadline. This only
reads anatomy and writes evidence; a right audit reflects the candidate,
never the native right context or the live model.
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
from geometry_tools.contact import pair_contacts
from geometry_tools.model import AtlasReader, Mesh
from geometry_tools.spatial import BVH

WORK = ROOT / '.local-models/geometry-workbench/circulatory-completion'
CONTEXT = WORK / 'clavicular-branch-bilateral-context-audit-v1/active-served-corridor-context.json'
FRAME = 'model:local-ta98:assembled-metres'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, report):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.partial')
    temporary.write_text(json.dumps(report, indent=2) + '\n')
    temporary.replace(path)


def reflect(mesh):
    vertices = mesh.vertices.copy()
    vertices[:, 0] *= -1
    attributes = {k: v.copy() for k, v in mesh.attributes.items()}
    if 'normals' in attributes:
        attributes['normals'][:, 0] *= -1
    return Mesh(vertices, mesh.faces[:, [0, 2, 1]], mesh.frame, attributes,
                dict(kind='audit-only-X0-reflection', original_geometry_sha256=mesh.geometry_sha))


def role(record, side):
    pid = record['id']
    suffix = '.l' if side == 'left' else '.r'
    if pid == 'ZA:Thoraco-acromial artery' + suffix:
        return 'source-parent-local-takeoff', 'declared parent; contact must remain in explicit takeoff patch'
    if pid == 'ZA:Subclavius muscle' + suffix:
        return 'transmuscular-recipient', 'muscular entry/course explicitly qualified; no lumen/perfusion proof'
    if 'Pectoralis major' in pid or 'pectoralis major' in pid:
        return 'pectoral-muscular-context', 'recipient muscular context, exact entry and distribution require anatomical review'
    if 'sternoclavicular' in pid.lower() and ('capsul' in pid.lower()):
        return 'terminal-capsular-recipient', 'capsular terminal entry is allowed only at its declared actual source selector'
    if record.get('system') == 'skeletal' or 'cartilage' in pid.lower():
        return 'hard-obstacle', 'fixed bone/cartilage; no crossing permitted'
    if record.get('system') in ('venous', 'nervous'):
        return 'vessel-nerve-obstacle', 'fixed vein/nerve; no crossing permitted'
    if record.get('system') == 'arterial':
        return 'other-artery', 'fixed artery; no undeclared crossing or junction'
    if record.get('system') in ('fascia', 'connective'):
        return 'sheet-ligament-context', 'record actual sheet contacts; no closed-volume assumption'
    if record.get('system') == 'respiratory':
        return 'respiratory-obstacle', 'fixed respiratory surface; no crossing permitted'
    return 'other-soft-context', 'actual contact requires anatomical review'


def audit(args):
    started = time.monotonic()
    binding = json.loads(CONTEXT.read_text())
    body_path = Path(binding['served_body'])
    catalogue_path = Path(binding['source_catalogue'])
    assert digest(body_path) == binding['served_body_sha256'], 'Served body changed; regenerate context bindings'
    assert digest(catalogue_path) == binding['source_catalogue_sha256'], 'Master changed; regenerate context bindings'
    body = json.loads(body_path.read_text())
    paths = {'main': args.course_dir / 'main-outer-course.npz',
             'muscular': args.course_dir / 'muscular-outer-course.npz'}
    courses = {name: Mesh.load(path, frame=FRAME) for name, path in paths.items()}
    if args.side == 'right':
        courses = {name: reflect(mesh) for name, mesh in courses.items()}
    trees = {name: BVH(mesh) for name, mesh in courses.items()}
    low, high = np.asarray(binding['roi_m'])
    if args.side == 'right':
        low[0], high[0] = -high[0], -low[0]
    reader = AtlasReader()
    hidden = set(body.get('defaultHidden', []))
    # Native context is selected independently on each side, including midline
    # records. No counterpart mapping or reflected ligament is substituted.
    records = [p for p in body['parts'] if p['id'] not in hidden
               and not p.get('defaultHidden') and not p.get('suppressed')
               and 'bounds' in p and np.all(np.asarray(p['bounds'])[0] <= high)
               and np.all(np.asarray(p['bounds'])[1] >= low)]
    spec_path = args.course_dir / 'course-spec.json'
    spec = json.loads(spec_path.read_text()) if spec_path.exists() else {}
    allowed_parent = spec.get('allowed_parent_face_ids')
    rows = []
    report = dict(schema='clavicular-course-native-context-actual-triangle-contacts/v2',
                  complete=False, side=args.side, frame=FRAME,
                  implementation_sha256=digest(__file__),
                  source_catalogue_sha256=digest(catalogue_path),
                  served_body_sha256=digest(body_path), context_binding_sha256=digest(CONTEXT),
                  course_npz_sha256={name: digest(path) for name, path in paths.items()},
                  course_geometry_sha256={name: mesh.geometry_sha for name, mesh in courses.items()},
                  course_spec_sha256=digest(spec_path) if spec_path.exists() else None,
                  candidate_transform='identity' if args.side == 'left' else 'X=0 reflection, reversed winding; audit only',
                  native_context_transform='identity; actual current independently selected native sources',
                  roi_m=[low.tolist(), high.tolist()], distance_threshold_m=.0002,
                  backend='Float64 literal BVH plus triangle_pairs; no normalization/remeshing',
                  broad_record_count=len(records), rows=rows, model_geometry_DB_modified=False)
    for index, record in enumerate(records):
        mesh = reader.read_part(body, record, frame=FRAME, catalogue_path=body_path,
                                catalogue_sha=binding['served_body_sha256'])
        triangles = mesh.triangles
        roi_faces = np.flatnonzero(np.all(triangles.min(axis=1) <= high, axis=1)
                                  & np.all(triangles.max(axis=1) >= low, axis=1))
        if not len(roi_faces):
            continue
        tree = BVH(mesh)
        kind, qualification = role(record, args.side)
        for name, course in courses.items():
            arrays, summary = pair_contacts(course, mesh, .0002, first_bvh=trees[name],
                                             second_bvh=tree, max_pairs=100000, max_tests=1000000)
            row = dict(context_part_id=record['id'], course=name, system=record.get('system'),
                       role=kind, qualification=qualification, summary=summary,
                       context_geometry_sha256=mesh.geometry_sha,
                       context_attributes_sha256=mesh.attributes_sha,
                       source={k: v for k, v in mesh.source.items()
                               if k not in ('source_vertex_ids', 'source_face_ids')})
            if len(arrays['pairs']):
                strict = np.isin(arrays['relation'], ['crossing', 'coplanar-overlap'])
                strict_context_ids = np.unique(arrays['pairs'][strict, 1])
                row.update(actual_context_face_ids=np.unique(arrays['pairs'][:, 1]).tolist(),
                           strict_context_face_ids=strict_context_ids.tolist(),
                           course_face_ids=np.unique(arrays['pairs'][:, 0]).tolist(),
                           witness_bounds_m=[arrays['first'].min(axis=0).tolist(), arrays['first'].max(axis=0).tolist()],
                           minimum_returned_triangle_distance_m=float(arrays['distance_m'].min()),
                           strict_witnesses=[dict(course_face=int(arrays['pairs'][j, 0]),
                                                  context_face=int(arrays['pairs'][j, 1]), relation=str(arrays['relation'][j]),
                                                  course_point_m=arrays['first'][j].tolist(),
                                                  context_point_m=arrays['second'][j].tolist())
                                             for j in np.flatnonzero(strict)[:12]])
                if kind == 'source-parent-local-takeoff':
                    root = course.vertices[:16].mean(axis=0)
                    row['maximum_contact_distance_from_first_profile_centre_m'] = float(np.linalg.norm(arrays['first'] - root, axis=1).max())
                    row['strict_faces_within_declared_parent_patch'] = (set(map(int, strict_context_ids)) <= set(allowed_parent)) if allowed_parent is not None else None
            rows.append(row)
        report.update(records_completed=index + 1, elapsed_seconds=time.monotonic() - started)
        save(args.output, report)
    strict = [dict(part_id=r['context_part_id'], course=r['course'], role=r['role'],
                   crossing=r['summary']['crossing_face_pairs'], coplanar=r['summary']['coplanar_overlap_face_pairs'],
                   bounds_m=r.get('witness_bounds_m')) for r in rows
              if r['summary']['crossing_face_pairs'] or r['summary']['coplanar_overlap_face_pairs']]
    forbidden = [r for r in strict if r['role'] in ('hard-obstacle', 'vessel-nerve-obstacle', 'other-artery', 'respiratory-obstacle')]
    report.update(complete=True, elapsed_seconds=time.monotonic() - started,
                  actual_context_count=len(rows) // len(courses),
                  all_pair_budgets_complete=all(r['summary']['complete_within_pair_budget'] for r in rows),
                  strict_contacts=strict, strict_forbidden_contacts=forbidden,
                  strict_forbidden_crossing_pairs=sum(r['crossing'] for r in forbidden),
                  strict_forbidden_coplanar_pairs=sum(r['coplanar'] for r in forbidden),
                  qualification='Exterior triangle contacts only; no lumen, containment, distribution or perfusion acceptance. Muscular and capsular contacts remain explicitly qualified. Independent native right context is used for right audit. No live anatomy changed.')
    save(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k != 'rows'}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--course-dir', type=Path, required=True)
    parser.add_argument('--side', choices=['left', 'right'], default='left')
    parser.add_argument('--output', type=Path, required=True)
    audit(parser.parse_args())
