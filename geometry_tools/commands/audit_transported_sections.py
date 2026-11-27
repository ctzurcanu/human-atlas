"""Re-cut a transported surface instead of reusing stale plane-edge weights.

python -m geometry_tools.commands.audit_transported_sections CANDIDATE --output REPORT
This is an independent sampled-section check, not a complete 3D acceptance.
Exact coordinates alone join contours; open cuts are never capped or repaired.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import shapely

from geometry_tools.course_io import load_arrays, load_course_spec, sha, write_json
from geometry_tools.algorithms.literal_plane_sections import literal_plane_sections


def section_polygons(vertices, faces, y):
    cut = literal_plane_sections(vertices, faces, [0, y, 0], [0, 1, 0])
    polygons = [shapely.polygons(q[:, [0, 2]]) for q in cut['loops']]
    valid = [bool(shapely.is_valid(p)) for p in polygons]
    return polygons, dict(
        closed=cut['closed_local_sections'], loops=len(polygons),
        invalidLoops=[i for i, ok in enumerate(valid) if not ok],
        openComponents=len(cut['open_components']),
        coplanarFaces=len(cut['coplanar_source_faces']),
        unresolvedFaces=len(cut['unresolved_source_faces']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('candidate')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    started = time.monotonic()
    folder, output = Path(args.candidate).resolve(), Path(args.output).resolve()
    if output.exists():
        raise ValueError('Use a new immutable report path')
    report_path = folder / 'report.json'
    candidate = json.loads(report_path.read_text())
    if candidate['schema'] not in ['proper-source-profile-transport/v1','tube-storage-preview/v1']:
        raise ValueError('Expected complete proper-profile transport candidate')
    spec, _, _, plan, course, *_ = load_course_spec(candidate['specPath'])
    meshes = []
    for row in candidate['rows']:
        path = folder / row['geometryFile']
        if sha(path) != row['sha256']:
            raise ValueError('Candidate mesh changed')
        meshes.append((row['partId'], load_arrays(path)))
    clearance = plan['physicalClearance_m']
    sections, invalid, open_cuts, contacts = [], 0, 0, 0
    for station in course['rows']:
        context = {key: shapely.from_wkb(value)
                   for key, value in station['contextWKB'].items()}
        rows, selected = [], []
        for pid, mesh in meshes:
            polygons, result = section_polygons(mesh['vertices'], mesh['faces'], station['y_m'])
            _, before = section_polygons(mesh['source_vertices'], mesh['faces'], station['y_m'])
            invalid += len(result['invalidLoops'])
            open_cuts += int(not result['closed'])
            checks = []
            for i, polygon in enumerate(polygons):
                if i in result['invalidLoops']:
                    continue
                gaps = {key: float(shapely.distance(polygon, shape))
                        for key, shape in context.items() if key != 'container'}
                gaps['container'] = float(shapely.distance(polygon, shapely.boundary(context['container'])))
                fits = bool(shapely.contains(context['container'], polygon)) and all(
                    shapely.disjoint(polygon, shape) for key, shape in context.items()
                    if key != 'container') and min(gaps.values()) >= clearance - 1e-8
                contacts += int(not fits)
                checks.append(dict(loop=i, clear=bool(fits), gaps_m=gaps))
                selected.append((pid, i, polygon))
            rows.append(dict(partId=pid, source=before, candidate=result, context=checks))
        mutual = []
        for i, (pid, loop, a) in enumerate(selected):
            for qid, other, b in selected[i + 1:]:
                gap = float(shapely.distance(a, b))
                ok = bool(shapely.disjoint(a, b)) and gap >= clearance - 1e-8
                contacts += int(not ok)
                mutual.append(dict(parts=[pid, qid], loops=[loop, other], gap_m=gap, clear=bool(ok)))
        sections.append(dict(index=station['index'], y_m=station['y_m'], parts=rows, mutual=mutual))
        if time.monotonic() - started > 18:
            raise TimeoutError('18s native re-cut audit limit')
    result = dict(schema='transported-actual-section-audit/v1',
                  candidateReportSha256=sha(report_path), specSha256=sha(candidate['specPath']),
                  sourceCodeSha256=sha(__file__), stationsChecked=len(sections),
                  invalidContourCount=invalid, openCutCount=open_cuts,
                  clearanceFailureCount=contacts, allSampledSectionsClear=not (invalid or open_cuts or contacts),
                  sections=sections, oldPlaneWeightsReused=False, capsCreated=0,
                  exactCoordinateJoinsOnly=True, full3DSurfaceAuditPending=True,
                  installed=False, anatomicalAcceptance=False,
                  elapsedSeconds=time.monotonic() - started)
    write_json(output, result)
    print(json.dumps({key: value for key, value in result.items() if key != 'sections'}))
    if not result['allSampledSectionsClear']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
