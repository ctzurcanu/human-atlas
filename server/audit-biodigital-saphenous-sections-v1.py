"""Actual axial calibre and inherited/new skin-crossing diagnosis; no fitting."""
import hashlib
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.oriented_rays import OrientedRayIndex

W = ROOT / '.local-models/geometry-workbench/biodigital-male-circulation-v1'


def axial_loops(triangles, y):
    """Intersect literal triangles, report every closed loop; never bridge gaps."""
    t = triangles[(triangles[:, :, 1].min(1) < y) & (triangles[:, :, 1].max(1) > y)]
    segments = []
    for tri in t:
        points = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            if (tri[i, 1] - y) * (tri[j, 1] - y) < 0:
                points.append((tri[i] + (y-tri[i, 1])/(tri[j, 1]-tri[i, 1])*(tri[j]-tri[i]))[[0, 2]])
        if len(points) == 2:
            segments.append(points)
    if not segments:
        return [], {'segments': 0, 'openOrBranchedVertices': 0}
    p, inverse = np.unique(np.round(np.asarray(segments).reshape(-1, 2), 10), axis=0, return_inverse=True)
    edges = inverse.reshape(-1, 2)
    edges = edges[edges[:, 0] != edges[:, 1]]
    graph = {}
    for a, b in edges:
        graph.setdefault(int(a), []).append(int(b))
        graph.setdefault(int(b), []).append(int(a))
    bad = {a: len(b) for a, b in graph.items() if len(b) != 2}
    if bad:
        return [], {'segments': len(segments), 'openOrBranchedVertices': len(bad)}
    unseen, loops = set(graph), []
    while unseen:
        first = min(unseen)
        ids, previous, current = [first], first, graph[first][0]
        while current != first:
            ids.append(current)
            nxt = next(x for x in graph[current] if x != previous)
            previous, current = current, nxt
            assert len(ids) <= len(graph)
        unseen -= set(ids)
        q = p[ids]
        area = abs(float(np.sum(q[:, 0]*np.roll(q[:, 1], -1)-q[:, 1]*np.roll(q[:, 0], -1))/2))
        loops.append({'area_m2': area, 'centroidXZ_m': q.mean(0).tolist(), 'vertices': len(q)})
    loops.sort(key=lambda x: -x['area_m2'])
    return loops, {'segments': len(segments), 'openOrBranchedVertices': 0}


def main():
    began = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate-stage', default='saphenous-envelope-v1')
    parser.add_argument('--output-name', default='actual-section-and-crossing-diagnosis-v1')
    args = parser.parse_args()
    assert all(Path(v).name == v and v not in ('.', '..') for v in vars(args).values())
    OUT = W/args.candidate_stage/args.output_name
    OUT.mkdir(exist_ok=False)
    ownership = json.loads((W/'lower-limb-joints-v4/candidate-v1/instance-ownership.json').read_text())['instances']
    by_index = {r['instanceIndex']: r for r in ownership}
    sources = []
    for stage in ('lower-limb-joints-v4', args.candidate_stage):
        with np.load(W/stage/'candidate-v1/venous.npz') as z:
            sources.append((z['vertices'], z['faces']))
    rows = []
    for part in ownership:
        if part['name'] not in [f'{side} {name}' for side in ('Left', 'Right') for name in ('femoral vein', 'great saphenous vein', 'femoral canal perforator veins')]:
            continue
        sections = []
        for y in np.linspace(.558371, .711371, 32):
            pair = []
            for v, f in sources:
                tri = v[f[part['candidateFaceStart']:part['candidateFaceStart']+part['candidateFaceCount']]]
                loops, topology = axial_loops(tri, y)
                pair.append({'loops': loops, 'topology': topology})
            areas = [sum(q['area_m2'] for q in p['loops']) for p in pair]
            ratio = areas[1]/areas[0] if min(areas) > 0 and len(pair[0]['loops']) == len(pair[1]['loops']) else None
            sections.append({'y_m': float(y), 'before': pair[0], 'after': pair[1], 'totalClosedSectionAreaRatio': ratio})
        ratios = [r['totalClosedSectionAreaRatio'] for r in sections if r['totalClosedSectionAreaRatio'] is not None]
        rows.append({'sourceId': part['sourceInstanceId'], 'name': part['name'], 'sections': sections,
                     'closedSectionAreaRatioMinMax': [min(ratios), max(ratios)] if ratios else None})
    evidence_path = W/args.candidate_stage/'moved-triangle-skin-regression-v1/actual-triangle-evidence.npz'
    with np.load(evidence_path) as z:
        before_faces, after_faces = set(z['before_pairs'][:, 0]), set(z['after_pairs'][:, 0])
        new_ids = np.asarray(sorted(after_faces-before_faces), dtype=int)
        tb, ta = z['before_moved_triangles'][new_ids], z['after_moved_triangles'][new_ids]
        owners, faces = z['source_instance_indices'][new_ids], z['source_face_ids'][new_ids]
    benchmark = json.loads((W/'body-envelope-fit-v2/current-envelope-ray-benchmark-v1/report.json').read_text())
    skin_path = Path(benchmark['skin_path'])
    assert hashlib.sha256(skin_path.read_bytes()).hexdigest() == benchmark['skin_sha256']
    boundary_path = Path(benchmark['boundary_evidence']['path'])
    boundary = np.load(boundary_path)
    boundary_report = json.loads((boundary_path.parent/'report.json').read_text())
    index = OrientedRayIndex(Mesh.load(skin_path), nthreads=2,
        uncertainty_segments=boundary['primary_segments_m'], uncertainty_component_ids=boundary['primary_boundary_component_index'],
        boundary_tolerance_m=boundary_report['geometricTolerance_m'], boundary_evidence=benchmark['boundary_evidence'])
    def witnesses(t):
        return np.concatenate([t, (t+np.roll(t, -1, axis=1))/2, t.mean(1)[:, None]], axis=1)
    queries = np.r_[witnesses(tb).reshape(-1, 3), witnesses(ta).reshape(-1, 3)]
    result, proof = index.classify(queries, np.asarray(benchmark['directions']),
        near_surface_m=benchmark['near_surface_m'], min_reliable_rays=benchmark['min_reliable_rays'])
    n = len(new_ids)*7
    cb, ca = result['classification'][:n].reshape(-1, 7), result['classification'][n:].reshape(-1, 7)
    crossed = []
    for i, packed in enumerate(new_ids):
        crossed.append({'movedFaceIndex': int(packed), 'sourceId': by_index[int(owners[i])]['sourceInstanceId'],
            'sourceFace': int(faces[i]), 'beforeWitnessClassification': cb[i].tolist(), 'afterWitnessClassification': ca[i].tolist(),
            'wasPreviouslyStrictlyCrossing': False, 'beforeAllSevenInside': bool(np.all(cb[i] == 1)),
            'beforeAllSevenOutside': bool(np.all(cb[i] == 0)), 'afterHasOutsideWitness': bool(np.any(ca[i] == 0))})
    report = {'schema': 'actual-axial-section-and-new-crossing-diagnosis/v1', 'sectionPlane': 'fixed world Y; no smoothing or proxy cylinder',
        'loopEndpointQuantization_m': 1e-10, 'lumenQualification': 'Closed source wall contours enclose section area; do not infer a second wall or patent lumen from a single contour.',
        'parts': rows, 'newCrossingFaces': crossed, 'newCrossingFacesPreviouslyAllInside': sum(x['beforeAllSevenInside'] for x in crossed),
        'newCrossingFacesPreviouslyAllOutside': sum(x['beforeAllSevenOutside'] for x in crossed),
        'rayProof': proof, 'skinSha256': benchmark['skin_sha256'], 'geometryModified': False, 'anatomicalAcceptance': False,
        'elapsedSeconds': time.perf_counter()-began}
    (OUT/'report.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savez(OUT/'classified-crossing-witnesses.npz', before_triangles=tb, after_triangles=ta, queries=queries, **result)
    print(json.dumps({'parts': [{k:r[k] for k in ('name', 'closedSectionAreaRatioMinMax')} for r in rows],
        'newCrossingsPreviouslyAllInside': report['newCrossingFacesPreviouslyAllInside'],
        'newCrossingsPreviouslyAllOutside': report['newCrossingFacesPreviouslyAllOutside'], 'elapsedSeconds': report['elapsedSeconds']}))


if __name__ == '__main__':
    main()
