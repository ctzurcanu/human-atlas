#!/usr/bin/env python3
"""Export measured uncovered interval segments without new spatial queries."""
import hashlib
import json
import time
from pathlib import Path

import numpy as np


def union(intervals):
    result = []
    for a, b in sorted(intervals):
        if result and a <= result[-1][1]:
            result[-1][1] = max(result[-1][1], b)
        else:
            result.append([float(a), float(b)])
    return result


def complement(intervals):
    end = 0.0
    for a, b in union(intervals):
        if a > end + 1e-12:
            yield end, a
        end = max(end, b)
    if end < 1 - 1e-12:
        yield end, 1.0


def main():
    start = time.monotonic()
    base = Path(".local-models/geometry-workbench/biodigital-male-circulation-v1/global-envelope-v1/fixed-skin/envelope-topology-v1")
    source = base / "geometric-seam-coverage-v1/report.json"
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    r = json.loads(source.read_text())
    output = base / "boundary-uncertainty-v1"
    output.mkdir(parents=True, exist_ok=True)
    arrays, summaries = {}, {}
    for scope, coverage_keys in {"primary": ["other_primary", "same_owner_nonincident_primary"], "rendered": ["all_nonincident_rendered"]}.items():
        segments, edge_ids, components, part_indices, part_ids, face_ids, source_faces, intervals = [], [], [], [], [], [], [], []
        for e in r["edges"]:
            xyz = np.asarray(e["rawEndpoints_m"], dtype=np.float64)
            covered = [pair for key in coverage_keys for pair in e["coverage"][key]["intervals"]]
            for a, b in complement(covered):
                segments.append([xyz[0] + a * (xyz[1] - xyz[0]), xyz[0] + b * (xyz[1] - xyz[0])])
                edge_ids.append(e["edgeIndex"])
                components.append(e["componentIndex"])
                part_indices.append(e["sourcePartIndex"])
                part_ids.append(e["sourcePartId"])
                face_ids.append(e["skinFaceId"])
                source_faces.append(e["originalSourceFaceId"])
                intervals.append([a, b])
        segments = np.asarray(segments, dtype=np.float64).reshape((-1, 2, 3))
        arrays.update({f"{scope}_segments_m": segments, f"{scope}_source_boundary_edge_index": np.asarray(edge_ids, dtype=np.int32), f"{scope}_boundary_component_index": np.asarray(components, dtype=np.int32), f"{scope}_source_part_index": np.asarray(part_indices, dtype=np.int32), f"{scope}_source_part_id": np.asarray(part_ids), f"{scope}_original_skin_face_id": np.asarray(face_ids, dtype=np.int64), f"{scope}_original_source_face_id": np.asarray(source_faces, dtype=np.int64), f"{scope}_source_edge_t_interval": np.asarray(intervals, dtype=np.float64)})
        component_scopes = []
        for ci in sorted(set(components)):
            points = segments[np.asarray(components) == ci].reshape((-1, 3))
            component_scopes.append({"componentIndex": ci, "bounds_m": [points.min(axis=0).tolist(), points.max(axis=0).tolist()], "segmentCount": int(sum(q == ci for q in components))})
        summaries[scope] = {"segments": len(segments), "distinctBoundaryEdges": len(set(edge_ids)), "totalUncoveredLength_m": float(np.linalg.norm(segments[:, 1] - segments[:, 0], axis=1).sum()), "componentScopes": component_scopes, "fullyCoveredHandNetworks5And79Omitted": not any(ci in (5, 79) for ci in components)}
    arrays["geometric_tolerance_m"] = np.asarray(r["geometricTolerance_m"])
    np.savez_compressed(output / "uncovered-boundary-segments.npz", **arrays)
    result = {"schema": "source-bound-uncovered-skin-boundary-intervals/v1", "inputReport": str(source), "inputReportSha256": digest, "skinInputBindings": r["inputBindings"], "intervalsAreActualRawSourceSegments": True, "noNewSpatialQueries": True, "geometricTolerance_m": r["geometricTolerance_m"], "primaryScope": "Use with literal 484-part primary-skin rays: union of other-primary and nonincident same-owner primary coverage; complement only", "renderedScope": "Use only if ray scene also contains the corresponding actual rendered appendage surfaces; excludes unused metric hair carrier", "classificationUse": "Per-component projected opening shadows, grazing flags, and local uncertainty; segments are not caps and are not observed skin clearance triangles", "globalClosednessCertified": False, "summaries": summaries, "elapsedSeconds": time.monotonic() - start}
    (output / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: {q: v for q, v in values.items() if q != "componentScopes"} for k, values in summaries.items()}))
    print(json.dumps({"elapsedSeconds": time.monotonic() - start, "sourceReportUnchanged": hashlib.sha256(source.read_bytes()).hexdigest() == digest}))


if __name__ == "__main__":
    main()
