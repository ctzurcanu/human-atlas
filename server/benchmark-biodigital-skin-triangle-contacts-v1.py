#!/usr/bin/env python3
"""Bounded actual vessel/skin triangle BVH benchmark; no containment or fit."""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.distance import exact_intersection_segment, exact_transverse_crossing, triangle_pairs
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-stage", required=True)
    parser.add_argument("--output", default="skin-triangle-benchmark-v1")
    parser.add_argument("--faces-per-instance", type=int, default=19)
    parser.add_argument("--max-faces", type=int, default=32768)
    parser.add_argument("--distance-m", type=float, default=.003)
    parser.add_argument("--max-tests", type=int, default=500000)
    parser.add_argument("--max-pairs", type=int, default=150000)
    parser.add_argument("--soft-seconds", type=float, default=32)
    args = parser.parse_args()
    if not 1 <= args.max_faces <= 32768 or not 1 <= args.faces_per_instance <= 19:
        raise ValueError("This benchmark is limited to 32768 faces / 19 per instance")
    started = time.monotonic()
    base = ROOT / ".local-models/geometry-workbench/biodigital-male-circulation-v1"
    candidate = base / args.candidate_stage / "candidate-v1"
    output = base / args.candidate_stage / args.output
    output.mkdir(exist_ok=False)
    fixed = base / "global-envelope-v1/fixed-skin"
    inputs = [candidate / "candidate-report.json", candidate / "instance-ownership.json", fixed / "skin.npz", fixed / "policy.json", fixed / "ownership.json", base / args.candidate_stage / "contract.json"]
    bindings = {str(p): sha(p) for p in inputs}
    contract = json.loads(inputs[5].read_text())
    candidate_report = json.loads(inputs[0].read_text())
    assert bindings[str(inputs[5])] == candidate_report["contractSha256"]
    rows = json.loads(inputs[1].read_text())["instances"]
    skin_z = np.load(fixed / "skin.npz")
    skin = Mesh(skin_z["vertices"], skin_z["faces"], str(skin_z["frame"]))
    skin_source = json.loads(str(skin_z["source_json"]))
    identity_frame_binding = None
    skin_owner = skin_z["face_owner_part_index"]
    skin_source_face = skin_z["face_owner_source_face_id"]
    skin_ownership = json.loads((fixed / "ownership.json").read_text())
    # The saved source ownership structure is included separately, unmodified.
    sample_triangles, sample_rows = [], []
    mesh_names = ["arterial", "venous", "arterial-remnant-companion"]
    for name in mesh_names:
        path, ancestry_path = candidate / (name + ".npz"), candidate / (name + "-ancestry.npz")
        bindings[str(path)], bindings[str(ancestry_path)] = sha(path), sha(ancestry_path)
        z, az = np.load(path), np.load(ancestry_path)
        vertices, faces = z["vertices"], z["faces"]
        source_faces, source_instances = az["candidate_face_source_face"], az["candidate_face_source_instance"]
        source_ids = az["source_instance_ids"]
        if str(z["frame"]) != skin.frame:
            # This candidate contract names the exact fixed skin bytes as its
            # recipient geometry. Bind the identity alias; do not assume it
            # from similar dimensions, and do not change either saved mesh.
            assert contract["globalBodyRegistration"]["targetSkinSha256"] == bindings[str(fixed / "skin.npz")]
            assert contract["target_catalogue_sha256"] == skin_source["catalogue_sha256"]
            assert str(z["frame"]) == "model:local-ta98:assembled-metres" and skin.frame == "active-male-metres"
            identity_frame_binding = {"candidateFrame": str(z["frame"]), "fixedSkinFrame": skin.frame, "mapping": np.eye(4).tolist(), "basis": "Hash-bound candidate contract explicitly uses these exact target skin bytes in unchanged active catalogue coordinates", "contract": str(inputs[5]), "contractSha256": bindings[str(inputs[5])], "geometryChanged": False}
        for row in rows:
            if row["candidateMesh"] != name + ".npz":
                continue
            start, count = row["candidateFaceStart"], row["candidateFaceCount"]
            ids = start + np.unique(np.linspace(0, count - 1, min(args.faces_per_instance, count)).astype(np.int64))
            sample_triangles.extend(vertices[faces[ids]])
            for face in ids:
                source_id = str(source_ids[source_instances[face]])
                assert source_id == row["sourceInstanceId"]
                sample_rows.append({"candidateMesh": name, "candidateFace": int(face), "sourceInstanceId": source_id, "sourceInstanceIndex": int(source_instances[face]), "sourceFace": int(source_faces[face]), "name": row["name"]})
        del z, az, vertices, faces, source_faces, source_instances, source_ids
    triangles = np.asarray(sample_triangles, dtype=np.float64)
    del sample_triangles
    if len(triangles) > args.max_faces:
        raise ValueError("Deterministic all-instance sample exceeds requested maximum")
    assert len(set(row["sourceInstanceId"] for row in sample_rows)) == 1724
    np.savez_compressed(output / "sample-ownership.npz", triangles_m=triangles, candidate_mesh=np.asarray([r["candidateMesh"] for r in sample_rows]), candidate_face=np.asarray([r["candidateFace"] for r in sample_rows]), source_instance_id=np.asarray([r["sourceInstanceId"] for r in sample_rows]), source_instance_index=np.asarray([r["sourceInstanceIndex"] for r in sample_rows]), source_face=np.asarray([r["sourceFace"] for r in sample_rows]))
    prep_seconds = time.monotonic() - started
    sample = Mesh(triangles.reshape((-1, 3)), np.arange(len(triangles) * 3).reshape((-1, 3)), skin.frame)
    a, b = BVH(sample), BVH(skin)
    b.save(output / "fixed-skin-bvh.npz")
    index_seconds = time.monotonic() - started - prep_seconds
    checkpoint = output / "checkpoint.json"
    def saved(status, **extra):
        checkpoint.write_text(json.dumps({"status": status, "candidateStage": args.candidate_stage, "sampledFaces": len(triangles), "elapsedSeconds": time.monotonic() - started, **extra}, indent=2) + "\n")
    saved("indices-ready")
    tested, returned, rational_replays = 0, 0, 0
    pieces, pair_chunks, literal_chunks = {}, [], []
    complete, stop_reason = False, None
    query_started = time.monotonic()
    for pairs in a.candidate_batches(b, args.distance_m, batch_size=2048):
        if tested + len(pairs) > args.max_tests:
            pairs = pairs[:args.max_tests - tested]
            stop_reason = "finite triangle-test budget"
        if not len(pairs):
            break
        result = triangle_pairs(triangles[pairs[:, 0]], b.triangles[pairs[:, 1]])
        tested += len(pairs)
        near = np.flatnonzero(result["distance_m"] <= args.distance_m)
        if returned + len(near) > args.max_pairs:
            near = near[:args.max_pairs - returned]
            stop_reason = "finite returned-pair budget"
        if len(near):
            kept = pairs[near]
            literal = np.zeros(len(near), dtype=bool)
            for i in np.flatnonzero(result["relation"][near] == "crossing"):
                literal[i] = exact_transverse_crossing(triangles[kept[i, 0]], b.triangles[kept[i, 1]])
                rational_replays += 1
            pair_chunks.append(kept)
            literal_chunks.append(literal)
            for key, value in result.items():
                pieces.setdefault(key, []).append(value[near])
            returned += len(near)
        saved("bounded-triangle-query-running", testedPairs=tested, returnedPairs=returned)
        if stop_reason:
            break
        if time.monotonic() - started >= args.soft_seconds:
            stop_reason = "soft saved-stage wall budget"
            break
    else:
        complete = True
    query_seconds = time.monotonic() - query_started
    pairs = np.concatenate(pair_chunks) if pair_chunks else np.empty((0, 2), dtype=np.int64)
    result = {key: np.concatenate(value) for key, value in pieces.items()}
    literal = np.concatenate(literal_chunks) if literal_chunks else np.empty(0, dtype=bool)
    if not len(pairs):
        result = triangle_pairs(np.empty((0, 3, 3)), np.empty((0, 3, 3)))
    strict = np.flatnonzero(literal)
    actual_source_face = np.asarray([r["sourceFace"] for r in sample_rows], dtype=np.int64)
    actual_instance_index = np.asarray([r["sourceInstanceIndex"] for r in sample_rows], dtype=np.int64)
    np.savez_compressed(output / "actual-triangle-pairs.npz", pairs=pairs, literal_source_strict_crossing=literal, vessel_candidate_face=np.asarray([r["candidateFace"] for r in sample_rows], dtype=np.int64)[pairs[:, 0]], vessel_original_source_face=actual_source_face[pairs[:, 0]], vessel_source_instance_index=actual_instance_index[pairs[:, 0]], skin_original_face=pairs[:, 1], skin_source_part_index=skin_owner[pairs[:, 1]], skin_original_source_face=skin_source_face[pairs[:, 1]], **result)
    strict_rows = []
    exact_segments, exact_errors, segment_pair_ids = [], [], []
    for i in strict:
        si, fi = pairs[i]
        strict_rows.append({**sample_rows[si], "sampleIndex": int(si), "skinFace": int(fi), "skinSourcePartIndex": int(skin_owner[fi]), "skinOriginalSourceFace": int(skin_source_face[fi])})
    # A bounded exact endpoint subset makes crossings independently reviewable.
    for i in strict[:128]:
        si, fi = pairs[i]
        segment = exact_intersection_segment(triangles[si], b.triangles[fi])
        assert segment["relation"] == "crossing"
        exact_segments.append(segment["endpoints"])
        exact_errors.append(segment["coordinate_error_upper_m"])
        segment_pair_ids.append(i)
    np.savez_compressed(output / "exact-crossing-segment-examples.npz", pair_row=np.asarray(segment_pair_ids, dtype=np.int64), endpoints_m=np.asarray(exact_segments, dtype=np.float64).reshape((-1, 2, 3)), conversion_error_upper_m=np.asarray(exact_errors, dtype=np.float64).reshape((-1, 2, 3)))
    counts = Counter(row["name"] for row in strict_rows)
    pair_count = len(pairs)
    near_thresholds = {str(d): {"facePairs": int(np.sum(result["distance_m"] <= d)), "distinctVesselSampleFaces": int(len(np.unique(pairs[result["distance_m"] <= d, 0])))} for d in [.0005, .002, .003]}
    full_count = sum(row["candidateFaceCount"] for row in rows)
    estimated_stages = int(np.ceil(full_count / 32768))
    report = {"schema": "actual-vessel-primary-skin-triangle-benchmark/v1", "candidateStage": args.candidate_stage, "identityFrameBinding": identity_frame_binding, "inputBindings": bindings, "sampling": "19 evenly spaced retained candidate/source face IDs per source instance; includes all 1724 instances; deterministic, deliberately not an unbiased whole-face prevalence estimate", "sampledVesselFaces": len(triangles), "fixedPrimarySkinFaces": len(skin.faces), "skinOwnershipInput": str(fixed / "ownership.json"), "sourceFaceOrderUnchanged": True, "nearDistance_m": args.distance_m, "broadPhase": "Complete Float64 triangle AABB BVH pair traversal expanded by Euclidean threshold; every candidate pair replayed on actual source triangles", "testedCandidatePairs": tested, "returnedNearPairs": pair_count, "completeForSelectedFacesWithinDistance": complete, "stopReason": stop_reason, "maxTests": args.max_tests, "maxPairs": args.max_pairs, "relations": dict(Counter(result["relation"].tolist())), "strictCrossingPairsLiteralSourceConfirmed": len(strict), "strictCrossingVesselFaces": len(np.unique(pairs[strict, 0])), "strictCrossingSkinFaces": len(np.unique(pairs[strict, 1])), "literalRationalReplays": rational_replays, "float64CrossingRejectedByLiteralReplay": int(np.sum((result["relation"] == "crossing") & ~literal)), "nearThresholds": near_thresholds, "strictByVesselName": [{"name": name, "sampledStrictCrossingPairs": count} for name, count in counts.most_common()], "strictWitnessOwnership": strict_rows, "pairRelationQualification": "crossing positives are rationally verified on literal original-coordinate triangles; touching/edge/coplanar and near cases remain separately recorded. Float64 negative predicates, degeneracy and skin incompleteness are not a global containment proof", "actualSkinOnly": True, "classifierCapsUsed": False, "insideOutsideClassified": False, "clearanceAccepted": False, "geometryModified": False, "installed": False, "timing": {"preparationSeconds": prep_seconds, "bvhBuildAndSaveSeconds": index_seconds, "querySeconds": query_seconds, "testedPairsPerQuerySecond": tested / query_seconds if query_seconds else None, "sampleFacesPerQuerySecond": len(triangles) / query_seconds if query_seconds else None, "elapsedSeconds": time.monotonic() - started}, "completeStagedStrategy": {"fullCandidateFaces": full_count, "maximumFacesPerSavedStage": 32768, "minimumDisjointFaceRangeStages": estimated_stages, "partition": "Disjoint contiguous face-ID ranges per candidate mesh, max32768 faces each; persist exact candidate mesh/ancestry hash and fixed skin BVH hash; do not use the per-instance benchmark sample as the full-scan partition", "indexReuse": str(output / "fixed-skin-bvh.npz"), "finiteJobBudget": "40s/1GiB per saved stage; stream pair results; if pair or time cap is reached, preserve exact completed face/subtree cursor, subdivide remaining stage and never label it complete", "strictPass": "0m BVH first for full strict-crossing inventory; separate0.5/2/3mm near-clearance passes under explicit pair budgets", "complement": "Combine actual-triangle crossing witnesses with independent vertex/centroid multi-ray inside and boundary uncertainty; Euclidean triangle distance alone is unsigned"}}
    (output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    assert all(sha(path) == value for path, value in bindings.items()), "Bound candidate or fixed-skin input changed"
    saved("benchmark-complete" if complete else "benchmark-saved-incomplete", testedPairs=tested, returnedPairs=pair_count)
    print(json.dumps({"output": str(output), "sampleFaces": len(triangles), "testedPairs": tested, "nearPairs": pair_count, "strictPairs": len(strict), "relations": report["relations"], "complete": complete, "timing": report["timing"], "topStrictVessels": report["strictByVesselName"][:10]}))


if __name__ == "__main__":
    main()
