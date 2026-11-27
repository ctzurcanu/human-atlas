"""Independent source/pose checks for the user-requested male heart-only fit.

Run through geometry_tools.runtime using the geometry Python environment.
This verifies a staged local installation, not anatomical packing acceptance.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from geometry_tools.common import atomic_json, digest_file
from geometry_tools.geometry_instance import instance_pose
from geometry_tools.model import AtlasReader

STAGE = ROOT / ".local-models/geometry-workbench/heart/male-fixed-vessel-heart-pose-v1"
SOURCE = ROOT / ".local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json"
SOURCE_SHA = "6526654a0e412947121518636d1b8004a934093f8a4cf459afcfa55c80a1b63c"
CANDIDATE_SHA = "a0d7b46e033d145617a8937e1aff489db90b7a2808aa76ec90192f155ad95518"
BUNDLE_SHA = "375bed821531d0bb6c85f5b330362161c9e50dcef2af40f670292844a03e3040"


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["before", "after"])
    parser.add_argument("--catalogue", type=Path, default=SOURCE)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert digest_file(STAGE / "atlas.json") == CANDIDATE_SHA
    assert digest_file(STAGE / "male/cardiac-label-sites.json") == BUNDLE_SHA
    candidate = json.loads((STAGE / "atlas.json").read_text())
    staged_bundle = json.loads((STAGE / "male/cardiac-label-sites.json").read_text())
    proof = json.loads((STAGE / "preview-proof.json").read_text())
    moving = {row["part"] for row in proof["physical_tissues"]}
    assert len(moving) == 33
    candidate_parts = {p["id"]: p for p in candidate["parts"]}
    assert len(candidate_parts) == len(candidate["parts"]) == 6938
    matrix = np.asarray(proof["matrix"], dtype=float)
    scale = float(np.cbrt(np.linalg.det(matrix[:3, :3])))
    rotation = matrix[:3, :3] / scale
    assert np.array_equal(matrix[3], [0, 0, 0, 1])
    assert scale > 0 and np.linalg.det(rotation) > 0
    assert np.max(abs(rotation.T @ rotation - np.eye(3))) < 1e-12
    bundle = staged_bundle["integration"]
    assert bundle["labeledTermCount"] == 113
    assert bundle["referenceTermCount"] == len(bundle["termGroups"]) == len(staged_bundle["concepts"]) == 71
    assert len(bundle["carrierPoses"]) == 22
    for pid, pose in bundle["carrierPoses"].items():
        assert pid in moving
        p = candidate_parts[pid]
        expected = {k: p.get(k, d) for k, d in [("geometryOrigin", [0, 0, 0]), ("geometryInstance", None), ("chunk", None), ("positions", None), ("indices", None), ("indexCount", None)]}
        assert pose == expected, pid
    report = dict(schema="independent-male-fixed-vessel-heart-install-check/v1", mode=args.mode,
                  candidate_sha256=CANDIDATE_SHA, bundle_sha256=BUNDLE_SHA,
                  uniform_scale=scale, physical_tissues=33, total_parts=6938,
                  controls=113, typed_references=71, carrier_poses=22,
                  anatomical_completion=False, packing_acceptance=False,
                  pending_cuff_interfaces=len(candidate["finiteCardiacIntegration"]["pendingRootInterfaces"]))
    assert report["pending_cuff_interfaces"] == 16
    actual = json.loads(args.catalogue.read_text())
    if args.mode == "before":
        assert digest_file(args.catalogue) == SOURCE_SHA
        old = {p["id"]: p for p in actual["parts"]}
        assert old.keys() == candidate_parts.keys()
        changed = {pid for pid in old if old[pid] != candidate_parts[pid]}
        assert changed == moving
        assert actual["chunks"] == candidate["chunks"]
        assert actual["concepts"] == candidate["concepts"]
        allowed = {"geometryOrigin", "geometryInstance", "bounds", "cardiacThoracicFit"}
        reader = AtlasReader()
        errors, bank_hashes = [], {}
        for pid in sorted(moving):
            before, after = old[pid], candidate_parts[pid]
            assert {k: v for k, v in before.items() if k not in allowed} == {k: v for k, v in after.items() if k not in allowed}, pid
            bmesh = reader.read_part(actual, before, frame="model:local-ta98:assembled-metres")
            amesh = reader.read_part(candidate, after, frame="model:local-ta98:assembled-metres")
            assert bmesh.source["bank_sha256"] == amesh.source["bank_sha256"]
            assert np.array_equal(bmesh.faces, amesh.faces)
            error = float(np.max(abs(amesh.vertices - (bmesh.vertices @ matrix[:3, :3].T + matrix[:3, 3]))))
            assert error < 2e-15, (pid, error)
            errors.append(error)
            bank_hashes[amesh.source["bank_sha256"]] = True
            br, bs = instance_pose(before)
            ar, ass = instance_pose(after)
            assert abs(ass / bs - scale) < 1e-14
            assert np.max(abs(ar - rotation @ br)) < 1e-14
        old_bundle = json.loads(args.catalogue.with_name("cardiac-label-sites.json").read_text())
        assert old_bundle["concepts"] == staged_bundle["concepts"]
        for before, after in zip(old_bundle["integration"]["termGroups"], bundle["termGroups"]):
            assert before["geometry"] == after["geometry"]
            assert before["anchor"]["partId"] == after["anchor"]["partId"]
            point = np.asarray(before["anchor"]["point"]) @ matrix[:3, :3].T + matrix[:3, 3]
            assert np.max(abs(point - after["anchor"]["point"])) < 2e-15
            assert before.get("points") is after.get("points") is None
        report.update(source_sha256=SOURCE_SHA, changed_parts=sorted(changed),
                      fixed_part_records_literal=len(old) - len(moving),
                      chunks_literal=True, concepts_literal=True, local_typed_geometry_literal=True,
                      geometry_banks_sha256=sorted(bank_hashes),
                      maximum_world_replay_error_m=max(errors))
    else:
        rows = {p["id"]: p for p in actual["parts"]}
        assert rows == candidate_parts, "Installed part records differ from independently audited heart-only candidate"
        assert actual["chunks"] == candidate["chunks"]
        assert actual["concepts"] == candidate["concepts"]
        actual_bundle_path = args.bundle or args.catalogue.with_name("cardiac-label-sites.json")
        assert json.loads(actual_bundle_path.read_text()) == staged_bundle
        report.update(installed_catalogue_sha256=digest_file(args.catalogue),
                      all_installed_part_records_equal_candidate=True,
                      installed_bundle_equal_candidate=True, chunks_literal=True, concepts_literal=True)
    report["elapsed_seconds"] = time.monotonic() - started
    atomic_json(args.output, report)
    print(json.dumps({k: v for k, v in report.items() if k not in {"changed_parts", "geometry_banks_sha256"}}, sort_keys=True))


if __name__ == "__main__":
    main()
