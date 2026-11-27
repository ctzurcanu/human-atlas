"""Initialize the complete local blood donor in the fixed male skeletal frame.

This creates a reversible, proper similarity only. It does not install or
certify anatomical fit, repair native defects, or assert vessel patency.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.fit import landmark_transform


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT / ".local-models/geometry-workbench/biodigital-male-circulation-v1")
    parser.add_argument("--source-subdir", default="source")
    parser.add_argument("--output-subdir", default="registration-v1")
    args = parser.parse_args()
    started = time.monotonic()
    base = args.root.resolve()
    output = base / args.output_subdir
    if output.exists():
        raise FileExistsError(output)
    landmarks_path = base / "landmarks/landmarks.json"
    inventory_path = base / "ta98-inventory/inventory.json"
    source_dir = base / args.source_subdir
    source_path = source_dir / "manifest.json"
    landmarks = json.loads(landmarks_path.read_text())
    inventory = json.loads(inventory_path.read_text())
    source = json.loads(source_path.read_text())
    target_path = Path(landmarks["target_catalog"])
    target_sha = sha(target_path)
    if target_sha != landmarks["target_catalog_sha256"]:
        raise ValueError("Current male catalogue differs from the measured landmark target")
    if sha(Path(landmarks["source_catalog"])) != landmarks["source_catalog_sha256"]:
        raise ValueError("Donor catalogue changed after measuring landmarks")
    for name, record in source["files"].items():
        if sha(source_dir / name) != record["sha256"]:
            raise ValueError("Imported source bank changed: " + name)
    indices = np.asarray(landmarks["initializer_indices"], dtype=int)
    x = np.asarray(landmarks["source_landmarks_native"], dtype=float)
    y = np.asarray(landmarks["target_landmarks_m"], dtype=float)
    # Freeze an explicit experiment before the numerical call. Scale converts
    # undeclared source coordinates to metres, rather than assuming centimetres.
    contract = dict(
        schema="biodigital-whole-blood-registration-contract/v1",
        mode="proper positive similarity",
        source_manifest_sha256=sha(source_path),
        source_landmarks_sha256=sha(landmarks_path),
        ta98_inventory_sha256=sha(inventory_path),
        target_catalogue_sha256=target_sha,
        moving_scope=f"all {source['counts']['instances']} selected vascular instances as one closure; cardiac anatomy and lymph excluded by explicit scope policy",
        fixed_scope="all target skeleton, muscles, organs and current circulation; no target mutation",
        source_unit_conversion="determined by landmark similarity; native units undeclared",
        positive_scale_bounds_metres_per_source_unit=[0.001, 0.1],
        initializer_indices=indices.tolist(),
        anatomy_gate="Initial correspondence residuals are diagnostic, not an installation acceptance tolerance",
    )
    output.mkdir(parents=True)
    save(output / "contract.json", contract)

    matrix, metrics = landmark_transform(x[indices], y[indices], similarity=True, scale_bounds=(.001, .1))
    scale = metrics["scale"]
    rotation = matrix[:3, :3] / scale
    transformed = x @ matrix[:3, :3].T + matrix[:3, 3]
    residuals = np.linalg.norm(transformed - y, axis=1)
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-12) or np.linalg.det(rotation) <= 0:
        raise ValueError("Registration must be a proper rotation with positive uniform scale")
    rows = [dict(label=row["label"], role=row["proposed_role"],
                 source_native=x[i].tolist(), target_m=y[i].tolist(),
                 fitted_m=transformed[i].tolist(), residual_m=float(residuals[i]))
            for i, row in enumerate(landmarks["landmarks"])]
    result = dict(
        schema="biodigital-whole-blood-similarity/v1",
        status="complete initialization; anatomical fit and installation unverified",
        matrix=matrix.tolist(), matrix_convention="row-major; column homogeneous vector; target_m = matrix @ source_native",
        rotation=rotation.tolist(), scale_metres_per_source_unit=scale,
        metrics=metrics, correspondences=rows,
        holdout_rms_m=float(np.sqrt(np.mean(residuals[np.setdiff1d(np.arange(len(x)), indices)] ** 2))),
        all_instances_share_one_transform=True, target_catalogue_sha256=target_sha,
        target_skeleton_changed=False, original_source_changed=False,
        source_native_defects_preserved=True, anatomical_acceptance=False, installed=False,
        elapsed_seconds=time.monotonic() - started,
    )
    save(output / "transform.json", result)

    # Source matrices remain immutable; this file is only a reversible placement.
    source_instances = json.loads((source_dir / "instances.json").read_text())["instances"]
    fitted = []
    for instance in source_instances:
        original = np.asarray(instance["worldMatrix"]).reshape((4, 4), order="F")
        world = matrix @ original
        fitted.append(dict(id=instance["id"], geometryId=instance["geometryId"],
                           fittedWorldMatrixColumnMajor=world.reshape(-1, order="F").tolist(),
                           reverseWinding=instance["reverseWinding"]))
    save(output / "fitted-instances.json", dict(instances=fitted, transform_sha256=sha(output / "transform.json")))
    if sha(target_path) != target_sha:
        raise ValueError("Target changed during registration")
    print(json.dumps(dict(output=str(output), scale=scale,
                          initializer_rms_mm=1000 * metrics["landmark_rms_m"],
                          holdout_rms_mm=1000 * result["holdout_rms_m"],
                          residual_mm={row["label"]:1000 * row["residual_m"] for row in rows},
                          elapsed_seconds=time.monotonic() - started)))


if __name__ == "__main__":
    main()
