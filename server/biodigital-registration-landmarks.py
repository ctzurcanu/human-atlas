#!/usr/bin/env python3
"""Measure named bone-region correspondences; never alter either atlas.

BioDigital coordinates retain native units. The source does not declare a metric
unit, so a raw-native-to-target-metres similarity must determine scale.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import struct
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from geometry_tools.model import AtlasReader


PAIRS = [
    ("C1", "atlas_first_cervical_vertebra", "ZA:Atlas (C1)", "initializer"),
    ("T1", "first_thoracic_vertebra", "ZA:Vertebra T1", "initializer"),
    ("T7", "seventh_thoracic_vertebra", "ZA:Vertebra T7", "initializer"),
    ("T12", "twelfth_thoracic_vertebra", "ZA:Vertebra T12", "initializer"),
    ("L5", "fifth_lumbar_vertebra", "ZA:Vertebra L5", "initializer"),
    ("sacrum", "sacrum", "ZA:Sacrum", "initializer"),
    ("manubrium", "manubrium_of_sternum", "ZA:Manubrium of sternum", "initializer"),
    ("sternal body", "body_of_sternum", "ZA:Body of sternum", "initializer"),
    ("left clavicle", "left_clavicle", "ZA:Clavicle.l", "initializer"),
    ("right clavicle", "right_clavicle", "ZA:Clavicle.r", "initializer"),
    ("left humerus", "left_humerus", "ZA:Humerus.l", "posture-sensitive holdout"),
    ("right humerus", "right_humerus", "ZA:Humerus.r", "posture-sensitive holdout"),
    ("left femur", "left_femur", "ZA:Femur.l", "posture-sensitive holdout"),
    ("right femur", "right_femur", "ZA:Femur.r", "posture-sensitive holdout"),
]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def resource_id(value, kind):
    if not value:
        return None
    first, remainder = value.split("/", 1)
    return first + "/" + kind + "/" + remainder


def local_matrix(transform):
    """Exact viewer/app.js semantics: Euler ZYX and T(t+p) R S T(-p)."""
    t = np.array([transform.get("translate", {}).get(k, 0) for k in "xyz"], float)
    p = np.array([transform.get("pivot", {}).get(k, 0) for k in "xyz"], float)
    s = np.array([transform.get("scale", {}).get(k, 1) for k in "xyz"], float)
    x, y, z = np.deg2rad([transform.get("rotate", {}).get(k, 0) for k in "xyz"])
    cx, cy, cz = np.cos([x, y, z])
    sx, sy, sz = np.sin([x, y, z])
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    matrix = np.eye(4)
    matrix[:3, :3] = rz @ ry @ rx @ np.diag(s)
    matrix[:3, 3] = t + p - matrix[:3, :3] @ p
    return matrix


def flatten(nodes, records, parent=None, ancestors=()):
    if parent is None:
        parent = np.eye(4)
    for node in nodes:
        tid = resource_id(node.get("transformId"), "transforms")
        transform = dict(records.get(tid, {}))
        transform.update({k: node[k] for k in ("translate", "rotate", "scale", "pivot") if k in node})
        matrix = parent @ local_matrix(transform)
        current = ancestors + ((node["objectId"], tid, transform),)
        if node.get("geometryId"):
            yield node, matrix, current
        yield from flatten(node.get("objects", []), records, matrix, current)


def decode_bda(raw):
    u32 = lambda offset: struct.unpack_from("<I", raw, offset)[0]
    version = u32(0)
    if version == 17:
        start = 44 + u32(4) + u32(8)
        p, n, ix = u32(12), u32(16), u32(20)
        vertices = np.frombuffer(raw, "<f4", p // 4, start).reshape(-1, 3).copy()
        faces = np.frombuffer(raw, "<u4", ix // 4, start + p + n).reshape(-1, 3).copy()
    elif version == 18:
        name_bytes, meta_bytes, header_bytes, index_offset, pos_bytes, normal_bytes, index_bytes = struct.unpack_from("<7I", raw, 4)
        decode_bytes, index_type, uv_sets = struct.unpack_from("<HBB", raw, 32)
        offset = header_bytes + name_bytes + meta_bytes
        decoder = np.frombuffer(raw, "<f4", decode_bytes // 4, offset).copy() if decode_bytes else None
        offset += decode_bytes + uv_sets * 16
        quantized = np.frombuffer(raw, "<u2", pos_bytes // 2, offset).reshape(-1, 3)
        # Match browser Float32Array assignment after double multiply/add.
        vertices = ((quantized.astype(float) * decoder[3:6] + decoder[:3]).astype("<f4")
                    if decoder is not None else quantized.astype("<f4"))
        if index_type not in (1, 2, 4):
            raise ValueError("invalid index type")
        faces = np.frombuffer(raw, {1: "u1", 2: "<u2", 4: "<u4"}[index_type], index_bytes // index_type, index_offset).reshape(-1, 3).copy()
    else:
        raise ValueError(f"unsupported BDA version {version}")
    if faces.max() >= len(vertices):
        raise ValueError("invalid source triangles; do not silently drop them for landmarks")
    return vertices, faces, version


def measure(vertices, faces):
    triangles = vertices[faces]
    twice_area = np.linalg.norm(np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]), axis=1)
    centroid = np.sum(triangles.mean(axis=1) * twice_area[:, None], axis=0) / twice_area.sum()
    used = np.unique(faces)
    return {
        "area_weighted_surface_centroid": centroid.tolist(),
        "used_vertex_centroid": vertices[used].mean(axis=0).tolist(),
        "bounds": [vertices[used].min(axis=0).tolist(), vertices[used].max(axis=0).tolist()],
        "surface_area_coordinate_units_squared": float(twice_area.sum() / 2),
        "used_vertex_count": len(used),
        "triangle_count": len(faces),
        "zero_area_triangle_count": int(np.count_nonzero(twice_area == 0)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-catalog", type=Path, default=Path("/Users/christiantzurcanu/Documents/dev/fetch/viewer/data/male/catalog.json"))
    parser.add_argument("--target-catalog", type=Path, default=ROOT / ".local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json")
    parser.add_argument("--output", type=Path, default=ROOT / ".local-models/geometry-workbench/biodigital-male-circulation-v1/landmarks")
    args = parser.parse_args()
    start = time.monotonic()
    if args.output.exists():
        raise ValueError("evidence output exists; use a fresh output to preserve previous measurements")
    source = json.loads(args.source_catalog.read_text())
    target = json.loads(args.target_catalog.read_text())
    entries = {n["objectId"]: (n, m, ancestors) for n, m, ancestors in flatten(source["sceneGraph"], source["records"])}
    reader = AtlasReader(root=ROOT)
    hidden = set(target.get("defaultHidden", []))
    evidence = []
    selected_arrays = {}
    for label, suffix, target_concept, role in PAIRS:
        source_id = "human_20_male_skeletal_system-" + suffix + "_ID"
        node, matrix, ancestors = entries[source_id]
        geometry_key = resource_id(node["geometryId"], "geometry")
        geometry = source["geometry"][geometry_key]
        pack = args.source_catalog.parent / geometry["file"]
        with pack.open("rb") as stream:
            stream.seek(geometry["offset"])
            raw = stream.read(geometry["length"])
        vertices, faces, version = decode_bda(raw)
        world = vertices.astype(float) @ matrix[:3, :3].T + matrix[:3, 3]
        # Reflections are present in source instancing; positions need only the
        # literal world matrix. Area centroids do not depend on winding.
        source_measure = measure(world, faces)
        target_parts = [p for p in target["parts"] if p.get("conceptId") == target_concept and p["id"] not in hidden and not p.get("suppressed")]
        if not target_parts or any(not p["id"].startswith("ZA:") or p["system"] != "skeletal" for p in target_parts):
            raise ValueError("target must be the active preferred original ZA skeletal assembly")
        target_v, target_f, offset, target_controls = [], [], 0, []
        for part in target_parts:
            mesh = reader.read_part(target, part, frame="model:local-ta98:assembled-metres", catalogue_path=args.target_catalog)
            target_v.append(mesh.vertices)
            target_f.append(mesh.faces + offset)
            offset += len(mesh.vertices)
            target_controls.append({"part_id": part["id"], "bank_path": mesh.source["bank_path"], "bank_sha256": mesh.source["bank_sha256"], "geometry_origin": mesh.source["geometry_origin"], "source_offset": mesh.source["source_offset"], "geometry_scale": mesh.source["geometry_scale"], "default_hidden": False})
        target_vertices, target_faces = np.concatenate(target_v), np.concatenate(target_f)
        target_measure = measure(target_vertices, target_faces)
        evidence.append({
            "label": label, "proposed_role": role,
            "evidence_type": "whole named bone region area-weighted exterior surface centroid; not a joint center or vertebral body center",
            "source": {"object_id": source_id, "name": node["displayName"], "geometry_key": geometry_key, "bda_version": version, "bda_slice": geometry, "bda_slice_sha256": hashlib.sha256(raw).hexdigest(), "world_matrix_row_major": matrix.tolist(), "world_matrix_column_major": matrix.T.ravel().tolist(), "world_matrix_determinant": float(np.linalg.det(matrix[:3, :3])), "nonidentity_transform_ancestry": [{"object_id": oid, "transform_record_key": tid, "transform": transform} for oid, tid, transform in ancestors if not np.array_equal(local_matrix(transform), np.eye(4))], "coordinate_frame": "original BioDigital scene world; hierarchy matrices applied exactly once; native units undeclared", **source_measure},
            "target": {"concept_id": target_concept, "parts": target_controls, "coordinate_frame": "current local-ta98 assembled metres; literal source Float32 * geometryScale + geometryOrigin + sourceOffset exactly once", **target_measure},
        })
        key = suffix
        selected_arrays[key + "_source_world_native"] = world
        selected_arrays[key + "_source_faces"] = faces
        selected_arrays[key + "_target_world_m"] = target_vertices
        selected_arrays[key + "_target_faces"] = target_faces
    source_points = np.array([e["source"]["area_weighted_surface_centroid"] for e in evidence])
    target_points = np.array([e["target"]["area_weighted_surface_centroid"] for e in evidence])
    initializer = np.array([e["proposed_role"] == "initializer" for e in evidence])
    report = {
        "schema": "biodigital-named-bone-region-correspondence-v1", "status": "measured initialization evidence; anatomical registration unverified; no installation",
        "source_catalog": str(args.source_catalog.resolve()), "source_catalog_sha256": sha(args.source_catalog),
        "target_catalog": str(args.target_catalog.resolve()), "target_catalog_sha256": sha(args.target_catalog),
        "source_unit_evidence": "No explicit metric unit declaration found in source catalog, viewer/app.js, viewer/README.md or male-complete/README.md. Native extent resembles centimeters but that interpretation is unverified. Fit raw source-native coordinates directly to target metres; do not apply a guessed height or unit conversion.",
        "target_height_evidence": target.get("heightAdjustment"),
        "method": "Area-weighted triangle centroid of corresponding complete named source bone and all active original ZA material patches for that named bone. Vertex and bounds diagnostics are retained. Whole-bone regional centroids are reproducible coarse registration landmarks; shape/material partition differences, source pose and sampling can bias them.",
        "frame_controls": "Hidden supplementary FJ/LOCAL skeletal variants excluded. All selected target skeletal records are default-visible ZA original model parts. No sourceOffset/origin applied twice. No geometry altered.",
        "landmarks": evidence,
        "source_landmarks_native": source_points.tolist(), "target_landmarks_m": target_points.tolist(),
        "initializer_indices": np.flatnonzero(initializer).tolist(),
        "initializer_source_centered_singular_values_native": np.linalg.svd(source_points[initializer] - source_points[initializer].mean(axis=0), compute_uv=False).tolist(),
        "initializer_target_centered_singular_values_m": np.linalg.svd(target_points[initializer] - target_points[initializer].mean(axis=0), compute_uv=False).tolist(),
        "elapsed_seconds": time.monotonic() - start,
    }
    args.output.mkdir(parents=True)
    np.savez(args.output / "selected-bone-regions.npz", **selected_arrays)
    (args.output / "landmarks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"output": str(args.output), "landmarks": len(evidence), "initializer_landmarks": int(initializer.sum()), "target_sha256": report["target_catalog_sha256"], "elapsed_seconds": report["elapsed_seconds"]}))
    for e in evidence:
        print(e["label"], "source_native", np.round(e["source"]["area_weighted_surface_centroid"], 6).tolist(), "target_m", np.round(e["target"]["area_weighted_surface_centroid"], 6).tolist())


if __name__ == "__main__":
    main()
