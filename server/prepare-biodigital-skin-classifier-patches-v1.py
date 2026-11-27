#!/usr/bin/env python3
"""Make bounded, source-bound classification patches; never rendered anatomy."""
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def cross(a, b, c):
    q, r = b - a, c - a
    return float(q[0] * r[1] - q[1] * r[0])


def triangulate(ids, vertices):
    """Ear clip the boundary's PCA projection without adding a centre vertex."""
    xyz = vertices[ids]
    _, _, vt = np.linalg.svd(xyz - xyz.mean(axis=0), full_matrices=False)
    uv = (xyz - xyz.mean(axis=0)) @ vt[:2].T
    area2 = sum(cross(np.zeros(2), uv[i], uv[(i + 1) % len(ids)]) for i in range(len(ids)))
    scale2 = max(float(np.ptp(uv, axis=0).max()) ** 2, 1e-20)
    eps = scale2 * 1e-12
    if abs(area2) <= eps:
        raise ValueError("degenerate projected cycle")
    orientation = 1 if area2 > 0 else -1
    n = len(ids)
    for i in range(n):
        a, b = uv[i], uv[(i + 1) % n]
        for j in range(i + 1, n):
            if j in (i, (i + 1) % n) or (j + 1) % n in (i, (i + 1) % n):
                continue
            c, d = uv[j], uv[(j + 1) % n]
            if cross(a, b, c) * cross(a, b, d) < -eps * eps and cross(c, d, a) * cross(c, d, b) < -eps * eps:
                raise ValueError("self-crossing projected cycle")
    left, result = list(range(n)), []
    while len(left) > 3:
        found = False
        for k, j in enumerate(left):
            i, l = left[k - 1], left[(k + 1) % len(left)]
            if orientation * cross(uv[i], uv[j], uv[l]) <= eps:
                continue
            others = [q for q in left if q not in (i, j, l)]
            if any(min(orientation * cross(uv[i], uv[j], uv[q]), orientation * cross(uv[j], uv[l], uv[q]), orientation * cross(uv[l], uv[i], uv[q])) >= -eps for q in others):
                continue
            result.append([ids[i], ids[j], ids[l]])
            left.pop(k)
            found = True
            break
        if not found:
            raise ValueError("ear clipping did not resolve cycle")
    if orientation * cross(uv[left[0]], uv[left[1]], uv[left[2]]) <= eps:
        raise ValueError("degenerate final projected triangle")
    result.append([ids[q] for q in left])
    result = np.asarray(result, dtype=np.int64)
    tri = vertices[result]
    areas = np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1) / 2
    if np.any(areas <= 1e-18):
        raise ValueError("zero area 3D patch triangle")
    return result, {"maximumPlaneResidual_m": float(np.abs((xyz - xyz.mean(axis=0)) @ vt[2]).max()), "area_m2": float(areas.sum()), "boundaryVertexCount": n}


def main():
    started = time.monotonic()
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixed-skin", type=Path, default=Path(".local-models/geometry-workbench/biodigital-male-circulation-v1/global-envelope-v1/fixed-skin"))
    args = parser.parse_args()
    base = args.fixed_skin
    topology = base / "envelope-topology-v1"
    output = topology / "classifier-patches-v1"
    output.mkdir(parents=True, exist_ok=True)
    paths = [base / "skin.npz", base / "policy.json", base / "ownership.json", topology / "proxy-and-source-evidence.npz", topology / "geometric-seam-coverage-v1/report.json", topology / "geometric-seam-coverage-v1/interval-witnesses.npz"]
    bindings = {str(p): digest(p) for p in paths}
    report = json.loads(paths[4].read_text())
    evidence = np.load(paths[3])
    skin = np.load(paths[0])
    vertices = evidence["tol_1um_proxy_vertices_m"]
    faces = evidence["tol_1um_proxy_faces"]
    boundary = evidence["tol_1um_boundary_edges"]
    boundary_faces = evidence["tol_1um_boundary_original_skin_face_ids"]
    assert len(faces) == 45248 and len(boundary) == 1042
    excluded = {32: "paired auricle attachment cycles need a bridge qualification", 33: "paired auricle attachment cycles need a bridge qualification", 40: "urogenital exterior scope unresolved; do not cap off genital attachment", 51: "paired auricle attachment cycles need a bridge qualification", 52: "paired auricle attachment cycles need a bridge qualification"}
    cap_faces, cap_components, cap_records, scopes = [], [], [], []
    for component in report["components"]:
        ci = component["componentIndex"]
        ids = evidence[f"tol_1um_boundary_component_{ci}_proxy_vertex_ids"]
        member = np.zeros(len(vertices), dtype=bool)
        member[ids] = True
        edge_ids = np.flatnonzero(member[boundary].all(axis=1))
        cov = component["coverage"]["other_primary"]
        item = {"componentIndex": ci, "sourcePartIds": [p["id"] for p in component["sourceParts"]], "bounds_m": component["bounds_m"], "simpleBoundaryCycle": component["simpleBoundaryCycle"], "otherPrimaryUncoveredLength_m": cov["uncoveredLength_m"], "sourceBoundaryEdgeIds": edge_ids.tolist()}
        scopes.append(item)
        if cov["uncoveredLength_m"] <= 1e-5:
            item["status"] = "geometrically covered; no patch"
            continue
        if not component["simpleBoundaryCycle"]:
            item["status"] = "unknown; branched uncovered boundary"
            continue
        if ci in excluded:
            item["status"] = "unknown; " + excluded[ci]
            continue
        if len(ids) > 128:
            item["status"] = "unknown; cycle exceeds fixed 128 vertex budget"
            continue
        try:
            following = {}
            for ei in edge_ids:
                u, v = boundary[ei]
                tri = faces[boundary_faces[ei]]
                direct = [(int(tri[k]), int(tri[(k + 1) % 3])) for k in range(3) if set((int(tri[k]), int(tri[(k + 1) % 3]))) == set((int(u), int(v)))]
                if len(direct) != 1 or direct[0][0] in following:
                    raise ValueError("source boundary direction is not one cycle")
                following[direct[0][0]] = direct[0][1]
            cycle = [min(following)]
            while following[cycle[-1]] != cycle[0]:
                cycle.append(following[cycle[-1]])
                if len(cycle) > len(ids):
                    raise ValueError("cycle walk exceeded fixed source vertex count")
            if len(cycle) != len(ids):
                raise ValueError("cycle walk omitted source vertices")
            # Opposite directed boundary closes the source-face orientation.
            cycle = list(reversed(cycle))
            tris, detail = triangulate(cycle, vertices)
            begin = len(cap_faces)
            cap_faces.extend(tris.tolist())
            cap_components.extend([ci] * len(tris))
            item.update({"status": "provisional classifier patch; unknown anatomy", "capFaceRange": [begin, len(cap_faces)], **detail})
            cap_records.append({"componentIndex": ci, "orderedProxyVertexIds": cycle, "sourceBoundaryEdgeIds": edge_ids.tolist(), "oppositeSourceWinding": True, **detail})
        except (ValueError, KeyError, np.linalg.LinAlgError) as error:
            item["status"] = "unknown; patch rejected: " + str(error)
    cap_faces = np.asarray(cap_faces, dtype=np.int64).reshape((-1, 3))
    cap_components = np.asarray(cap_components, dtype=np.int32)
    exact_rep = evidence["exact_vertex_representative_original_skin_vertex"]
    representative = exact_rep[evidence["tol_1um_representative_exact_vertex"]]
    np.savez_compressed(output / "classifier-patches.npz", vertices=vertices, faces=cap_faces, cap_component_index=cap_components, vertex_representative_original_skin_vertex=representative, vertex_source_part_index=skin["vertex_owner_part_index"][representative], vertex_source_vertex_id=skin["vertex_owner_source_vertex_id"][representative])
    np.savez_compressed(output / "source-plus-patches.npz", vertices=vertices, faces=np.vstack((faces, cap_faces)), face_kind=np.concatenate((np.zeros(len(faces), dtype=np.uint8), np.ones(len(cap_faces), dtype=np.uint8))), face_original_skin_face_id=np.concatenate((np.arange(len(faces), dtype=np.int64), np.full(len(cap_faces), -1, dtype=np.int64))), cap_component_index=np.concatenate((np.full(len(faces), -1, dtype=np.int32), cap_components)), source_face_owner_part_index=skin["face_owner_part_index"], source_face_owner_source_face_id=skin["face_owner_source_face_id"])
    proposal = {"schema": "source-bound-skin-classifier-patches/v1", "status": "minimal provisional patches prepared; globally closed proxy is not certified", "inputBindings": bindings, "sourceFaceCount": len(faces), "sourceTriangleOrderUnchanged": True, "classifierOnlySourceSnapTolerance_m": 1e-6, "renderedSourceUnchanged": True, "newRenderedTissue": False, "fullyCoveredHandNetworksPatched": False, "cyclePatches": cap_records, "componentScopes": scopes, "capFaceCount": len(cap_faces), "triangulation": "Single PCA projection and bounded ear clipping; source boundary vertices only; no centre, relaxation, search, or source repair", "validationLimits": {"readWallSeconds": 20, "numericalWallSeconds": 40, "memoryMiB": 1024, "maximumCycleVertices": 128}, "classifierContract": {"distanceTarget": str(base / "skin.npz"), "patchesMustNeverCountAsObservedTissueOrClearance": True, "nearNormalSignIsNotAnInsideTest": True, "rayUse": "All-intersection rays on literal source skin, geometric hit deduplication, oriented crossings and multi-ray agreement; use these patch triangles only as closure/uncertainty witnesses unless a separate full certification succeeds", "uncertainty": ["Any tested ray crossing a synthetic patch is closure dependent and is excluded from an unqualified inside/outside decision", "Keep open-boundary shadows for every unpatched or rejected component, including branched networks and urogenital component 40", "If all usable directions are excluded, disagree, graze a boundary, or have unresolved duplicate/coplanar crossings, return unknown", "Boundary-loop AABBs are useful local warning regions; AABB separation alone is not a proof that an arbitrary ray classification is independent of closure", "No global watertightness, self-intersection, nesting, or material-discard validation is claimed"]}, "recommendedNextStage": "Use the existing 8192 witness multi-ray benchmark now; no whole-skin topology rebuild. Compare raw-skin results with patch-assisted results, preserving unknown labels and patch-hit IDs.", "elapsedSeconds": time.monotonic() - started}
    (output / "report.json").write_text(json.dumps(proposal, indent=2) + "\n")
    after = {str(p): digest(p) for p in paths}
    assert after == bindings, "bound source evidence changed during patch preparation"
    summary = {"patchCycles": len(cap_records), "capFaces": len(cap_faces), "sourceFaces": len(faces), "crownPatched": any(c["componentIndex"] == 41 for c in cap_records), "unpatchedUncoveredComponents": [c["componentIndex"] for c in scopes if c["status"].startswith("unknown")], "sourceInputsUnchanged": True, "globallyClosedCertified": False, "elapsedSeconds": time.monotonic() - started}
    (output / "validation.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
