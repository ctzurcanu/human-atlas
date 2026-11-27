"""Partition a completed superior-gluteal wall at witnessed paired source rings.

Run with the workspace runtime supervisor.  ``--stage prepare`` saves the
inexpensive, source-bound arrangement; ``--stage pack`` validates and packs that
saved arrangement.  No catalogue, installation, or anatomical acceptance is
performed here.  The whole blood classifier is never partitioned or capped.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from geometry_tools.common import ROOT, atomic_json, atomic_npz, digest_file
from geometry_tools.model import Mesh
from geometry_tools.algorithms.finite_tissue_boolean import oriented_key
from geometry_tools.algorithms.tube_continuity import annular_rim, boundary_loops
from geometry_tools.tissue import normal_words, pack


NAMES = ("trunk", "superficial", "deep", "superior", "inferior")
TERMS = {name: f"A12.2.15.{number:03}" for name, number in zip(NAMES, range(13, 18))}
PART_IDS = {name: "TA98GLUTEAL:male:left:" + name for name in NAMES}
PART_IDS["trunk"] = "ZA:Superior gluteal artery.l"
# The first end is the anatomical upstream end of each native component.
CUTS = ((4, 1, "trunk", "superficial"), (5, 0, "trunk", "deep"),
        (6, 0, "deep", "superior"), (7, 0, "deep", "inferior"))
PORT_OWNERS = {(0, 0): "trunk", (1, 1): "trunk", (2, 0): "superficial",
               (3, 1): "superficial", (6, 1): "superior", (7, 1): "inferior"}
ALLOWED_ANCESTRY = {"trunk": {0, 1, 4, 5, 8}, "superficial": {2, 3, 4, 9},
                    "deep": {5, 6, 7, 10}, "superior": {6}, "inferior": {7}}
EXPECTED_GENUS = dict(zip(NAMES, (3, 2, 2, 1, 1)))


def local_path(value, relative_to=ROOT):
    path = Path(value)
    return (path if path.is_absolute() else relative_to / path).resolve()


def edge_incidence(faces):
    edges = defaultdict(list)
    for fi, face in enumerate(faces):
        for a, b in zip(face, np.roll(face, -1)):
            a, b = int(a), int(b)
            edges[tuple(sorted((a, b)))].append((fi, a, b))
    return edges


def closed_edge_check(faces):
    edges = edge_incidence(faces)
    if not all(len(rows) == 2 and rows[0][1:] == rows[1][1:][::-1]
               for rows in edges.values()):
        raise ValueError("The actual wall must be closed with opposite directed edge uses")
    used = np.unique(faces)
    euler = len(used) - len(edges) + len(faces)
    return edges, dict(vertices=len(used), edges=len(edges), faces=len(faces),
                       euler=int(euler), genus=(2 - euler) / 2)


def oriented_counter(faces):
    return Counter(oriented_key(face) for face in faces)


def native_corner_charts(positions, encoded_positions, surfaces, owner_masks,
                         source_components, source_faces, sources):
    """Keep actual parent charts, including explicitly transported v2 sources."""
    words = normal_words(encoded_positions, surfaces)
    uvs = {name: np.zeros((len(faces), 3, 2), dtype="<f4")
           for name, faces in surfaces.items()}
    replay = {}
    counts = {str(i): dict(literal_corners=0, barycentric_corners=0) for i in range(8)}
    for name in NAMES:
        sid = name + ":outer"
        comp = source_components[owner_masks[name]]
        parent_faces = source_faces[owner_masks[name]]
        barycentrics = np.full((len(comp), 3, 3), np.nan)
        for fi, (face, component, source_face) in enumerate(zip(surfaces[sid], comp, parent_faces)):
            component, source_face = int(component), int(source_face)
            if component >= 8:
                continue  # Only the declared, newly generated fork envelopes.
            source = sources[component]
            if source_face < 0 or source_face >= len(source.faces):
                raise ValueError(("Invalid source face ancestry", component, source_face))
            ids = source.faces[source_face]
            triangle = source.vertices[ids]
            edge = triangle[1:] - triangle[0]
            gram = edge @ edge.T
            face_normal = np.cross(encoded_positions[face[1]] - encoded_positions[face[0]],
                                   encoded_positions[face[2]] - encoded_positions[face[0]])
            for corner, point in enumerate(positions[face]):
                same = np.flatnonzero(np.all(triangle == point, axis=1))
                if len(same):
                    weights = np.eye(3)[same[0]]
                    normal = source.attributes["normals"][ids[int(same[0])]].copy()
                    counts[str(component)]["literal_corners"] += 1
                else:
                    xy = np.linalg.solve(gram, edge @ (point - triangle[0]))
                    weights = np.r_[1 - xy.sum(), xy]
                    residual = np.linalg.norm(weights @ triangle - point)
                    if weights.min() < -1e-8 or weights.max() > 1 + 1e-8 or residual > 1e-11:
                        raise ValueError(("Changed corner lacks an explicit source chart", name,
                                          component, source_face, corner, float(residual)))
                    value = weights @ source.attributes["normals"][ids].astype(float)
                    length = np.linalg.norm(value)
                    if not length:
                        raise ValueError("Zero interpolated native normal")
                    normal = np.rint(value / length * 32767).astype("<i2")
                    counts[str(component)]["barycentric_corners"] += 1
                if normal @ face_normal <= 0:
                    raise ValueError(("Nonoutward native chart after position encoding", name,
                                      component, source_face, corner))
                words[sid][fi, corner] = normal
                uvs[sid][fi, corner] = weights @ source.attributes["uvs"][ids]
                barycentrics[fi, corner] = weights
        replay[name + "_component"] = comp
        replay[name + "_source_face"] = parent_faces
        replay[name + "_barycentrics"] = barycentrics
    return words, uvs, replay, counts


def prepare(source_directory, output_directory, origin_override):
    began = time.monotonic()
    source = source_directory
    output = output_directory
    if output.exists():
        raise FileExistsError(output)
    plan_path = source / "plan.json"
    plan = json.loads(plan_path.read_text())
    inventory = local_path(plan["source"])
    base = inventory.parent
    paths = [local_path(p) for p in plan.get(
        "source_component_paths", [str(base / f"component-{i:02}.npz") for i in range(8)])]
    if len(paths) != 8:
        raise ValueError("Exactly eight explicit source component charts are required")
    profiles_path = local_path(plan.get("source_profiles_path", str(base / "literal-source-course-graphs-v1.npz")))
    if digest_file(inventory) != plan["source_sha256"]:
        raise ValueError("The declared source inventory changed")
    if "source_component_paths" in plan:
        expected_components = plan.get("source_component_sha256")
        if not isinstance(expected_components, list) or len(expected_components) != 8:
            raise ValueError("Eight declared construction component hashes are required")
        for i, (path, expected) in enumerate(zip(paths, expected_components)):
            if digest_file(path) != expected:
                raise ValueError(("Construction component differs from the declared whole-wall plan", i, str(path)))
    if "source_profiles_path" in plan:
        if digest_file(profiles_path) != plan.get("source_profiles_sha256"):
            raise ValueError("Construction profiles differ from the declared whole-wall plan")
    if "inlet_correction_evidence" in plan:
        inlet_evidence = local_path(plan["inlet_correction_evidence"])
        if digest_file(inlet_evidence) != plan.get("inlet_correction_evidence_sha256"):
            raise ValueError("Declared inlet correction evidence changed")
    sources = [Mesh.load(path) for path in paths]
    for i, mesh in enumerate(sources):
        if set(("normals", "uvs")) - mesh.attributes.keys():
            raise ValueError(("Native source lacks its chart", i))
        if mesh.attributes["normals"].dtype != np.dtype("int16"):
            raise ValueError(("Literal integer native normal words required", i))
        if not np.isfinite(mesh.attributes["uvs"]).all():
            raise ValueError(("Nonfinite source UV", i))
    whole = Mesh.load(source / "finite-wall.npz")
    if any(mesh.frame != whole.frame for mesh in sources):
        raise ValueError("All sources and the completed wall must use the same explicit frame")
    proof = json.loads((source / "evidence.json").read_text())
    if not proof.get("passed") or proof.get("geometry_sha256") != whole.geometry_sha:
        raise ValueError("A completed source-bound whole-wall candidate is required")
    data = np.load(source / "surfaces-and-ancestry.npz", allow_pickle=False)
    profiles = np.load(profiles_path, allow_pickle=False)
    positions = data["vertices"]
    outer, inner = data["outer_faces"], data["inner_faces"]
    port_data = json.loads((source / "ports.json").read_text())["ports"]
    if {(p["component"], p["end"]) for p in port_data} != set(PORT_OWNERS) or len(port_data) != 6:
        raise ValueError("The source family must have exactly its six witnessed external mouths")
    annuli = [data[f"annulus_{i}" if f"annulus_{i}" in data else f"annulus_{i}_faces"] for i in range(6)]
    face_sets = [outer, inner] + annuli
    full = np.concatenate(face_sets)
    if not np.array_equal(positions, whole.vertices) or not np.array_equal(full, whole.faces):
        raise ValueError("Ancestry arrangement must replay the actual complete wall literally")
    edges, whole_edges = closed_edge_check(full)
    used = np.unique(full)
    lookup = {tuple(positions[i]): int(i) for i in used}
    if len(lookup) != len(used):
        raise ValueError("Partitioning requires the whole wall's canonical shared controls")
    fraction = float(plan["inner_radial_fraction"])
    if not 0 < fraction < 1:
        raise ValueError("Finite positive source wall depth required")
    cuts = {}
    blocked = set()
    for component, upstream_end, parent, child in CUTS:
        p = profiles[f"c{component}_profiles"]
        if p.ndim != 3 or p.shape[1:] != (24, 3):
            raise ValueError(("Twenty-four corresponding actual profile controls required", component))
        centers = p.mean(1)
        arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(centers, axis=0), axis=1))]
        candidates = range(len(p) - 2, 0, -1) if upstream_end else range(1, len(p) - 1)
        found = None
        for k in candidates:
            polygons = np.array([p[k], centers[k] + fraction * (p[k] - centers[k])])
            try:
                ids = np.array([[lookup[tuple(point)] for point in polygon] for polygon in polygons])
            except KeyError:
                continue
            edge_set = {tuple(sorted((int(a), int(b))))
                        for polygon in ids for a, b in zip(polygon, np.roll(polygon, -1))}
            if len(edge_set) == 48 and all(e in edges for e in edge_set):
                found = k, ids, edge_set
                break
        if found is None:
            raise ValueError(("No intact paired source ring beyond the declared junction", component))
        k, ids, edge_set = found
        if blocked & edge_set:
            raise ValueError("Distinct ownership cuts must not reuse a source edge")
        blocked.update(edge_set)
        cuts[child] = dict(component=component, source_station=int(k), upstream_end=upstream_end,
                           distance_from_upstream_m=float(arc[-1] - arc[k] if upstream_end else arc[k]),
                           outer=ids[0], inner=ids[1], parent=parent, child=child)
    rows, cols = [], []
    for edge, uses in edges.items():
        if edge not in blocked:
            a, b = uses[0][0], uses[1][0]
            rows.extend((a, b))
            cols.extend((b, a))
    number, labels = connected_components(coo_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(full), len(full))).tocsr(), directed=False)
    if number != 6:
        raise ValueError(("Expected four externally rimmed owners and two deep wall sheets", int(number)))
    offsets = np.cumsum([0] + [len(faces) for faces in face_sets])
    label_owner = {}
    port_owners = []
    for i, port in enumerate(port_data):
        owner = PORT_OWNERS[(port["component"], port["end"])]
        region = np.unique(labels[offsets[2 + i]:offsets[3 + i]])
        if len(region) != 1:
            raise ValueError("An external thickness annulus crosses an ownership boundary")
        region = int(region[0])
        if region in label_owner and label_owner[region] != owner:
            raise ValueError("Unrelated branch mouths remain in the same tissue region")
        label_owner[region] = owner
        port_owners.append(owner)
    if len(label_owner) != 4 or set(label_owner.values()) != {"trunk", "superficial", "superior", "inferior"}:
        raise ValueError("External mouth ownership does not establish four distinct wall regions")
    remaining = set(range(number)) - set(label_owner)
    if len(remaining) != 2:
        raise ValueError("The deep common stem requires its actual inner and outer sheets")
    label_owner.update({region: "deep" for region in remaining})
    owner_index = {name: i for i, name in enumerate(NAMES)}
    face_owner = np.array([owner_index[label_owner[int(region)]] for region in labels])
    surfaces, cells, masks = {}, {name: [] for name in NAMES}, {}
    for name in NAMES:
        masks[name] = face_owner[:len(outer)] == owner_index[name]
        for sheet, faces, start, component_key in (
                ("outer", outer, 0, "outer_component"),
                ("inner", inner, len(outer), "inner_component")):
            mask = face_owner[start:start + len(faces)] == owner_index[name]
            ancestry = set(map(int, np.unique(data[component_key][mask])))
            if not ancestry <= ALLOWED_ANCESTRY[name]:
                raise ValueError(("Connectivity ownership conflicts with the declared source graph", name, sheet, ancestry))
            sid = name + ":" + sheet
            surfaces[sid] = faces[mask]
            if not len(surfaces[sid]):
                raise ValueError(("An owner lacks a finite wall sheet", sid))
            cells[name].append(dict(surface=sid, side=1, neighbor=None))
    for i, owner in enumerate(port_owners):
        sid = "external-mouth:" + str(i)
        surfaces[sid] = annuli[i]
        cells[owner].append(dict(surface=sid, side=1, neighbor=None))
    interfaces = []
    for child, cut in cuts.items():
        parent = cut["parent"]
        outer_loop = next((loop for loop in boundary_loops(surfaces[child + ":outer"])
                           if set(loop) == set(cut["outer"])), None)
        inner_loop = next((loop for loop in boundary_loops(surfaces[child + ":inner"])
                           if set(loop) == set(cut["inner"])), None)
        if outer_loop is None or inner_loop is None:
            raise ValueError(("The exact paired source cut is not the child's boundary", child))
        permutation = [int(np.flatnonzero(cut["outer"] == index)[0]) for index in outer_loop]
        sid = "shared-wall-ring:" + child
        surfaces[sid] = annular_rim(outer_loop, cut["inner"][permutation])
        cells[child].append(dict(surface=sid, side=1, neighbor=parent))
        cells[parent].append(dict(surface=sid, side=-1, neighbor=child))
        interfaces.append({key: value for key, value in cut.items() if key not in ("outer", "inner")}
                          | dict(surface=sid, blood_lumen_open=True, reciprocal_thickness_annulus=True))
    cell_topology = {}
    combined = Counter()
    for name, uses in cells.items():
        faces = np.concatenate([surfaces[u["surface"]] if u["side"] == 1
                                else surfaces[u["surface"]][:, (0, 2, 1)] for u in uses])
        _, topology = closed_edge_check(faces)
        if topology["genus"] != EXPECTED_GENUS[name]:
            raise ValueError(("Unexpected finite owner topology", name, topology))
        cell_topology[name] = topology
        combined.update(oriented_counter(faces))
    cancelled = Counter()
    for key, count in combined.items():
        opposite = oriented_key((key[0], key[2], key[1]))
        if count > combined[opposite]:
            cancelled[key] = count - combined[opposite]
    if cancelled != oriented_counter(full):
        raise ValueError("Reciprocal ownership faces do not replay the original whole wall exactly")
    origin = np.array(origin_override, float) if origin_override is not None else whole.bounds.mean(0)
    encoded = (positions - origin).astype("<f4").astype(float) + origin
    words, uvs, chart_replay, chart_counts = native_corner_charts(
        positions, encoded, surfaces, masks, data["outer_component"], data["outer_source_face"], sources)
    arrays = {"positions": positions}
    surface_map, word_map, uv_map = {}, {}, {}
    for i, (sid, faces) in enumerate(surfaces.items()):
        surface_map[sid], word_map[sid], uv_map[sid] = f"faces_{i}", f"words_{i}", f"uvs_{i}"
        arrays[f"faces_{i}"], arrays[f"words_{i}"], arrays[f"uvs_{i}"] = faces, words[sid], uvs[sid]
    output.mkdir(parents=True)
    arrangement = output / "prepared-arrangement.npz"
    atomic_npz(arrangement, **arrays)
    chart_archive = output / "outer-chart-ancestry.npz"
    atomic_npz(chart_archive, **chart_replay)
    spec = dict(schema="superior-gluteal-finite-owner-pack-spec/v1", source=str(arrangement),
                frame=whole.frame, origin=origin.tolist(), surfaces=surface_map, corner_words=word_map,
                corner_uvs=uv_map, cells=cells, part_ids=PART_IDS, max_pairs=250000, max_tests=2000000)
    atomic_json(output / "pack-spec.json", spec)
    source_files = [plan_path, inventory, profiles_path, source / "finite-wall.npz",
                    source / "surfaces-and-ancestry.npz", source / "ports.json", source / "evidence.json",
                    source / "closed-blood-classifier.npz", *paths]
    if "inlet_correction_evidence" in plan:
        source_files.append(inlet_evidence)
    prepared = dict(schema="superior-gluteal-five-owned-cells-prepared/v1",
                    source_directory=str(source), source_geometry_sha256=whole.geometry_sha,
                    source_catalogue_sha256=plan["source_catalogue_sha256"], terms=TERMS, part_ids=PART_IDS,
                    interfaces=interfaces, sheet_regions=int(number), cells=cell_topology,
                    whole_wall_topology=whole_edges, whole_wall_cancellation_exact=True,
                    blood_classifier_unchanged=True, central_lumen_disks=0,
                    native_chart_counts=chart_counts,
                    native_chart_sources=[dict(component=i, path=str(path), sha256=digest_file(path),
                                               geometry_sha256=mesh.geometry_sha, attributes_sha256=mesh.attributes_sha,
                                               original_component_path=str(base / f"component-{i:02}.npz"),
                                               chart_basis=("literal-original-component" if path == (base / f"component-{i:02}.npz").resolve()
                                                            else "explicit-construction-component-chart; transport proof remains separately required"))
                                          for i, (path, mesh) in enumerate(zip(paths, sources))],
                    literal_corner_count_basis="The explicit component chart used to construct this candidate; corrected charts do not claim unchanged pre-correction normal words.",
                    explicit_source_component_paths="source_component_paths" in plan,
                    source_files={str(path): digest_file(path) for path in source_files},
                    prepared_files={str(path): digest_file(path) for path in
                                    (arrangement, chart_archive, output / "pack-spec.json")},
                    producer_sha256=digest_file(Path(__file__)),
                    preparation_passed=True, packing_and_decoded_checks_pending=True,
                    source_boundary_support="The parent retains each rounded fork and the proximal child support up to the first intact paired native ring.",
                    anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic() - began)
    atomic_json(output / "prepared.json", prepared)
    return prepared


def pack_prepared(source_directory, output_directory):
    began = time.monotonic()
    prepared = json.loads((output_directory / "prepared.json").read_text())
    if not prepared.get("preparation_passed") or Path(prepared["source_directory"]) != source_directory:
        raise ValueError("A matching saved arrangement is required")
    if prepared["producer_sha256"] != digest_file(Path(__file__)):
        raise ValueError("Prepared arrangement belongs to a different producer revision")
    for path, expected in (prepared["source_files"] | prepared["prepared_files"]).items():
        if digest_file(Path(path)) != expected:
            raise ValueError(("A bound input changed after preparation", path))
    spec = json.loads((output_directory / "pack-spec.json").read_text())
    data = np.load(spec["source"], allow_pickle=False)
    result = pack(data["positions"], {k: data[v] for k, v in spec["surfaces"].items()},
                  spec["cells"], spec["part_ids"], spec["frame"], output_directory / "finite-cells",
                  origin=spec["origin"], corner_words={k: data[v] for k, v in spec["corner_words"].items()},
                  corner_uvs={k: data[v] for k, v in spec["corner_uvs"].items()},
                  max_pairs=spec["max_pairs"], max_tests=spec["max_tests"],
                  checkpoint_directory=output_directory / "validation-checkpoints")
    if not result["passed"]:
        raise ValueError(("Decoded finite owner checks failed", result["checks"]["failures"]))
    for name, cell in result["checks"]["cells"].items():
        components = cell["topology"]["components"]
        if len(components) != 1 or components[0]["closed_orientable_genus"] != EXPECTED_GENUS[name]:
            raise ValueError(("Emitted owner changed topology", name, components))
    evidence = dict(schema="superior-gluteal-five-owned-cells-packed/v1", source_geometry_sha256=prepared["source_geometry_sha256"],
                    source_catalogue_sha256=prepared["source_catalogue_sha256"], finite_cells=5,
                    represented_unique_terms=list(TERMS.values()), part_ids=PART_IDS,
                    interfaces=prepared["interfaces"], whole_wall_cancellation_exact=True,
                    blood_classifier_unchanged=True, all_emitted_cells_checked=True,
                    native_chart_counts=prepared["native_chart_counts"], bank_sha256=result["bank_sha256"],
                    bank_bytes=result["total_bytes"], finite_owner_checks_passed=True,
                    upstream_recipient_context_and_render_acceptance_separate=True,
                    anatomical_acceptance=False, installed=False, elapsed_seconds=time.monotonic() - began)
    atomic_json(output_directory / "evidence.json", evidence)
    return evidence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-directory", required=True)
    parser.add_argument("--output-directory", required=True)
    parser.add_argument("--stage", choices=("prepare", "pack", "all"), default="all")
    parser.add_argument("--origin", type=float, nargs=3)
    args = parser.parse_args()
    source, output = local_path(args.source_directory), local_path(args.output_directory)
    if args.stage == "pack" and args.origin is not None:
        parser.error("Position origin is fixed by the saved prepare stage")
    if args.stage in ("prepare", "all"):
        result = prepare(source, output, args.origin)
    if args.stage in ("pack", "all"):
        result = pack_prepared(source, output)
    print(json.dumps(dict(stage=args.stage, output_directory=str(output),
                          preparation_passed=result.get("preparation_passed"),
                          finite_owner_checks_passed=result.get("finite_owner_checks_passed"),
                          elapsed_seconds=result["elapsed_seconds"], installed=False)))


if __name__ == "__main__":
    main()
