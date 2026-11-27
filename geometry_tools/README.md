# Local geometry workbench

Run from `/Users/christiantzurcanu/Documents/dev/human-atlas`:

```sh
.local-models/geometry-tools-env/bin/python -m geometry_tools --help
```

The dedicated PostgreSQL 17 database is `human_atlas_geometry`, using pgvector
and the 3D `cube` extension. `ATLAS_GEOMETRY_DSN` selects another **local**
database. No original mesh arrays are copied into PostgreSQL. Generated BVHs,
witnesses, candidates, job evidence and Blender requests live under
`.local-models/geometry-workbench`; `ATLAS_GEOMETRY_CACHE` overrides this for
ordinary numerical tools. The Blender queue uses the project's default path.
All coordinates and distances use metres and explicit frames.

Reusable mesh algorithms live in [`algorithms/`](algorithms/): tube sections
and junctions, finite walls, shared boundaries, Boolean provenance, constrained
fitting and storage repairs. Their regression fixtures live in
`tests/algorithms/` and run with the workbench tests below. New code imports
`geometry_tools.algorithms.MODULE`. The former `geometry/` contains compatibility
imports for frozen repair scripts; each forwards to the same canonical module,
without maintaining another implementation. Jobs hash nested algorithm sources.
Clarabel is the existing corridor-fitting solver and is included in the isolated
runtime requirements.

| Tool | Commands | Skill |
| --- | --- | --- |
| Mesh registry and TA98 inventory | `registry-import`, `register`, `export`, `inventory` | [Registry](../skills/atlas-mesh-registry/SKILL.md) |
| Persistent spatial index | `index`, `neighbors`, `closest` | [Spatial queries](../skills/atlas-mesh-spatial/SKILL.md) |
| Actual near-contact patches | `contact` | [Contact analysis](../skills/atlas-mesh-contact/SKILL.md) |
| Complete ownership and constraints | `constraints` | [Dependencies](../skills/atlas-geometry-constraints/SKILL.md) |
| Controlled assembly fitting | `fit` | [Fitting](../skills/atlas-geometry-fit/SKILL.md) |
| Shared finite tissues | `tissue` | [Tissue construction](../skills/atlas-tissue-construction/SKILL.md) |
| Tube course, wall and mouth repair | `tube` | [Tube/junction tools](../skills/atlas-tube-junction/SKILL.md) |
| Blender review, jobs, guarded local install | `review`, `job`, `install` | [Review and installation](../skills/atlas-blender-install/SKILL.md) |
| Filtered geometric shape retrieval | `search` | [Shape search](../skills/atlas-shape-search/SKILL.md) |
| Open3D native surface queries and diagnostics | `surface-query`, `open3d-status`, `audit --open3d` | [Spatial queries](../skills/atlas-mesh-spatial/SKILL.md) |

Each skill is also linked into the personal Codex skills directory. The project
files are the maintained sources. Command-specific `--help` and these skills
describe the actual JSON/NPZ inputs. Use `--report NEW.json` before a subcommand
to retain complete evidence; large report output is summarized on stdout.

Example on the installed female data:

```sh
.local-models/geometry-tools-env/bin/python -m geometry_tools contact \
  RENALDETAIL:female:left:capsule TA98PART:female:Suprarenal:left:capsule \
  --model local-female-ta98 --distance 0.002
```

The query returns closest points on triangles, source face IDs/barycentrics,
patch face sets, crossings and a checksum-verified reusable artifact. Touched
face area is only an upper bound. Optional adaptive area bounds retain their
uncertainty and work limits. Proximity never creates an anatomical attachment.

The fitting tools prefer proper rigid/similarity transforms, expand explicit
dependency closures and reject fixed skeleton members. Tube transport acts on
all corresponding walls together. Finite tissues share one encoded bank and
reciprocal normal signs; emitted Float32 topology, intersections, occupied
volumes and local wall intervals are checked independently.

Jobs bind actual mesh/attribute identities, input files, constraints and code.
`catalogue-candidate` jobs prepare and audit emitted catalogue banks before local
installation. Existing open reference sheets require a source-bound
`open_reference_refits` policy that replays their original parent triangles,
fixed controls, UVs and bounded local displacement. New finite tissues retain
the closed-wall gate.

Passed work is reused; concurrent duplicate jobs and unbounded retries are
rejected. Blender review runs a data-only queue in the existing window, creates
separate diagnostic scenes, retains source hashes, and verifies original
objects/transforms. It never launches Blender or saves its main file.

For assembled finite tissues, review specs may set `assembled_tissues: true`
and give each mesh its `part_id` and actual `tissue_draw_groups`. Only interfaces
to explicitly visible neighbours are hidden. Section clipping retains source
normals; diagnostic cuts remain uncapped. Tissue pack specs may provide
`corner_words` and `corner_uvs` array maps. UV seams remain distinct render
charts even when they share the same physical coordinates.

The local installer requires a complete classified TA98 inventory, declared
assembly/interfaces, independent emitted-geometry checks and actual reviewed
bank artifacts. It preserves old banks, protects unchanged items/skeleton and
updates visit bindings, hierarchy and forest through the existing reversible
revision adapter. Snapshot rollback rejects subsequent edits. No public atlas
was changed while developing these tools; the installer has transaction and
failure regression coverage. Local urinary revisions use this guarded installation path.

Numerical passes are not anatomical acceptance. Source terminology mappings,
generic shells, descriptor matches, open sheets and sampled lumen radii have
the explicit limitations recorded by their tools. Follow the repository's
anatomical manual and inspect actual surfaces/sections before accepting anatomy.

Open3D 0.20.0 runs on the CPU in the isolated Python environment. Its import is
lazy, so registry/database commands do not load the native library or GUI.
`surface-query MESH QUERIES.npz --output NEW.npz` accepts registered mesh IDs or
local mesh NPZs. Queries contain `points` [N,3] and a scalar `frame`; `--rays`
instead takes `origins` and `directions`. Frames must match. Reuse `SurfaceIndex`
in Python for repeated batches. Output retains original triangle IDs,
barycentric witnesses and source hashes. Ray misses have a false hit mask and
NaNs in NPZ arrays. No anatomy is modified.

Native raycasting uses an origin-shifted, normalized Float32 copy. Nearest-point
witnesses are recalculated on original Float64 triangles; fast distances remain
upper bounds until `--certify` performs the existing BVH search. Certified batch
queries are also available through `BVH.closest_points(backend="open3d")` and
`winding_contains(distance_backend="open3d")`. The existing BVH remains the
default for final checks: native seeding plus full certification can be slower
than the original search. A source triangle that collapses in Float32
uses the Float64 fallback in automatic mode. `audit --open3d` adds independent
C++ manifold/intersection diagnostics alongside the original exact predicates.
Native intersection candidates are replayed through the source-triangle
predicates so near-coplanar false positives remain explicitly distinguished.
Native diagnostics and approximate rays never replace geometry acceptance.

`fit` also supports `mode: "icp"`. Specify `registration_source` inside the moving
closure, `registration_target` (registered reference or mesh NPZ), `target_frame`,
and `initial_matrix` or corresponding landmarks. `initial_similarity: true`
permits landmark-initialized uniform scale; ICP itself refines only rigid pose.
Its required `icp` controls are `max_correspondence_m`, decreasing positive
`voxel_sizes_m`, `max_correction_m` and `max_rotation_degrees`; optional controls
include `sample_count`, `seed`, `iterations`, `minimum_fitness` and `method`
(`point-to-point` or `point-to-plane`). Use only an anatomically justified target
and initialization. The result applies one common transform to the complete
dependency closure and retains attributes, fixed-structure blockers and all
existing geometric/constraint checks. A failed registration raises an error;
no atlas installation occurs. Jobs bind target sources and library versions.

Algorithm references: [raycasting queries](https://www.open3d.org/docs/release/python_api/open3d.t.geometry.RaycastingScene.html)
and [ICP registration](https://www.open3d.org/docs/release/tutorial/pipelines/icp_registration.html).

Run all regression tests, including PostgreSQL transactions that roll back:

```sh
ATLAS_GEOMETRY_TEST_DB=1 .local-models/geometry-tools-env/bin/python \
  -m unittest discover -s geometry_tools/tests -v
```

The isolated Python runtime is reproducible from `requirements.txt` using
Python 3.12. The original atlas assets and `.gitignore` are unchanged.
