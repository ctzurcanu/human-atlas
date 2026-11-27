---
name: atlas-mesh-contact
description: Measure actual near-contact patches, crossings and gaps between local
  anatomy meshes; distinguish proximity from a declared anatomical interface.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

`contact FIRST SECOND --model MODEL --distance 0.002 --max-pairs 250000` tests actual triangles through cached BVHs. The NPZ artifact retains both face IDs, closest world points, barycentric coordinates, normal relationship and crossing/coplanar/touching classification. `--self` requires identical triangle sequences.

Use `--area-error 0.000001 --area-max-depth 8` when patch area matters. Report its lower/upper bounds and unresolved uncertainty; the sum of touched source-face areas is only an **upper bound**, not the contact area. Depth/work/pair truncation must remain visible. Float64 predicates and numerical margins are not exact arithmetic or anatomical tolerances.

A gap is never welded or called a join automatically. Crossing counts alone do not establish tissue penetration: distinguish wall, lumen and occupied tissue, then review actual surfaces and sections. Retain the cache artifact rather than rerunning the same numerical solve.
