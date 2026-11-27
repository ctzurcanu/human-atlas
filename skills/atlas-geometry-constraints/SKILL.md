---
name: atlas-geometry-constraints
description: Manage explicit anatomical dependencies, fixed structures and interface
  constraints before fitting or moving a Human Atlas assembly.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

- `constraints import-forest FOREST.json --model MODEL` binds the existing complete ownership forest to the actual catalogue revision. It rejects stale or incomplete forests and fixes skeletal items by policy.
- `constraints closure PART ... --model MODEL` returns the entire move/scale closure, including retained donor members and fixed blockers. Ownership and explicit `co_move` declarations propagate; proximity never does.
- `constraints add SPEC.json --model MODEL` stores `kind`, `first`, optional `second`, `evidence.source` and `parameters`. Kinds: fixed, shared-interface, disjoint, nested, adjacent-unjoined, sliding, open-conduit, clearance.
- `constraints relations --model MODEL --kind anatomical:innervates --mesh PART` retrieves preserved forest relationships without inferring geometry. `constraints refine ID SPEC.json` adds reviewed geometric parameters to a pending contract and retains its previous declaration.
- `constraints check --model MODEL` checks declarations with geometric parameters; imported declarations lacking reviewed face/ring/clearance parameters stay pending and cannot pass.

For shared interfaces provide `first_faces`/`second_faces`; triangles must be literally equal with opposite winding. For conduits provide ordered `first_ring`/`second_ring`. For clearance/sliding provide `minimum_m` and optional `maximum_m`. `disjoint` with `occupancy: solid` requires closed volumes and checks containment as well as surface crossings. Nested relations require a closed container. A matched collar does not prove an unobstructed whole lumen.
