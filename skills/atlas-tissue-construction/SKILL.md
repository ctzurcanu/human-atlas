---
name: atlas-tissue-construction
description: Construct finite anatomical wall layers or explicitly owned tissue cells
  with a single shared geometry bank and reciprocal interfaces.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

- `tissue layers SOURCE.npz --names OUTER_LAYER INNER_LAYER --frame FRAME --output NEW_DIRECTORY` consumes explicit `fields[K,N,3]` and carrier `faces[F,3]`. K interfaces make K-1 finite layers; the names must match that count. Fields run outside to inside.
- `tissue pack SPEC.json --output NEW_DIRECTORY` consumes a source NPZ with `positions` and named surface-face arrays. The spec supplies `source`, `frame`, `surfaces` (surface name to array key), `cells` (surface/side/neighbor references) and unique `part_ids`; optional `corner_words` maps surfaces to explicit Int16 normal arrays.

Actual open carrier boundaries receive annular thickness rims. No disk spans a lumen. Shared interfaces are stored once, reused with opposite orientation and identical normal words/sign inversion. The output includes a renderable packed bank, catalogue records, source arrangement, decoded cell NPZs, byte counts and checks on the actual emitted Float32 geometry.

Reject reversed/collapsed prisms, bad topology/self-crossings, nonpositive tissue volume or occupied-volume overlap. Engine volume checks supplement independent triangles/topology. The tool supplies no anatomical thickness, regional boundaries, materials, microscopic anatomy or implicit solid openings. Establish and review those inputs under AGENTS.md first. Do not call a generic concentric shell a completed TA98 layer.
