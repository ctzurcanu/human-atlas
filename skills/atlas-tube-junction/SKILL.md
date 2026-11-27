---
name: atlas-tube-junction
description: Round a tube course, transport corresponding inner/outer walls, audit
  patency or join explicitly matched anatomical mouths without plugging the lumen.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

Use `tube path SPEC.json --output NEW.npz` for a rounded guide: `frame`, metre `waypoints`, `method: fillet`, `bend_radius_m`, optional `spacing_m`; `method: centripetal` accepts `chord_tolerance_m` instead.

`tube sweep SPEC.json --output NEW.npz` takes `source`, `path`, `ring_size`, optional `frame`, `keep_source_roll` and `tolerances` (the `TubeTolerances` fields in `geometry_tools/algorithms/tube_sections.py`). Source NPZ: `fields[K,N,3]`, `faces[F,3]`, `frame`, optional `geometryOrigin`. Controls must be ordered into corresponding complete station rings. All layers receive the same proper station transform; no radius reduction hides a tight-turn failure. Float32 export is rechecked.

`tube join SPEC.json --output NEW.npz` takes `first`, `second`, ordered `first_mouth`/`second_mouth` indices and explicit position/area/tangent tolerances. Matched mouths share controls; unrelated controls stay literal. `tube audit SOURCE.npz --frame FRAME` checks topology, every wall interval and actual interface/self intersections.

A guide-wall crossing check covers complete guide segments. Its sampled radius is not a continuous inscribed-radius proof. Whole physiological flow still needs adjoining organs, inlet/outlet and wall-versus-lumen review. Never leave an inlet tube continuing inside its recipient or close an opening with a solid plug.
