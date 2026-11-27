---
name: atlas-geometry-fit
description: Create rigid, uniform-scale, constrained-translation or explicit cage
  candidates for local anatomy while preserving dependencies, normals and fixed skeleton
  geometry.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

`fit SPEC.json --output NEW_DIRECTORY` writes candidates and checks; it never installs them. Required spec fields: `model`, `moving` (mesh references), and `mode`.

Prefer `rigid` or `similarity` with corresponding `source_landmarks` and `target_landmarks`, optional positive `weights` and `scale_bounds`. At least three noncollinear points determine a proper rotation. Reflection and silent anisotropic scale are rejected. `translation` takes `translation_m`, optional `translation_bounds_m`, and reviewed `halfspaces` (`normal`, `minimum_dot_m`). Linearized clearance is only a proposal; check actual resulting triangles.

Use `cage` only for an explicit longitudinal cage spec accepted by `geometry_tools.algorithms.longitudinal_organ_cage.LongitudinalOrganCage`. One field acts on the whole closure; normals use the inverse transpose. Its continuous positive Jacobian does not prove discrete mesh quality. Tube walls should use rigid section transport instead of cage shear.

`icp` uses Open3D to refine an explicit `initial_matrix` or landmark initialization against `registration_target` (registered reference or mesh NPZ). Declare `registration_source` inside the moving closure and `target_frame`. The `icp` object requires `max_correspondence_m`, decreasing `voxel_sizes_m`, `max_correction_m` and `max_rotation_degrees`; optional `method` is `point-to-point` or `point-to-plane`. Sampling is area weighted and reproducible. ICP preserves initialized uniform scale, transforms the entire closure together and still runs all geometry/constraint gates. It requires an anatomically justified target and does not prove anatomical correspondence or install candidates. Jobs bind target geometry and native library versions.

`protected` adds explicit blockers. The complete closure, UVs, packed normal signs, source faces, fixed context and all touched constraints must survive. Failed/pending checks produce a diagnostic candidate, never accepted anatomy. A cross-model transfer needs explicit `target_frame`, landmarks/scale and target context constraints; similarity retrieval is not registration. Review emitted coordinates and sections before guarded installation.
