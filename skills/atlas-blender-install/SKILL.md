---
name: atlas-blender-install
description: Queue surface/section review in the existing Blender window, track reproducible
  geometry jobs and install a reviewed local TA98 revision with rollback.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

`review prepare SPEC.json` takes `model` or explicit `frame`, `meshes` (references or `{mesh,name,color}` entries), optional `views` and bounded `resolution`. A view has `name`, `direction`, and optionally diagnostic `section: {normal,point,clear_outer}`. It returns a request ID and bootstrap Python line. Run that line through native computer control in the **existing** Blender console; never launch another instance, factory-reset, open a replacement main file or save the user's main file. The data-only queue renders separate scenes and verifies original objects/transforms. `review status ID` verifies image hashes. Inspect the actual images; diagnostic sections are explicitly uncapped.

For an assembled finite-cell preview, `assembled_tissues: true` plus each mesh's `part_id` and actual `tissue_draw_groups` hides only internal charts whose explicit neighbour is also displayed. It does not infer a tissue attachment. Diagnostic clipping preserves imported normal fields.

`job run SPEC.json` accepts `kind` plus `parameters`, optional `meshes`/`model`. Kinds: audit, contact, fit, tissue-layers, tissue-pack, tube-path, tube-sweep, tube-join, catalogue-candidate. A catalogue-candidate job runs the install preparation spec as its parameters and binds the actual audited output banks. Inputs, source meshes, code and results are hash-bound. `job status ID` reports durable progress. Reuse passed results; duplicate running jobs are rejected. `--retry` is explicit, bounded to three attempts. Recover a crashed owner only with `job recover ID --reason TEXT`; a live owner blocks recovery.

`job review ID EVIDENCE.json` requires the current `input_signature`, `reviewer`, `accepted_for_local_preview: true`, `remaining_limits`, and actual inspected surface/section `views` with paths/hashes. This is a recorded review decision, not automatic anatomical acceptance.

`install prepare SPEC.json --output NEW_DIRECTORY` checks a complete `proposed_catalogue`, `model`, `expected_catalogue_sha256`, exact `changed_parts`, source-bound classified `inventory`, `assembly.whole_organ`, `assembly.adjoining_interfaces`, and `job_id`. It protects unchanged items/skeleton and independently checks emitted replacements. `install apply PROPOSAL/install.json --name UNIQUE_LOCAL_REVISION` requires reviewed job artifacts binding every emitted bank. It updates local catalogues, visit source bindings, TA98 hierarchy and ownership forests through the existing revision adapter. `install rollback ID` rejects later edits and restores saved files. Never change .gitignore, stage, push or deploy as part of this skill. Reimport the resulting catalogues into the registry after installation or rollback.
