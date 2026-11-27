---
name: atlas-shape-search
description: Search indexed Human Atlas meshes by semantic metadata and pgvector shape
  similarity to find candidate source geometry for further anatomical review.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

`search --query kidney --model local-female-ta98 --side left` first filters actual source metadata. Filters include `--model`, `--sex`, `--side`, `--term`, `--system` and `--license`. Unknown metadata is not guessed. Hidden/suppressed donors require `--include-hidden`.

`search --similar INSTANCE_SHA --model TARGET_MODEL --limit 10` ranks the filtered target cohort by a deterministic 64-value surface-area/D2/normal/moment descriptor. Use the full source instance SHA for cross-model queries; a part name is resolved against `--model`. PostgreSQL provides a pgvector HNSW index, while the current filtered query ranks exactly to avoid losing small-cohort results after approximate filtering.

Descriptors are invariant to rigid motion, reflection and uniform scale up to floating/sampling arithmetic. This intentionally loses orientation/size information. Compare dimensions, handedness, side, source frame and anatomical landmarks before choosing a donor. Similarity is retrieval, never registration, TA98 identity, mesh quality, acceptable placement or permission to redistribute source models.
