---
name: atlas-mesh-registry
description: Index, find or decode local Human Atlas mesh sources, metadata and authoritative
  TA98 descendant inventories. Use before replacing geometry or tracing a stale source.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

- `status` checks PostgreSQL and current revisions; `db-init` creates the workbench schema in the selected dedicated database.
- `registry-import CATALOGUE --model MODEL` hashes the JSON **and every referenced binary bank**. Include hidden/suppressed source variants. World coordinates apply `geometryOrigin` and `sourceOffset` exactly once.
- `export PART --model MODEL --output NEW.npz` retains indexed source controls, raw normal words/signs, UVs and colors. `register NEW.npz --part-id ID --name NAME` registers a candidate without changing an atlas.
- `inventory TA98_CODE --model MODEL --classifications FILE.json` inventories every descendant from the local authority SQLite cache. The classification map uses `whole-organ`, `tissue-region`, `wall-layer`, `opening`, `boundary-curve` or `microscopic-structure`. Unclassified terms and missing distinct geometry remain explicit. Save the inventory and satisfy the repository's checked-in anatomical report requirement before reconstruction.

An exact term mapping proves identity only. It does not prove a distinct layer or microscopic representation. Unknown source/license and source registration remain unknown; never infer redistribution rights. Reimport after source changes rather than bypassing a stale-bank error.
