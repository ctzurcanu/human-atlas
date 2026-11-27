---
name: atlas-mesh-spatial
description: Find nearby mesh candidates, build or reuse triangle BVHs and measure
  closest surface points in local Human Atlas coordinates.
---

Run commands from `/Users/christiantzurcanu/Documents/dev/human-atlas` with `.local-models/geometry-tools-env/bin/python -m geometry_tools`. Use `--report /absolute/new-file.json` **before** the subcommand when retaining evidence. Inputs use metres and explicit coordinate frames. PostgreSQL defaults to the dedicated local `human_atlas_geometry` database; `ATLAS_GEOMETRY_DSN` may select another local database. Source arrays stay file-backed.

- `neighbors PART --model MODEL --distance 0.002` uses PostgreSQL 3D bounding boxes to shortlist candidates. Include `--include-hidden` only when donor variants are relevant. A truncated shortlist is incomplete.
- `index PART --model MODEL` builds a checksum-verified persistent BVH for the actual indexed triangle sequence.
- `closest FIRST SECOND --model MODEL` gives triangle IDs, barycentric witnesses and Euclidean distance. Different frames require an explicit reviewed registration before comparison.
- `surface-query MESH QUERIES.npz --output NEW.npz --certify` uses Open3D proposals and verifies nearest points against original Float64 triangles/BVH. Queries require `points` [N,3] and matching scalar `frame`. Without `--certify`, distances are upper bounds only. `--rays` takes `origins`/`directions` and gives approximate first hits, with explicit miss masks. Reuse `SurfaceIndex` in Python for batch queries. Automatic certified batches fall back to Float64 when thin triangles collapse in the native normalized Float32 copy.
- `audit MESH --open3d` retains the existing exact geometry gates and adds independent C++ diagnostics. `open3d-status` checks installation. Native diagnostics/raycasts cannot establish wall thickness, patency or anatomical acceptance alone.

Bounds select candidates; they do not prove contact, occupied space or attachment. Use the contact tool for patches/crossings and the constraint registry for declared mechanical relationships. Never move anatomy from an AABB result alone. Cache identity includes actual geometry, tree parameters and implementation hashes.
