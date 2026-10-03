# Embryo 3month — GitHub Pages preparation

The model selector uses **Embryo 3month**. The stable URL parameter remains `model=embryo`. Production and development load `public/models/embryo-3month/atlas.json`; there is no dependency on a development-only `/local-models/` route.

The publication collection preserves 166 original CS23 meshes, their 1,072,606 triangles and original material colors. Original terminal Skin, Umbilical arteries and Umbilical vein geometry is restored in this collection. Eight HRA placental/cord objects remain separately selectable. Contact/overlap between the independent source surfaces does not establish a complete tissue assembly, continuous vessel lumens or stage-matched anatomy.

CS23 source is © 2016 Department of Medical Biology, Academic Medical Center, University of Amsterdam, Bernadette S. de Bakker and the atlas project team, CC BY-NC-ND 4.0. HRA context is Kristen Browne and Heidi Schlehlein / HuBMAP, CC BY 4.0. The source licenses remain separate. Credits and license links are in both `ATTRIBUTION.md` and `public/ATTRIBUTION.md`, plus the model-specific notice. The source stage remains CS23, 56–60 days after fertilization; the requested display name does not establish a three-month specimen.

## Build and verification

The source-preserving public collection is already staged in `public/models/embryo-3month/`, with content-addressed gzip chunks and `release-manifest.json`. All files needed by Pages are ordinary public assets; no ignored review assets must be forced into Git.

Run:

```sh
npm run check
npm run test:explosion-plane
npm run test:embryo-pages
npm run test:model-download
VITE_BASE_PATH=/human-atlas/ npm run build
```

The existing `.github/workflows/pages.yml` runs these checks, uses the repository Pages base path, uploads `dist`, and deploys on a push to `main` or manual workflow dispatch. No push or deployment was performed as part of preparation. The release includes `.nojekyll` and content-versioned geometry in the asset manifest. Gzip buffers decode through the native browser API or the bundled fallback, including hosts that already decode Content-Encoding.

The five-model explode-plane regression covers orbit, upside-down views and camera roll. Exploded clusters occupy one shared plane perpendicular to the POV. Selected clusters rotate about their own centers so a multiple selection cannot swing their centers off that plane. The rule is also recorded in `AGENTS.md`.

## Reproducible local source preparation

`node server/prepare-embryo-pages.mjs --stage-public` reproduces the collection from the existing ignored local study and context, checks the pre-cut study hash, and writes a release manifest. This requires the existing local source evidence; CI consumes the staged public files instead. The source study hash is `5bddc6c1392898cbf06e3225103d8b1bcd897e8350bead065966c4bb98bd9630`.

The optional `npm run build:embryo-pages-preview` stages a separate local preview package. Local preview public assets use hard links to avoid duplicating the large atlas resources. No ignore rules are changed.

## Coordinate reference and camera behavior

The embryo's model-local coordinate origin is the spinal axis at the umbilical attachment level. `coordinateFrame.origin` stores that point in the preserved native buffer coordinates; source geometry is not deformed to set the reference. The current point is inferred from the recorded cord attachment and a section through the source spinal cord. Its specimen landmark correspondence remains provisional. Evidence is in `reports/embryo-cs23/coordinate-origin.json`.

A drag on visible anatomy rotates about its picked surface point. An empty-space drag uses the umbilical-level spinal origin. Both permit complete turns without pole limits or a jump at drag start. Explode starts at the exact current view and fits the current moving extent as necessary. Full Implode restores that same original view, including pan, roll, zoom and projection offset, even when the original view was cropped. Those actions retain the viewing angle and the shared plane perpendicular to POV. Selection, component/relationship inspection, sections and depth edits never activate fitting. Tests: `npm run test:camera-intent`, `npm run test:explosion-plane`, `npm run test:explosion-viewport`, and `node --experimental-strip-types scripts/validate-camera.mjs`.
