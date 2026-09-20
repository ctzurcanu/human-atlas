# Human Atlas

An interactive 3D anatomy explorer built with React, Three.js, and shadcn/ui. The default **Male · detailed** reference combines **4,391 selectable pieces** from Z-Anatomy and Open 3D Model, adapted by Brian Pridgen. Switch to the original BodyParts3D male reference (**2,234 selectable meshes**, **3,432 named concepts**) and the Human Reference Atlas female reference (**1,038 selectable meshes**, **1,253 named concepts**).

**[Explore the live demo](https://human-atlas-seven.vercel.app)**

## Explore

- Choose male or female anatomy from the selector below the title.
- Orbit, zoom, and select structures directly on the body.
- Toggle individual systems or use skeleton and organ presets.
- Move from assembled anatomy to a spaced inventory of every visible piece.
- Search anatomical names and source identifiers.
- Isolate a selected structure and read its details.
- Use compact controls and detail panels on mobile.

## Run locally

Requires Node.js 22.13 or newer. No API keys or accounts are needed.

```sh
npm ci
npm run dev
```

Open http://localhost:3016. To build the static site, run `npm run build`; the output is in `dist/`.

## Validate

```sh
npm run check
node scripts/validate-atlas.mjs
node scripts/validate-atlas.mjs atlas-female.json
node scripts/validate-atlas.mjs atlas-male-detail.json
node scripts/validate-interactions.mjs
npm run build
```

Validation covers mesh buffers, names and concept membership, nonoverlapping exploded layouts at desktop and mobile aspect ratios, search and inspection contracts, and tap-versus-drag handling. Browser interaction checks have exercised selection, system controls, search, isolation, rotation, and 390×844, 320×568, and 844×390 layouts. Phone controls stay clear of the exploded inventory, and isolated structures fit the space above or beside the detail panel. Physical-device performance and real multitouch hardware have not been tested.

## Anatomy data

The default detailed male reference includes 564 muscle/tendon pieces, 537 nervous-system pieces, 717 connective-tissue pieces, and 680 muscle attachment markers. It imports all 769 Open 3D Model structures present in the source catalog, including the adapted upper limbs and brachial plexus. Attachment markers, surface regions, 64 fascia pieces, and 27 source-designated schematic details have separate layers, hidden by default. Schematic geometry retains its source placement and is not presented as anatomically positioned in the default body. Counts represent modeled pieces, not unique anatomical organs; coverage and granularity differ from BodyParts3D. The original male reference remains selectable.

The geometry is adapted from Brian Pridgen’s Anatomy Atlas; its viewer software is not copied. Geometry is hash-pinned, transformed to the stage, simplified within a 0.2% relative error limit per structure, and repacked. The detailed reference has about 4.3 million triangles and a 52 MB compressed download. Some source catalog entries have no identifiable geometry and are omitted from search instead of providing empty selections. Their names are recorded in the manifest. Mirrored anatomy retains its provenance.

**Detailed male geometry has mixed licenses**, including CC BY-SA 4.0, CC BY-NC 4.0 (kidneys), and CC BY-NC-SA 4.0 (inner-ear components). The noncommercial restrictions apply to those components and outputs containing them. The original MIT application license does not relicense the anatomy. Preserve [the complete source notice](public/MALE-ANATOMY-NOTICE.md) and [our adaptation credits](public/ATTRIBUTION.md).

The original male viewer uses **BodyParts3D 4.0**, an adult male reference anatomy, licensed **CC BY 4.0**. It does not represent every human structure or variation. Individual source meshes are distinct from named concepts, which may group multiple meshes. Descriptions distinguish general system context from individual organ explanations.

The female viewer uses the **Human Reference Atlas / HuBMAP 3D Reference Organ Set for Female v1.10**, licensed **CC BY 4.0**. It includes whole-body surface, selected organs, and female reproductive anatomy, plus 76 female-source lower-limb muscles from Andreassen et al. (2023). There are 90 muscle pieces in total after replacing two duplicate thigh muscles. Upper-body muscle and skeleton coverage remains incomplete; this is not a complete counterpart to the male atlas. No male-derived meshes are used. The lower-limb study surfaces are approximately fitted to HRA bones (17–36 mm source fit residuals), not a validated single-scan assembly. Eight pregnancy reference pieces are available in a separate layer, hidden by default. The two collections have different coverage.

Geometry is simplified for browser performance while retaining every source mesh. The packaged male model contains 2,288,268 triangles and downloads approximately 33 MB of compressed geometry. Full credits, source links, and adaptation details are in [ATTRIBUTION.md](public/ATTRIBUTION.md).

This is an educational explorer, not a diagnostic or surgical tool.

## How it works

Geometry is merged into batches. Per-structure GPU textures control translation, visibility, and selection, while component geometry supports accurate picking. Exploded layouts pack only the visible pieces. Rendering updates when the scene changes; orbit controls remain responsive without thousands of separate draw calls.

The optional WebMCP tools expose anatomy search and inspection in compatible browsers. The visible interface works without them.

## Rebuilding geometry

The repository includes browser-ready geometry. Rebuilding it is optional: obtain the official BodyParts3D OBJ archive and English metadata tables, prepare the joined concepts and display-system mappings, run `scripts/convert-anatomy.py`, then `node scripts/optimize-anatomy.mjs` and `node scripts/compress-models.mjs`. Simplification uses a 0.2% relative error limit per structure.

## Deploy

Import this repository into Vercel as a Vite project. The included `vercel.json` configures `npm ci`, `npm run build`, and the `dist` output directory. It can also be served by a static host.

## License

Original application code is released under the [MIT License](LICENSE). **The anatomy data has its own CC BY 4.0 license**; preserve the attribution when redistributing it. Third-party dependencies retain their respective licenses.

Issues and pull requests are welcome. Please include reproduction steps and browser/device details for interaction problems.

To reproduce the expanded female assets, run `node scripts/upgrade-female.mjs`. It downloads a pinned CC BY geometry adaptation, excludes all male-derived meshes, removes duplicate rectus femoris surfaces, and repacks female-only buffers. See the attribution file for original datasets and the geometry adaptation credit.

Rebuild the detailed male reference with `node scripts/import-male-detail.mjs [CACHE_DIRECTORY]`. The importer uses the catalog and SHA-256 manifest in `scripts/data/`, fetches only public geometry, rejects changed assets, resolves unique naming differences, and excludes duplicate or unidentified geometry.

### Display and model quality

The interface and 3D stage follow the operating system’s light/dark preference, including changes while the viewer is open. Wheel zoom follows the pointer; selecting a structure makes it the orbit pivot, including when other anatomy remains visible.

The model menu includes **Male · full resolution + skin**, preserving all 11,647,803 triangles in the imported named source meshes (124 MB compressed). The existing detailed option retains 4,428,317 triangles for lighter devices. Both have the same named anatomy and component licenses. The full-resolution option initially shows the derived outer body surface. The Systems panel’s skin opacity slider can fade or remove it; this surface is not a histological skin-layer model. Rendering still uses the viewer’s system colors, not the source’s textures.

Regenerate full resolution using `node scripts/import-male-detail.mjs /tmp/male-atlas-source --full`. Validate with `node scripts/validate-atlas.mjs atlas-male-full.json`.

### Dissection, selection sets, and view URLs

Open **Study** to choose body regions and camera directions, keep multiple structures solid against adjustable transparent context, peel system layers, or move axial/sagittal/coronal clipping planes. **Add to selection set** works in search; Shift-click adds meshes directly. **Focus selection** zooms in while retaining context. **Isolate selection set** removes all context. Undo/redo covers dissection and selection changes; continuous slider changes are grouped. Selected pieces survive broad layer peeling, but remain subject to cross-sections.

Labels attach to actual triangle centroids, moving with the model. Dashed leaders and “behind” identify anchors obscured by selected geometry. Up to 20 labels are drawn simultaneously; every selected item remains listed in Study. Cuts expose open mesh boundaries without synthesizing tissue interiors. Regional filtering uses source memberships where provided and spatial bounds for older datasets. The whole-body surface is hidden in regional views.

The platform and floor have been removed. Zoom is cursor-centered, supports close inspection, and limits zoom-out relative to the current fitted view.

Use **Share view URL** to get a URL containing the model, selected IDs, camera position/target/projection offset, visible layers, region, opacity, hidden pieces, peeling and section. The URL can also be authored directly, for example:

`http://localhost:3016/?model=male-detail&select=Stomach&select=Pancreas&view=back&context=0.15`

`select` accepts a concept ID, source ID, or exact display name and can be repeated. `view` accepts `three-quarter`, `front`, `back`, `side` (patient-left), `right`, `superior`, and `inferior`. `focus=1` frames a selection when a camera is not supplied. `cut=axial|sagittal|coronal` and `slice=0..1` control sections; `flip=1` reverses the retained side. Unknown names are ignored. URLs are self-contained and do not upload views to a server.

Source coverage is audited by `node scripts/validate-source-coverage.mjs` against the downloaded GLBs and reference metadata. The import now applies 70 source node aliases, including laterality corrections, and whole-structure representation groups. All 4,391 geometry-backed source records are available, including 254 tendon/cartilage patches stored in parent mesh materials. Canonical aliases and whole-structure handles suppress duplicate geometry, matching the reference’s 4,350 physical representations. Catalogue names without direct geometry are listed in `public/SOURCE-COVERAGE.json`; no geometry is fabricated for them. Only factual anatomy metadata and licensed geometry are imported; no reference viewer code is included.

Additional checks: `node --experimental-strip-types scripts/validate-study.mjs` and `node --experimental-strip-types scripts/validate-camera.mjs`.
