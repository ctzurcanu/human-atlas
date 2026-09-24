# Human Atlas

An interactive 3D anatomy explorer built with React, Three.js, and shadcn/ui. The default **Male · detailed** reference has **4,391 selectable pieces**. Switch to the standard male reference (**2,234 selectable meshes**, **3,432 named concepts**) or the female reference (**1,038 selectable meshes**, **1,253 named concepts**).

**[Explore the live demo](https://ctzurcanu.github.io/human-atlas/)**

## Explore

- Choose male or female anatomy from the selector below the title.
- Orbit, zoom, and select structures directly on the body.
- Toggle individual systems or use the All, Skin, Muscles, Vessels, Skeleton, and Organs presets.
- Move from assembled anatomy to a spaced inventory of every visible piece.
- Search anatomical names and source identifiers.
- Isolate a selected structure and read its details.
- Download a clean PNG of the current 3D view, without the interface.
- Save named views in this browser and reopen them from the share panel.
- Use compact controls and detail panels on mobile.

## Run locally

Requires Node.js 22.13 or newer. No API keys or accounts are needed.

```sh
npm ci
npm run dev
```

Open http://localhost:3016. To build the static site, run `npm run build`; the output is in `dist/`.

### Use the local Zygote male and female models

The `fetch/content` mirror can be imported into a development-only model directory:

```sh
npm run import:local-models
npm run dev
```

The importer reads `/Users/christiantzurcanu/Documents/dev/fetch/content` by default. Pass a different mirror directory after `--` if needed, for example `npm run import:local-models -- /path/to/content`. In the local viewer, choose **Male · local Zygote** or **Female · local Zygote**, or open `http://localhost:3016/?model=local-male&select=Stomach&focus=1&context=0.18&skin=0` and the corresponding `local-female` URL. The MCP server also accepts `model="local-male"` and `model="local-female"` while the local development server is running.

Converted geometry and source texture maps are written to `.local-models/`, which is Git-ignored and served only by the Vite development server. They are not copied into `dist` or published to GitHub Pages. Local models open with skin hidden so the textured anatomy is visible; use the Skin / body surface slider to show the textured skin.

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

The default detailed male reference includes 564 muscle/tendon pieces, 537 nervous-system pieces, 717 connective-tissue pieces, and 680 muscle attachment markers. Attachment markers, surface regions, 64 fascia pieces, and 27 schematic details have separate layers, hidden by default. Schematic geometry is not presented as anatomically positioned in the default body. Counts represent modeled pieces, not unique anatomical organs; coverage and granularity differ between models. The standard male reference remains selectable.

Geometry is hash-pinned, transformed to the stage, simplified within a 0.2% relative error limit per structure, and repacked. The detailed reference has about 4.3 million triangles and a 52 MB compressed download. Some catalogue entries have no identifiable geometry and are omitted from search instead of providing empty selections. Their names are recorded in the manifest.

The standard male option represents an adult male reference anatomy. It does not represent every human structure or variation. Individual meshes are distinct from named concepts, which may group multiple meshes. Descriptions distinguish general system context from individual organ explanations.

The female option includes whole-body surface, selected organs, and female reproductive anatomy, plus 76 female lower-limb muscle surfaces. There are 90 muscle pieces in total after replacing two duplicate thigh muscles. Upper-body muscle and skeleton coverage remains incomplete; this is not a complete counterpart to the male atlas. No male-derived meshes are used. The lower-limb surfaces are approximately fitted to the bones (17–36 mm fit residuals), not a validated single-scan assembly. Eight pregnancy reference pieces are available in a separate layer, hidden by default. The two collections have different coverage.

Geometry is simplified for browser performance while retaining every source mesh. The packaged standard male model contains 2,288,268 triangles and downloads approximately 33 MB of compressed geometry.

This is an educational explorer, not a diagnostic or surgical tool.

## How it works

Geometry is merged into batches. Per-structure GPU textures control translation, visibility, and selection, while component geometry supports accurate picking. Exploded layouts pack only the visible pieces. Rendering updates when the scene changes; orbit controls remain responsive without thousands of separate draw calls.

The optional WebMCP tools expose anatomy search and inspection in compatible browsers. The visible interface works without them.

## Rebuilding geometry

The repository includes browser-ready geometry. Rebuilding it is optional: obtain the original male OBJ archive and English metadata tables, prepare the joined concepts and display-system mappings, run `scripts/convert-anatomy.py`, then `node scripts/optimize-anatomy.mjs` and `node scripts/compress-models.mjs`. Simplification uses a 0.2% relative error limit per structure.

## Deploy

GitHub Pages publishes this project at https://ctzurcanu.github.io/human-atlas/ using `.github/workflows/pages.yml`. Pushes to `main` run type checking, build the static site, and deploy `dist`. The workflow can also be run manually. Repository Settings → Pages → Source must be set to **GitHub Actions**.

The workflow uses the Pages deployment base path for scripts, model catalogues, binary chunks, and icons. Shared-view URLs preserve this path. Local development continues at `/`.

To reproduce the Pages build locally:

```sh
VITE_BASE_PATH=/human-atlas/ npm run build
npx vite preview --base /human-atlas/ --port 4173
```

The `dist` directory can also be served by another static host; set `VITE_BASE_PATH` to its deployment path when building.

## Embed the viewer

Open a view, choose the model and layers you want, then use **Share and save this view**. Select which controls should appear in the iframe and choose **Copy embed code**. The generated iframe opens the saved view with only those controls. For example, this one includes Systems and a link to the full viewer, but excludes Explode:

```html
<iframe src="https://ctzurcanu.github.io/human-atlas/?model=female&amp;embed=1&amp;ui=systems%2Copen" title="Human Atlas interactive anatomy viewer" loading="lazy" style="width:100%;height:600px;border:0" allowfullscreen></iframe>
```

You can use `model=male-detail`, `male-full`, or `male` instead of `female`. Set a fixed or responsive height on the iframe; 600px gives the controls room, and the minimum supported height is 340px. The `ui` parameter is a comma-separated list of visible controls: `model`, `search`, `study`, `systems`, `camera`, `explode`, `details`, `open`, and `download`. Use `ui=` for a bare viewer. If `ui` is omitted, the default controls are model, search, systems, details, and open; Study, Camera controls, Explode, and PNG download are off by default. The Systems button works independently of Explode. Named views are stored only in the browser where they were saved; the generated URL and iframe code can be shared separately.

## Use with an MCP client

Human Atlas includes a local stdio MCP server. It searches the packaged anatomy catalogues and creates focused views of the deployed viewer. Run `npm ci`, then add this server to your MCP client's configuration (replace the path with your checkout's absolute path):

```json
{
  "mcpServers": {
    "human-atlas": {
      "command": "node",
      "args": ["/absolute/path/to/human-atlas/mcp/server.mjs"]
    }
  }
}
```

The server provides `search_anatomy` for names and IDs and `show_anatomy` for an interactive view. For example, call `show_anatomy` with `{"structure":"Stomach","model":"male-detail","controls":["systems","open"]}` to show Systems without Explode. Omit `controls` for the default interface, or pass `[]` for a bare viewer. To select a precise side or variant, use the ID returned by `search_anatomy`, such as `DETAIL:Sternocostal head of pectoralis major muscle.l`. The view tool returns a deployed URL and copyable iframe code. MCP Apps-capable clients can display the interactive viewer inline; other clients can open the URL or use the iframe HTML. The MCP process runs locally; GitHub Pages hosts the viewer only.

Run `npm run test:mcp` to verify tool calls, URL selection, and the UI resource.

## License

The application code is distributed under the [GNU General Public License, version 3](LICENSE) (GPL-3.0-only).

Issues and pull requests are welcome. Please include reproduction steps and browser/device details for interaction problems.

To reproduce the expanded female assets, run `node scripts/upgrade-female.mjs`. It downloads a pinned geometry adaptation, excludes all male-derived meshes, removes duplicate rectus femoris surfaces, and repacks female-only buffers.

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
