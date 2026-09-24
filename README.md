# Human Atlas

An interactive 3D anatomy explorer built with React, Three.js, and shadcn/ui. The default **Male · detailed** reference is built directly from Z-Anatomy and Open 3D Model publisher files: **4,405 selectable surfaces**, **3,801 named concepts**, and **5,764,617 triangles**. The standard BodyParts3D male and female references are also available.

**[Explore the live demo](https://ctzurcanu.github.io/human-atlas/)**

## Explore

- Choose male or female anatomy from the selector below the title.
- Orbit, zoom, and select structures directly on the body.
- Toggle individual systems or use the All, Skin, Muscles, Vessels, Skeleton, and Organs presets.
- Move from assembled anatomy to a spaced inventory of every visible piece.
- Search anatomical names and source identifiers.
- Isolate a selected structure and read its details.
- Download a clean PNG of the current 3D view, without the interface.
- After positioning the model, select the bookmark icon at the top right, name the view, and select **Save view** in the panel. Reopen saved views from that same panel. They are stored in this browser's localStorage.
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
node scripts/validate-primary-female.mjs public/models/atlas-hra-female.json
node scripts/validate-primary-atlas.mjs public/models/atlas-z-anatomy.json
node scripts/validate-interactions.mjs
npm run build
```

Validation covers mesh buffers, names and concept membership, nonoverlapping exploded layouts at desktop and mobile aspect ratios, search and inspection contracts, and tap-versus-drag handling. Browser interaction checks have exercised selection, system controls, search, isolation, rotation, and 390×844, 320×568, and 844×390 layouts. Phone controls stay clear of the exploded inventory, and isolated structures fit the space above or beside the detail panel. Physical-device performance and real multitouch hardware have not been tested.

## Anatomy data

The detailed male reference retains each original Z-Anatomy object and its material-defined surfaces, including tendon and cartilage patches. Original Open 3D Model upper-limb files contribute additional nerves, vessels, muscles, hand bones, ligaments, and sheaths; the right-side additions are mirrored to the left and labeled as such in the manifest. Blender reference lettering is excluded. Muscle attachments and fascia have separate layers, hidden by default. Counts represent selectable surfaces, not unique organs; coverage and granularity differ between models. The standard male reference remains selectable.

The original Z-Anatomy and Open 3D Model archives and source files are SHA-256 verified. The direct importer bakes source transforms, converts Z-up to the viewer's Y-up coordinates, splits material surfaces, and packs raw and gzip buffers without simplifying the source triangles. Original material colors are recorded as `sourceColor`; a subdued anatomy palette is used for display. The Z-Anatomy Blender archive contains no anatomical image textures. The Open 3D Model upper-limb GLB contains embedded images, but the current batched viewer uses its material colors without rendering those images.

The standard male option represents an adult male reference anatomy. It does not represent every human structure or variation. Individual meshes are distinct from named concepts, which may group multiple meshes. Descriptions distinguish general system context from individual organ explanations.

The female option has 1,030 selectable meshes and 1,234 concepts. It combines the original Human Reference Atlas united-female v1.10 GLB, eight pelvic bone surfaces from official HRA v1.5, and 74 original University of Denver Visible Human Female lower-limb muscle STLs. The eight placenta-group meshes are in the separate Embryo model. The other two muscle STLs duplicate HRA's rectus femoris meshes and are omitted. Upper-body muscles and skeleton are still incomplete; this is not a complete counterpart to the male atlas. The separately sourced lower-limb muscles are fitted to the HRA bones, with measured alignment errors in the manifest.

Geometry is simplified for browser performance while retaining every source mesh. The packaged standard male model contains 2,288,268 triangles and downloads approximately 33 MB of compressed geometry.

This is an educational explorer, not a diagnostic or surgical tool.

## How it works

Geometry is merged into batches. Per-structure GPU textures control translation, visibility, and selection, while component geometry supports accurate picking. Exploded layouts pack only the visible pieces. Rendering updates when the scene changes; orbit controls remain responsive without thousands of separate draw calls.

The optional WebMCP tools expose anatomy search and inspection in compatible browsers. The visible interface works without them.

## Rebuilding geometry

The repository includes browser-ready geometry. Rebuild the detailed male assets directly from the publishers' archives with Blender 4.3 or newer:

```sh
npm run import:primary-male
```

Set `BLENDER` to the Blender executable if it is not on `PATH`. The script fetches and hash-verifies the [original Z-Anatomy archive](https://github.com/Z-Anatomy/Models-of-human-anatomy) and [Open 3D Model upper-limb GLB](https://anatomytool.org/open3dmodel-create), opens both with source scripts disabled, runs the two direct importers, and validates every atlas buffer. It does not read geometry, catalogues, aliases, or labels from another anatomy viewer. To stage the output elsewhere, run `node scripts/import-primary-male.mjs --out /path/to/output`.

Rebuild the standard male reference with `npm run import:primary-standard`. This fetches and verifies the official [BodyParts3D OBJ archive and four metadata tables](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html), joins the IS-A and PART-OF concepts, applies Human Atlas's separately maintained element-to-system mapping, converts and optimizes the meshes, and validates the output before installing it. The tested rebuild matches every shipped standard-male geometry chunk byte for byte.

Rebuild the female reference from the [official HRA united-female v1.10 GLB](https://purl.humanatlas.io/ref-organ/united-female/v1.10) and [v1.5 GLB](https://purl.humanatlas.io/ref-organ/united-female/v1.5) with `npm run import:primary-female`. The importer verifies both GLB hashes, reads node names, hierarchy, transforms, ontology IDs, and material colors, then simplifies each mesh within a 0.2% relative error bound. Both GLBs contain no image textures. Put the [University of Denver original female STL ZIP](https://digitalcommons.du.edu/visiblehuman/1/) at `.local-models/source/Final 3D STL Models-stl-female.zip`, install NumPy for your Python interpreter, then run `npm run import:primary-female-muscles`. The second script verifies the ZIP hash, uses its original femur, tibia, and fibula STLs to fit each leg to HRA, and adds 74 nonduplicate original muscle meshes. The source ZIP stays outside `public` and the transform and bone-fit measurements are recorded in the manifest. No intermediary adaptation is used.

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

The server provides `search_anatomy` for names and IDs and `show_anatomy` for an interactive view. For example, call `show_anatomy` with `{"structure":"Stomach","model":"male-detail","controls":["systems","open"]}` to show Systems without Explode. Omit `controls` for the default interface, or pass `[]` for a bare viewer. To select a precise side or variant, use the ID returned by `search_anatomy`, such as `ZA:Sternocostal head of pectoralis major muscle.l`. The view tool returns a deployed URL and copyable iframe code. MCP Apps-capable clients can display the interactive viewer inline; other clients can open the URL or use the iframe HTML. The MCP process runs locally; GitHub Pages hosts the viewer only.

Run `npm run test:mcp` to verify tool calls, URL selection, and the UI resource.

## License

The application code is distributed under the [GNU General Public License, version 3](LICENSE) (GPL-3.0-only).

Issues and pull requests are welcome. Please include reproduction steps and browser/device details for interaction problems.

The selected detailed male and female models use the direct-source importers above. The former intermediary-derived binaries, manifests, and import scripts have been removed from this checkout.

### Display and model quality

The interface and 3D stage follow the operating system’s light/dark preference, including changes while the viewer is open. Wheel zoom follows the pointer; selecting a structure makes it the orbit pivot, including when other anatomy remains visible.

**Male · detailed** and **Male · full resolution + skin** use the same full-triangle direct Z-Anatomy build with original Open 3D Model upper-limb additions. The full resolution option initially shows the original source's segmented outer body surface. The Systems panel's skin opacity slider can fade or remove it; this surface is not a histological skin-layer model. The original publisher material colors are retained in the manifest, and the viewer uses its own calmer anatomy palette.

### Dissection, selection sets, and view URLs

Open **Study** to choose body regions and camera directions, keep multiple structures solid against adjustable transparent context, peel system layers, or move axial/sagittal/coronal clipping planes. **Add to selection set** works in search; Shift-click adds meshes directly. **Focus selection** zooms in while retaining context. **Isolate selection set** removes all context. Undo/redo covers dissection and selection changes; continuous slider changes are grouped. Selected pieces survive broad layer peeling, but remain subject to cross-sections.

Labels attach to actual triangle centroids, moving with the model. Dashed leaders and “behind” identify anchors obscured by selected geometry. Up to 20 labels are drawn simultaneously; every selected item remains listed in Study. Cuts expose open mesh boundaries without synthesizing tissue interiors. Regional filtering uses source memberships where provided and spatial bounds for older datasets. The whole-body surface is hidden in regional views.

The platform and floor have been removed. Zoom is cursor-centered, supports close inspection, and limits zoom-out relative to the current fitted view.

Use **Share view URL** to get a URL containing the model, selected IDs, camera position/target/projection offset, visible layers, region, opacity, hidden pieces, peeling and section. The URL can also be authored directly, for example:

`http://localhost:3016/?model=male-detail&select=Stomach&select=Pancreas&view=back&context=0.15`

`select` accepts a concept ID, source ID, or exact display name and can be repeated. `view` accepts `three-quarter`, `front`, `back`, `side` (patient-left), `right`, `superior`, and `inferior`. `focus=1` frames a selection when a camera is not supplied. `cut=axial|sagittal|coronal` and `slice=0..1` control sections; `flip=1` reverses the retained side. Unknown names are ignored. URLs are self-contained and do not upload views to a server.

Direct-source coverage is checked by `scripts/validate-primary-atlas.mjs` and `scripts/validate-primary-female.mjs`. They verify source hashes, part and concept identities, indices, bounds, and every raw/gzip chunk. Meshes absent from the publishers' files are not fabricated.

Additional checks: `node --experimental-strip-types scripts/validate-study.mjs` and `node --experimental-strip-types scripts/validate-camera.mjs`.
