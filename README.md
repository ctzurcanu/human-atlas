# Human Atlas

An interactive 3D anatomy explorer built with React, Three.js, and targeted publisher geometry: **5,209 selectable surfaces**, **4,669 named concepts**, and **10,391,079 triangles** in the default **Male · detailed** model. The standard BodyParts3D male and partial female references are also available. Confirmed duplicate imports are excluded for bowel, pancreatic duct, deferent ducts, dorsal penile veins, hand and wrist bones, and selected upper-body vessels.

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

### Browser data cache

Downloaded model catalogues, compressed geometry (including VR), textures, built-in hierarchies and transcript branches are reused from IndexedDB on later visits. The small file-version index is stored in localStorage; model files exceed localStorage's approximately 5 MiB capacity. Each build publishes `asset-manifest.json` with content hashes, so changed files refresh automatically while unchanged files remain cached. Local development uses file sizes and modification times. The cache is limited to 512 MiB and evicts the least recently used files; browsers can also reclaim this storage. Blocked or full storage falls back to ordinary downloads. Connect traffic and external hierarchy URLs are not cached.

Add **`?rc=1`** to the viewer URL to reset the cache, or **`&rc=1`** when there are already query parameters, for example `http://localhost:3016/?model=cell&rc=1`. The flag is consumed once before loading; other query parameters and saved-view fragments remain intact. Only this atlas deployment's cached files and version index are cleared. Saved views, slides and connection preferences are preserved. Run `npm run test:cache` to verify persistence, file updates, reset, eviction, aborts, damaged files and storage failures.

### Connect: host a live view

In **Advanced tools → Connect**, choose **Host** to prepare a waiting lobby. The Copy icon beside Host copies the invitation URL. The next 150px row contains the host address, Play, and connected guest count. Localhost and loopback addresses cannot be shared: the viewer discovers a non-loopback LAN IP or the local ngrok HTTPS tunnel when available, and you can enter another IP or HTTPS address. Use the HTTPS/ngrok address when guests need secure WebSockets. In Guest mode, the next row contains the invitation URL input and Play. Pasting an invitation for another viewer opens that viewer and joins in one Play click. Guests connect and wait; pressing Host **Play** includes exactly the guests online at that moment, and guests arriving later wait for the next session. Names default to the IP seen by the server, with local connections using a non-loopback IP. You can edit waiting guest names before Play, and guests can keep their own labels for previous hosts.

During a session, the host's model, selections, layers, sections, hierarchy, camera movement, exploded-piece rotation, and visible viewer panels are shared. Guests follow the host without changing the host's view. Collapsed Connect shows a hollow button containing **H** or **G**. **H** reopens Advanced tools while the broadcast continues, so the host can choose Slides or Quiz. In the expanded Connect panel, the actual **Stop** replaces Play beside the host address. The panel lists session guests with a Remove button for each; removing one guest keeps everyone else connected, and the removed identity cannot rejoin until that session ends. Later arrivals are marked Waiting. A guest's **G** leaves the connection; the host's actual **Stop** ends it for every guest and returns the host to the setup panel. The host invitation remains available for the next lobby. Closing a tools panel keeps the connection running. Brief connection interruptions retry automatically; a host missing for 60 seconds ends the session.

Connect's relay uses `/atlas-connect`. HTTPS viewers use **secure WebSockets (`wss://`)**; localhost development uses its local relay. For a built viewer with the relay, run `npm run build` followed by `npm start` behind HTTPS. GitHub Pages serves the viewer but cannot run a WebSocket server. On GitHub Pages, enter the HTTPS address of a running atlas relay in the Host field, then leave the field to prepare the lobby. Invitations keep the GitHub Pages viewer URL and include the relay address, so guests use the same server automatically. The relay permits this project's GitHub Pages origin and its own origin; other browser origins are rejected. Ngrok is one way to expose the relay, but any HTTPS server running it works. Set the repository variable `VITE_CONNECT_RELAY_URL` for a default relay in Pages builds, or open the viewer with `?relay=https%3A%2F%2Fyour-relay.example`. The manually entered address is remembered on GitHub Pages. Initial failures stop retrying after 10 seconds with an actionable error. Session invitations are private, random codes, and only their host can broadcast state. Run `npm run test:connect` to exercise WSS sessions, roles, joining, reconnection, names, and Stop behavior.

### Cloudflare Connect hosting

For a hosted relay without ngrok, follow [Cloudflare Connect setup](cloudflare/SETUP.md). The deployed relay is `https://human-atlas-connect.ctzurcanu.workers.dev`, already configured as the GitHub Pages viewer's default in `app/connect-config.ts`; publish the updated viewer to activate it. The Worker uses SQLite-backed Durable Objects and idle WebSocket hibernation. Run `npm run connect:cloudflare:check` before `npm run connect:cloudflare:deploy`.

### Portable view framing

New saved links and slides use `frame=` with dimensionless coordinates: direction, a target relative to the displayed anatomy bounds, an anchor within the usable viewing area, relative distance, and up direction. Live Connect snapshots and camera messages use the same frame. Screen anchors and direction components are 0–1 fractions; target fractions can extend beyond that range when panning outside the anatomy. The local viewport, visible panels, advanced tools, and aspect ratio determine the actual projection on each device. The visual center is kept above half the viewport height, with advanced tools reserved below it. Resize and panel changes preserve the current normalized framing. Legacy `camera=` links remain readable and adapt using a canonical 16:9 source projection because their original screen dimensions were never saved.

Supported HTTPS WebXR browsers expose VR and AR entry buttons. AR uses passthrough and optional DOM overlay controls. XR presentation placement uses the device's projection while the headset or phone continues to control tracked movement. New iframe snippets include `xr-spatial-tracking`. Physical device validation is still required for headset and AR browser behavior. `npm run test:framing` checks framing across aspect ratios, anatomy scales, panel reservations, axial orientation, and XR transform behavior.

### Quest VR

Open the viewer over HTTPS in Meta Quest Browser, wait for **Preparing VR…** to become **Enter VR**, then enter the headset view. Triggers select structures and the left stick moves closer or farther.

Quest uses separate simplified geometry, one anatomy draw per eye, Lambert lighting, and hashed opacity. Hidden and fully transparent parts are excluded from the draw. Source textures and detailed filled section caps are omitted in this mode; section clipping remains available. The detailed male VR model has 387,438 triangles and retains all 5,208 part IDs. Desktop, mobile, and other browsers keep the original model and materials, and leaving VR restores the original view.

After changing a public anatomy model, regenerate its VR assets with `npm run build:vr` and check them with `npm run test:vr`. The generated `.vr.json` and `.vr.bin.gz` assets are published alongside the source models. For local development, `?vr-preview=1` renders the Quest profile without a headset; this preview flag is disabled in production. Browser preview and geometry tests do not measure physical headset frame rates.

### Local Anatomy Atlas GLB model


Systems, Regions, and Depth use one shared set of intermediate branches across models. Empty branches disappear. Verified anatomical parents follow TA98; model subdivisions and origin/insertion markers retain their own atlas IDs beneath verified TA98 parents. Unresolved structures remain available in their source and regional browsing groups. Source region groups add finer skin branches when available; mammary tissue, subcutaneous fat, and lacrimal glands stay outside Skin. Run `npm run test:hierarchies` to check that every available model part appears exactly once in each hierarchy. [TA98-HIERARCHY-AUDIT.md](TA98-HIERARCHY-AUDIT.md) lists every unresolved mapping, parent-only placement, cellular component, and suppressed mesh. Rebuild the concordance with `npm run build:terminology` and regenerate the report with `npm run audit:ta98`.

The download and converted data stay in `.local-models/`, which is ignored by Git and excluded from `dist/`. 

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

### Additional hierarchies

The final **v** menu in Layers lists **Genes**, **Cell Types**, **Physiology**, **Dermatomes and Myotomes**, **Drugs**, **Physical Exercise**, then **Chakras**. Genes retain the full HGNC family hierarchy; Cell Types retain the human Cell Ontology subtype hierarchy with a separate developmental Lineage view. Both offer search, reciprocal marker links, and selectable modeled anatomy leaves. Mapped genes can load their transcript/protein branches on demand. See [BIOLOGICAL-HIERARCHIES.md](BIOLOGICAL-HIERARCHIES.md) for sources, scope, counts and review links. Run `npm run build:biology` and `npm run test:biology` to rebuild and validate the pinned biological data. **Physiology** organizes documented structure–function associations into nested processes, with links to relevant Cell Types and Genes. It appears immediately after Cell Types and loads only when selected. See [PHYSIOLOGY-HIERARCHY.md](PHYSIOLOGY-HIERARCHY.md) for scope, function links, geometry coverage and rebuild commands; all scientific references are in [ATTRIBUTION.md](ATTRIBUTION.md). **Dermatomes and Myotomes** contains all 31 spinal levels with left/right branches, source-linked sensory territories, and representative overlapping muscle associations. Root geometry varies by model; true dermatome skin patches are not present in the source models. See [DERMATOMES-MYOTOMES.md](DERMATOMES-MYOTOMES.md) for linked levels, source references, and geometry gaps. Run `npm run build:spinal-roots` and `npm run test:spinal-roots` to regenerate and validate this hierarchy. It loads only when selected. **Drugs** retains all five ATC levels supplied by ChEMBL, drug records, curated mechanisms and human target genes, with links to Genes, Cell Types and available anatomy. It loads only when selected. See [DRUGS-HIERARCHY.md](DRUGS-HIERARCHY.md) for scope and build commands; credits are in [ATTRIBUTION.md](ATTRIBUTION.md). **Physical Exercise** groups source exercises and variants by movement type, with primary/secondary muscle roles and documented joint, tendon, nerve, respiratory, cardiac, vestibular and pelvic-organ associations. It loads only when selected. See [PHYSICAL-EXERCISE-HIERARCHY.md](PHYSICAL-EXERCISE-HIERARCHY.md) for scope, linked review gaps and rebuild commands; all source credits are in [ATTRIBUTION.md](ATTRIBUTION.md). **Chakras** is based on the [source outline](https://hackmd.io/21uN7y70SPOLKgK_mYZ-mA). Its **Other** branch contains every model part not mapped to a chakra, grouped by system. The deployable JSON file is [public/assets/chakras.json](public/assets/chakras.json); `npm run build` copies it to `dist/assets/chakras.json`. You can load another hierarchy from a URL in the same menu. Its server must allow the viewer's origin to fetch the JSON. Share links retain custom guest URLs, and Explode follows the selected guest tree.

A guest file uses `{"schema":"human-atlas-hierarchy/v1","id":"example","name":"Example","nodes":[{"name":"Group","children":[{"name":"Heart","matches":["Heart"]}]},{"name":"Other","unmapped":true}]}`. `partIds` in a v2 graph selects exact meshes without expanding their parent concepts; `matches` contains exact atlas concept names or IDs; `matchBases` contains exact concept names ignoring a trailing `(left)` or `(right)`; `systems` may group all meshes in a source system. The optional top-level `unmapped` branch collects remaining model parts. Nodes without matching anatomy remain visible as labels but have no checkbox. Run `npm run build:chakras` after updating anatomy catalogues or the TA98 concordance. This expands [scripts/chakras-outline.json](scripts/chakras-outline.json) with spinal/vagus branches, modeled innervation targets, and the requested system associations. The derivation runs offline, so opening Chakras in VR only resolves name lookups. Shared structures may belong to multiple chakras; Explode still assigns each mesh once. Run `npm run test:chakras` to check all available catalogues, all 31 spinal levels bilaterally, innervation membership, and the Explode layout. See [CHAKRA-HIERARCHY.md](CHAKRA-HIERARCHY.md) for the mapping rules.

The detailed male reference retains each original Z-Anatomy object and its material-defined surfaces, including tendon and cartilage patches. Original Open 3D Model upper-limb files contribute additional nerves, vessels, muscles, hand bones, ligaments, and sheaths; the right-side additions are mirrored to the left and labeled as such in the manifest. Blender reference lettering is excluded. Muscle attachments and fascia have separate layers, hidden by default. Counts represent selectable surfaces, not unique organs; coverage and granularity differ between models. The standard male reference remains selectable.

The original Z-Anatomy and Open 3D Model archives and source files are SHA-256 verified. The direct importer bakes source transforms, converts Z-up to the viewer's Y-up coordinates, splits material surfaces, and packs raw and gzip buffers without simplifying the source triangles. Original material colors are recorded as `sourceColor`. The Z-Anatomy Blender archive contains no anatomical image textures, and its packaged atlas has no muscle UV maps; the viewer projects the verified Open 3D Model muscle texture onto those muscle surfaces with a three-axis mapping. This projection is an approximate visual treatment, not a source-authored Z-Anatomy texture. For Open 3D Model surfaces, the importer preserves UV seams and the embedded base-color and normal textures used by the included meshes, including the original muscle tiles. Those texture files are copied directly from the verified publisher GLB and their hashes are recorded in the manifest. Two hand ligaments in the source lack UV maps and remain palette-colored.

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

The server provides `get_anatomy_options` for valid model, system, depth, region, and control IDs; `search_anatomy` for names and IDs; and `show_anatomy` for an interactive view. For example, `{"structure":"Stomach","model":"male-detail","hierarchy":"depth","depthHidden":["skin"],"explode":0.5,"controls":["systems","explode","open"]}` selects the stomach in Depth and opens the hierarchy halfway. Omit `structure` to show the whole model, or pass `structures` to select several concepts. The view tool can also set visible systems, hidden pieces, region, camera, labels, skin opacity, isolation, and one or two sections. It returns a deployed URL and copyable iframe code. MCP Apps-capable clients can display the interactive viewer inline; other clients can open the URL or use the iframe HTML. The MCP process runs locally; GitHub Pages hosts the viewer only.

In browsers that support WebMCP, the open viewer also exposes `get_anatomy_view`, `set_anatomy_view`, `act_on_anatomy_view`, and tools to list, save, open, or delete browser-local views. These live tools can change the current selection, Systems/Regions/Depth tab, visible systems and depth layers, hidden pieces, Explode, cuts, camera, opacity, labels, isolation, and rotation. Actions include undo, redo, reset, focus, hiding or clearing the selection, and PNG download. Browser-local saved views and PNG capture require the open browser; the standalone MCP server only generates shareable views.

Run `npm run test:mcp` to verify tool calls, URL selection, and the UI resource.

## License

The application code is distributed under the [GNU General Public License, version 3](LICENSE) (GPL-3.0-only).

Issues and pull requests are welcome. Please include reproduction steps and browser/device details for interaction problems.

The selected detailed male and female models use the direct-source importers above. The former intermediary-derived binaries, manifests, and import scripts have been removed from this checkout.

### Display and model quality

The interface and 3D stage follow the operating system’s light/dark preference, including changes while the viewer is open. Wheel zoom follows the pointer; selecting a structure makes it the orbit pivot, including when other anatomy remains visible.

**Male · detailed** and **Male · full resolution + skin** use the same full-triangle direct Z-Anatomy build with original Open 3D Model upper-limb additions. The full resolution option initially shows the original source's segmented outer body surface. The Systems panel's skin opacity slider can fade or remove it; this surface is not a histological skin-layer model. The original publisher material colors are retained in the manifest, and the viewer uses its own calmer anatomy palette.

### Dissection, selection sets, and view URLs

Open **Study** to choose body regions and camera directions, keep multiple structures solid against adjustable transparent context, or peel system layers. Open **Sections** for a skin preview with section planes and one or two tabs. The default first cut is axial at 38%. Select a tab to edit its cut; enabled cuts apply together to the same checked anatomy. Place the active cut on the main model or move its position slider, and choose an axial, sagittal, coronal, or oblique plane. Turn and Tilt set an oblique angle. Scroll over the model to move the active cut; Ctrl-scroll zooms. Closed mesh intersections receive tissue-colored caps; open surfaces and thin coverings receive narrow cut bands. Skin cuts have a narrow thickness, and bone cuts show a cortical rim around a seamless illustrated cancellous-bone texture. The bone texture is a generated visual aid, not an image of the source bone or a histology scan. Every cut face follows the same selection, context, skin, and source opacity as its mesh. The detailed male stomach uses the aligned BodyParts3D publisher surface because the Z-Anatomy stomach has a large opening in its anterior wall. **Add to selection set** works in search; Shift-click adds meshes directly. **Focus selection** zooms in while retaining context. **Isolate selection set** removes all context. Undo/redo covers dissection and selection changes; continuous slider changes are grouped. Selected pieces survive broad layer peeling, but remain subject to cross-sections.

Cut faces use exactly tileable illustrated images in `public/models`; these depict tissue appearance but are not volumetric scans of the source anatomy. `npm run build:tissue-textures` regenerates the muscle surface and the axial, coronal, and sagittal muscle cuts plus pink lung, liver, yellow spleen, CNS, and general organ cuts. Bone retains its cancellous-bone tile and adipose its own tile. Untextured muscle exteriors receive a red muscle surface rather than source overlay colors. The spleen has a separate warm yellow treatment instead of the lymphatic group's green. Digestive canals, heart chambers, nasal and buccal cavities, trachea, and vessels show a wall cut with an empty lumen; solid tissues retain a filled cut. Exterior source textures remain intact where supplied. Reference GLB material primitives of one muscle or lung lobe are combined before cut geometry is calculated, so shared material boundaries do not become false open cuts. Caps are added as model chunks load and refreshed when a concept's later material pieces arrive. `npm run audit:meshes` writes a per-model source-geometry report under `reports/mesh-audit/`, covering boundaries, winding, duplicate faces, normals, missing mapped UVs, and coarse out-of-body bounds. Imports run the relevant audit automatically. The audit flags defects for review; it does not silently rewrite the publisher's geometry. `npm run test:sections` checks section geometry, and `npm run test:section-visuals -- <screenshot-directory>` compares the fixed browser views listed in `tests/section-visual/cases.json` against reviewed baselines at 1280×720. Capture each view after its model and caps have loaded before comparing. For one reference concept, `node scripts/inspect-cap-group.mjs 'REF:Rectus abdominis muscle.l' 1.05257` reports primitive and joined cut topology.

Labels attach to actual triangle centroids, moving with the model. Up to 20 labels are drawn simultaneously; every selected item remains listed in Study. Open meshes remain narrow bands rather than invented solid tissue. Regional filtering uses source memberships where provided and spatial bounds for older datasets. The whole-body surface is hidden in regional views.

The platform and floor have been removed. Zoom is cursor-centered, supports close inspection, and limits zoom-out relative to the current fitted view.

Use **Share view URL** to get a URL containing the model, selected IDs, camera position/target/projection offset, visible layers, region, opacity, hidden pieces, peeling and section. The URL can also be authored directly, for example:

`http://localhost:3016/?model=male-detail&select=Stomach&select=Pancreas&view=back&context=0.15`

`select` accepts a concept ID, source ID, or exact display name and can be repeated. `view` accepts `three-quarter`, `front`, `back`, `side` (patient-left), `right`, `superior`, and `inferior`. `focus=1` frames a selection when a camera is not supplied. `cut=axial|sagittal|coronal|oblique` and `slice=0..1` control the first section; oblique cuts also accept `azimuth=-180..180` and `elevation=-90..90` in degrees. `flip=1` reverses the retained side. Add `sections=2` and the corresponding `cut2`, `slice2`, `flip2`, `azimuth2`, and `elevation2` parameters for a second section; `sectionTab=2` selects it. Unknown names are ignored. URLs are self-contained and do not upload views to a server. New links generated by the viewer use compact mesh positions (`p`, `s`, `h`) and a compact layer mask (`l`); manually authored links and previously saved links with anatomical names still open. Dense hidden sets use ranges, while scattered sets use a bitset. A fingerprint prevents a short link from selecting unrelated anatomy if the model's part list changes; in that case its indexed selection and hidden set need to be recreated against the updated model. Use explicit `select` and `hide` IDs for links that must survive a model catalogue rebuild.

Direct-source coverage is checked by `scripts/validate-primary-atlas.mjs` and `scripts/validate-primary-female.mjs`. They verify source hashes, part and concept identities, indices, bounds, and every raw/gzip chunk. Meshes absent from the publishers' files are not fabricated.

Additional checks: `node --experimental-strip-types scripts/validate-study.mjs` and `node --experimental-strip-types scripts/validate-camera.mjs`.
