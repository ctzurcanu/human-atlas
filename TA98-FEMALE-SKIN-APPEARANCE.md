# Female skin and scalp-hair appearance

Current appearance: UV-mapped thick 1 mm brush stubble, with actual scalp skin visible between hairs. No skin or facial geometry changes. Preserve all native component ownership.

One existing scalp-hair source mesh, kept separate from skin. Skin may select its TA98 hair child as a surface appearance component.

| TA98 | Parent | Name | Representation |
|---|---|---|---|
| A15.2.07.037 | A16.0.00.014 | eyelashes | missing distinct anatomical representation |
| A16.0.00.014 | A16.0.00.002 | hairs | visible scalp-hair envelope only |
| A16.0.00.015 | A16.0.00.014 | downy hair | missing distinct anatomical representation |
| A16.0.00.016 | A16.0.00.014 | hairs of head | missing distinct anatomical representation |
| A16.0.00.017 | A16.0.00.014 | eyebrows | missing distinct anatomical representation |
| A16.0.00.018 | A16.0.00.014 | beard | missing distinct anatomical representation |
| A16.0.00.019 | A16.0.00.014 | hairs of tragus | missing distinct anatomical representation |
| A16.0.00.020 | A16.0.00.014 | hairs of vestibule of nose | missing distinct anatomical representation |
| A16.0.00.021 | A16.0.00.014 | axillary hairs | missing distinct anatomical representation |
| A16.0.00.022 | A16.0.00.014 | pubic hairs | missing distinct anatomical representation |
| A16.0.00.024 | A16.0.00.014 | arrector muscle of hair | missing distinct anatomical representation |
| A16.0.00.025 | A16.0.00.014 | hair streams | missing distinct anatomical representation |
| A16.0.00.026 | A16.0.00.014 | hair whorls | missing distinct anatomical representation |
| A16.0.00.027 | A16.0.00.014 | hair crosses | missing distinct anatomical representation |

Skin/body/head UVs and geometry remain unchanged. The brush-cut envelope replaces the rejected shortened-bob hairstyle. It uses new hair-only topology and repeat coordinates; this is a cosmetic envelope, not a complete hair anatomy repair. Native skin and facial geometry remain unchanged.

The generated open-eye head atlas and seamless body detail remain local. The original body anatomical markings remain in their UV positions. Skin hue changes are material appearance only. No hierarchy IDs are removed, no ignore rules changed and no public deployment performed.

## Brush cut replacement plan

The rejected shortened bob is replaced by a new cosmetic scalp envelope derived from the unchanged native head surface. Its smooth frontal/temporal/nape boundary excludes face, eyes and ears. Roots sit 0.4 mm outside the scalp; the brush tips lie 10 mm along the scalp normals. The root and tip surfaces join only at their common hairline. This is a visible scalp appearance, not individual strands or follicles. The same Hair ID and Skin membership remain; new hair indices and repeat coordinates replace the old hairstyle only. Head/body coordinates, normals, UVs and interfaces remain exact. Validate closure, one component, self-crossings and scalp intersections before local installation and inspect the rendered result.

Historical 1 cm cosmetic preview (superseded): `outputs/female-brush-hair.png`. Brush shell: 448 vertices, 892 triangles, one closed genus-zero component, zero self-crossings and zero scalp intersections. Runtime catalogue validation and TypeScript compilation passed. The generated texture is used through the built-in image generation mode, optimized to 512 × 512 JPEG; no strand geometry or new skin UV attributes are stored.

## Bristle appearance correction

Replace the rejected opaque 1 cm cap with a submillimeter scalp root envelope and 768 deterministic short crossed-card bristles. Each bristle starts at the actual scalp normal and ends approximately 1 cm outward; roots use the same scalp carrier. This is a visible coarse hairstyle, not microscopic hair anatomy. Cards are intentionally open appearance sheets. Keep the Skin/head mesh, hairline and all native face UVs exact; check the closed root shell separately and verify root/tip placement. Only the Hair geometry ID is replaced; old immutable buffers remain stored but unused. Eyelashes use separate native eyelid-rooted fibers, visible at Depth 0. Their follicles and detailed coats remain missing.

## Superseded UV-cap appearance

The 1 cm cap and the 1 mm crossed-card version are superseded. The active Hair ID uses a 1 mm offset scalp carrier with an authored spherical UV chart, seam-corner duplicates and a transparent 512-pixel stubble tile. Texture gaps expose the unchanged native head skin; no skin color is painted into the hair texture. Coverage dithering preserves scalp gaps when the tiny stubs minify, rather than losing the whole crown through a hard mipmap cutoff. CPU picking samples a cached alpha mask so transparent gaps select the skin beneath. The old immutable hair buffers remain stored but are not active.

The active Skin/head source IDs and geometry are unchanged. Separate head/body material gains and a smooth neck calibration match the face/body transition. Native eyelid UV openings are narrowed to half-width 0.017 and half-height 0.0072; original painted eye remnants around the openings use nearby skin samples. Transparent holes expose two compact original ocular exteriors sharing a 20,196-byte buffer (960 triangles each), with procedural surface appearance and no ocular UV/image files. Their complete anatomical coats, optic-nerve interfaces and clinical packing are unaccepted; see the 131-term eye inventory.

Per the latest request, the small eyelash fibers are hidden by default and removed from the Skin assembly; they remain stored for recovery. Internal checked anatomy now remains present beneath opaque Skin at Depth 0 in male and female models. Depth testing controls visibility. This removes the prior unrequested opaque-skin shortcut rather than concealing source placement defects.

Eye-opening/visibility, runtime byte/gzip, TypeScript, camera-intent, explosion-plane and explosion-viewport checks pass. All appearance changes remain local; original assets, ignore rules and source anatomy are preserved.

Image mode: built-in image generation. Final UV tile: `.local-models/ta98-review/female-skin-appearance/hair-stubble-uv-v2.png`. Full prompt in `reports/ta98-female/skin-appearance-prompts.md`. Preview: `outputs/female-uv-stubble-final.png`.

## Current measured brush hair

The oversized transparent UV marks and offset cap are superseded. The same Hair ID now displays 56,597 instances of one closed six-sided shaft. Each shaft is 1 mm long and 0.07 mm in diameter. Area-weighted native scalp roots have mean nearest-neighbor spacing 0.694 mm (median 0.678 mm); this is approximate spacing, not a regular tiled grid. Root normals set the local shaft direction. Circumferential shaft U and root-to-tip V map a procedural fiber appearance independently of the scalp chart, so scalp UV stretch cannot enlarge a hair. Root UVs follow the carrier for subtle color variation. Actual native skin remains visible between shafts. No bitmap cap is drawn.

The root field uses 1,471,524 shared bytes; there is one prototype rather than 56,597 separately stored meshes or transform matrices. The closed prototype has 24 triangles per shaft; total submitted hair triangles are approximately 1.36 million. This increases GPU geometry work compared with the rejected flat cap. It is not a complete anatomical follicle representation. Old hair assets stay preserved and unused. Installation is local in all three female catalogues and the female rebuild hook reapplies it after appearance installation. Closed-shaft dimensions/UV/topology, runtime integrity, rendered close-up and camera/explosion invariants are checked before completion.


## Revised female density and scalp extent

Latest requested dimensions: shaft length remains 1 mm, diameter becomes 0.1 mm, and mean neighboring-root spacing targets 0.6 mm. Rebuild the native scalp carrier with a lower occipital/nape boundary (1.601 m rather than 1.621 m), a smooth temporal descent, the original frontal hairline, and a lateral ear exclusion guard. Crown, temples, above/behind ears and nape belong to the one Hair ID; face and auricles remain native skin. The original scalp/head vertices, normals, UVs and anatomy remain unchanged. The carrier is a hidden picking/support surface, not a displayed textured cap. Every displayed shaft uses its own cylindrical UV chart and the shared root field. Keep all older metric fields/carriers preserved. Validate carrier topology/self/skin contacts, field dimensions/spacing, and actual anterior/lateral/posterior renders before accepting the appearance.


The user requested a further nape descent. The occipital boundary is now 1.566 m with a smooth rise to the frontal hairline at 1.666 m; ear exclusions remain. Requested diameter 0.1 mm, mean pitch 0.6 mm and length 1 mm are retained. Eye apertures are narrowed to half-width 0.0155 and half-height 0.0056, with tapered eyelid corners, shallow lower contour and mirrored corner tilt; CPU selection and GPU visibility use the same mask. Duplicate added globe parts are retired; the unchanged original native eye exterior/cap surfaces are the visible pair. Earlier appearance sections are historical and superseded by these latest changes.


Final current metric field: 98,762 shared shafts, 0.600536 mm mean nearest-root spacing, 0.1 mm diameter, 1 mm length. Lower nape coverage includes an auricular-base exclusion to keep the ears clear; the earlier lower trial without that exclusion was rejected for self/skin crossings. The accepted carrier has zero detected self/skin crossings. Native eyes remain the sole pair and smaller eyelid masks remain active. Four targeted dimension/visibility tests, runtime preservation, TypeScript and camera/explosion checks pass. Actual posterior, lateral and frontal appearance is reviewed locally; this is cosmetic appearance evidence, not full anatomical acceptance of the head.


The carrier-edge hairstyle is superseded by direct root sampling on the unchanged native head triangles. Root acceptance evaluates a continuous smoothstep curve with curved frontal/nape profiles and smooth ear exclusion, rather than the clipped triangle-edge polyline. Current female field has 98,370 roots, 0.598854 mm mean spacing, unchanged 1 mm length and 0.1 mm diameter. Native skin source hashes remain exact. Live final curved-edge review remains pending because the browser/computer-use previews became unavailable during the male/female appearance turn.


Preview access recovered. Actual separate frontal/lateral review confirms the continuous female hairline and preserved native eyes. The user-controlled prior view was not moved; it is cropped behind the inspector, so a separate temporary review was used. Zero pending chunks and no browser errors. Images: outputs/female-curved-brush-front.png and female-curved-brush-side.png. Cosmetic appearance verification does not establish full scalp/head anatomical completion.
