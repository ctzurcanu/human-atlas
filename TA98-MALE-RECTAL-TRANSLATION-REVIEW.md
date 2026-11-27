# Male rectal group: literal translation review

The current Blender trial uses the intact native distal shared wall, replacing
the earlier deformed review candidates. Its group moves 20 mm posteriorly and
8 mm superiorly, with zero lateral displacement, zero rotation and unit scale.
These distances are an initial placement trial, not accepted anatomical values.

`reports/ta98-male/male-rectal-translation-only-replay.json` proves that every
source NPZ array except the common origin is unchanged and that the wall binary
is byte-identical. Both complete current external anal sphincter groups receive
the same translation. Sacrum, coccyx, pubic symphysis, prostate and the three
bladder regions retain their original positions, indices and normals.

The focused outer-surface screen still finds 38 triangle intersections with the
sacrum, 49 with the prostate and 91 with the bladder fundus. These counts are
triangle-pair evidence, not numbers of separate anatomical defects. This trial
is not anatomically accepted or installed in the main male TA98 catalogue.
Remote sigmoid, pelvic-floor, perineal and neurovascular interfaces remain
unfinished; moving the two sphincter groups does not establish all attachments.

Local rendered review:
<http://localhost:3016/caecocolic-wall-review.html?assembly=male-native-distal&course=translated&context=pelvis>.
The initial view includes the complete wall and pelvic guides. Explicit
**Center pelvis** fits that same complete extent. Visibility and section changes
preserve the camera; this was verified on the rendered page by exact equality
of its recorded position, orientation, target and zoom before and after toggles.
The review route does not replace historical `course=smooth` geometry or modify
the main atlas. Blender renders are retained in
`outputs/male-rectal-translation-only-{lateral,sagittal}.png`.

The bounded rigid comparison in
`reports/ta98-male/male-rectal-rigid-translation-screen.json` tests twelve
placements against the same fixed six pelvic surfaces. The source has 466
triangle pairs; 20 mm back / 8 mm up has 178. At 28 mm back / 4 mm up the total
falls to 77, but sacral pairs rise from 38 to 65. At 32 mm back / 4 mm up,
coccygeal intersections appear as well. None of these trials clears every
surface. A lower aggregate count is therefore not sufficient to accept a move.

The contact-owner report identifies the actual source regions: at 20/8,
all 38 sacral and 91 bladder-fundus pairs involve the rectal region; the 49
prostate pairs involve the unclassified distal region. At 28/4 the sacral
contacts extend into the sigmoid region. The contact bounds and exact triangle
pairs are retained, so subsequent placement review can inspect these locations
without guessing from the organ bounding box or altering the source shape.
No lower-scoring trial has automatically replaced the Blender or main model.

## Installed source distinction

The current main catalogue uses active part `FJ2571` (BodyParts3D), with 1,538
triangles. Its forest owner contains that part alone; the native distal shared
wall is a separate, uninstalled alternative. The installed source is a single
closed genus-zero surface with no boundary loops and no modeled lumen. Its
current fixed-pelvis screen finds 1,287 sacral and 181 coccygeal triangle pairs,
with zero pairs against the three bladder regions and prostate. These results
must not be substituted for the native trial's results or described as two
active rectums.

`reports/ta98-male/installed-rectum-pelvis-screen.json` binds this evidence to
the current catalogue. `scripts/blender-review-installed-male-rectum.py`
renders it with the exact current fixed pelvic guides, one rectum only, in a
separate Blender scene. The actual sagittal rendering visibly confirms the
sacral overlap. No source positions or catalogue data change during that review.
This prevents further fitting work from confusing the installed BodyParts3D
object with the native wall trial or promoting a closed-solid source as a
completed hollow rectal wall.

## Sphincter placement correction under review

The user's observation of the displaced sphincter is confirmed. Translating
the unrelated source muscles with the tube did not establish an attachment.
The geometric rim/contour anchor trial requires a common relative muscle
translation of `[0, -4.928462, 40.664617]` mm before the existing rectal-group
translation. It applies no rotation, scale or source-mesh edits. The anchor
aligns the muscle's lower extent with the literal distal rim height and its
sagittal centre with the actual tube section at the muscle's middle height.
These are provisional geometric anchors, not accepted anal-verge landmarks.

The resulting source shapes still cross the tube: 94 left and 101 right
triangle pairs. This prevents treating the centred position as a completed
sphincter fit or installing it in the main atlas. The render at
`outputs/translated-distal-sphincter-alignment.png` shows a separate Blender
alignment scene that shares all original mesh datablocks, adds no mesh vertices
or triangles, and preserves the previous review scene and main catalogue.
Evidence: `reports/ta98-male/translated-distal-sphincter-alignment.json` and
`translated-distal-sphincter-alignment-render.json`.

## Supplied source-reference review

The supplied 2024 source plan contains a six-piece male and three-piece female
external sphincter assembly, each alongside a ten-piece anal-canal assembly.
These counts are source pieces, not proof of correspondence to the named deep,
superficial and subcutaneous TA98 portions. Assembly candidates are kept
distinct from the two alternative individually named canal source structures.

`scripts/blender-review-sphincter-source-reference.py` reads the exact SHA-bound
source geometry into separate reference-only Blender scenes. Both lateral and
sagittal source-relative views were rendered and inspected. The source muscles
surround a substantial length of the canal wall; the current floating two-piece
muscle and its simple bounds-based alignment do not reproduce that assembly
relationship. No source geometry is installed, redistributed, or used as a
replacement main mesh. The reference scene does not superimpose another pelvis
on the current body.

Evidence is in `reports/ta98-male/sphincter-canal-source-reference.json`, with
source hashes, every selected node, strip-degenerate counts, source-relative
bounds and render hashes for both sexes. Clinical metric calibration,
registration and exact muscle-subdivision identities remain unverified.
Reference images:
`outputs/{male,female}-sphincter-canal-reference-{lateral,sagittal}.png`.

## Original coarse cuff around the distal carrier

The construction uses the already inventoried whole external muscle, with the
named subcutaneous/superficial/deep portions still missing. Its inner surface
encloses actual sections of the unchanged distal carrier; its outer surface
adds finite illustrative muscle thickness. Annular rims join those surfaces
at both ends without a disc across the passage. Upper/lower allocation and
puborectal, perineal, skin and internal-sphincter interfaces remain provisional.
The whole cuff follows the same untranslated-source/shared-group coordinate
convention as the rectal wall, with no rotation of that group.

`scripts/construct-male-distal-muscle-cuff.py` creates an original 32-sided,
17-ring envelope. Current source-muscle height and one sparse relative-width
measurement from the supplied male reference inform illustrative dimensions;
no reference vertex, face, UV or image is copied. This produces 1,088 vertices
and 2,176 triangles, with no UV data or texture images. Appearance uses a
procedural circumferential muscle field.

The complete cuff is one closed genus-one component, with no open/nonmanifold
edges, duplicate/degenerate faces, self-intersections or tube intersections.
Independent sixteen-section replay confirms that the actual carrier remains
inside the cuff's void throughout the modeled interval. These checks establish
a usable coarse construction, not clinical EAS anatomy. Fixed-pelvis checks
still find 164 cuff/prostate triangle pairs at the present carrier placement,
with zero against sacrum, coccyx and the three bladder regions. It is therefore
not installed or anatomically accepted; the carrier's pelvic fit remains open.

Local preview:
<http://localhost:3016/caecocolic-wall-review.html?assembly=male-native-distal&course=translated&context=pelvis&muscle=cuff>.
The displaced two-part source alternative is omitted from this preview. The
cuff has its own visibility checkbox and uses the same live, tissue-only
section-cap bank as the seven wall intervals. Its visibility, section position
and section enablement preserve the recorded camera exactly. Actual Blender
lateral/sagittal renders and the browser preview were inspected; browser error
logs are empty. Evidence is in `male-distal-muscle-cuff{,-replay,-render}.json`
under `reports/ta98-male/`. The original review and main catalogues remain intact.

## Coupled local geometry correction plan

The original zero-rotation group pose remains the starting point. Further
correction explicitly changes the mesh: a shared smooth sagittal field carries
all eight wall interfaces and the cuff. Lateral and vertical coordinates,
face ownership, shared source controls, the skeleton and urinary structures
stay fixed. The full inventoried distal carrier is repaired together; no wall
layer, opening, muscle subdivision or boundary is represented by a duplicate
outer surface. The cuff remains a coarse whole-muscle construction.

The first field clears sacral, prostatic and bladder contacts, but creates 120
posterior rectal/coccygeal triangle pairs between heights 0.8915 and 0.9226 m.
The next trial preserves the anterior contour at those heights and narrows the
local front-to-back profile, following actual carrier sections. Its sectionwise
map is injective analytically, but discrete topology, all finite tissue cells,
all eight interface pairings, cuff void, fixed-pelvis contacts, rendered
sections and full-body packing remain required checks. Finite tissue rays may
need local direction correction after the complete exterior has changed.
No extra source vertices, source faces, UVs or texture images are planned.

This numerical fitting plan does not establish clinical rectal diameter,
mesorectal clearance or EAS subdivisions. Proximal digestive junctions,
pelvic-floor/perineal/skin insertions and neurovascular attachments remain
unaccepted. No candidate replaces a main catalogue without those remaining
assembly checks and anatomical review.

The corrected candidate passes all 36 wall-interface pair checks and seven
closed finite-shell checks, with one genus-one component per shell. All sixteen
checked tube sections stay inside the cuff void. Wall and cuff each have zero
surface intersections with sacrum, coccyx, prostate, bladder body, bladder fundus
and bladder neck. No source faces, ownership labels, controls, UVs or images
were added. The exterior's lateral/vertical coordinates are unchanged; local
inner tissue-ray direction corrections are explicitly recorded and do change
some inner coordinates.

Actual Blender and browser views were inspected. The result has a visibly
narrower front-to-back profile near the coccyx; its clinical shape and diameters
remain unaccepted. These source-body surface checks are not full-body packing
or evidence that remote interfaces fit. Section enablement and movement preserve
the browser camera exactly, with no browser error logs. The main models and
the previous straight-translation candidate remain unchanged.

Corrected local preview:
<http://localhost:3016/caecocolic-wall-review.html?assembly=male-native-distal&course=corrected&context=pelvis>.
Evidence: `coupled-distal-interface-correction{,-replay,-render}.json` under
`reports/ta98-male/`, plus the corresponding source-bound construction,
validation, export and Blender-review scripts. This is an explicit geometry
correction, not a claim that the original mesh was merely translated.

### Original source sphincters at the cuff centre (user-directed placement)

The original left and right external anal sphincter meshes share one translation:
−0.157 mm lateral, +5.572 mm superior, +6.890 mm anterior. Their combined
bounding-box centre now coincides with the corrected coarse cuff centre. Source
faces, normals, scale, orientation and relative placement are preserved. The
corrected review displays these two source muscles by default; the coarse cuff
is an unchecked comparison overlay. Regenerating the corrected review invokes
this placement explicitly. The main atlas remains unchanged.

This establishes the requested position, not anatomical acceptance of the fit:
104 left and 111 right triangle pairs intersect the current provisional tube.
The broad wall correction also failed the full current-body packing screen
(7,339 → 13,074 triangle pairs; 23 newly contacted neighbours), so it must not be
promoted to the main model. The localized correction remains under validation.
Evidence: `reports/ta98-male/original-sphincters-at-cuff.json`.

### External sphincter fit supersedes centre-only placement

The user's next instruction requires both source halves outside the canal.
Centre-only placement did not satisfy that requirement. Bounded canal-deformation
trials eliminated muscle intersections but inverted wall prisms and were rejected.
A common fit of the two source muscles instead preserves their topology, centre,
normal thickness and bilateral relationship: radial scale 1.7 in their original
ring plane and about 11 degrees rotation towards the canal's transverse plane.
The rectal wall, skeleton and urinary organs are unchanged. Actual source-muscle
surface intersections are zero against the canal, sacral bone/cartilage, coccyx,
prostate and all three bladder parts. The muscles render opaque with procedural
muscle colour; the coarse cuff remains hidden. This corrects the local review's
external presentation, not clinically accepted dimensions or whole-body packing.
Evidence: `external-original-sphincters-fit.json`; local runtime screenshot
`outputs/external-sphincters-browser.png`. The source-bound export repeats this
fit after regeneration, rather than restoring the inadequate centre-only pose.

The localized wall replay completed 36 boundary, seven finite-shell and 16
cuff-void checks. Its current-body screen completed 430 parts in 315 complete
owner groups, but found 12 newly contacted neighbours and remains unaccepted.
The sphincter external fit is separate from that wall trial.

### Full fitted-sphincter neighborhood screening

The complete current dependency forest expands the fitted muscle pair's bounds
to 73 active parts in 56 complete source-owner groups. The source-bound screen
finds 493 → 771 triangle pairs and six newly contacted neighbours: both fasciae
latae, both glutei maximi and both iliococcygei. Pubo-analis contacts also increase.
Consequently the narrow pelvic pass does not establish acceptable packing and
this fitted assembly remains excluded from main-model installation. Triangle
pairs are not penetration volume and do not distinguish shared attachment tissue.
Related structures must be reviewed as whole groups rather than moved silently.

The local corrected review now exposes both native pubo-analis and both native
iliococcygeus meshes as optional context. Four complete owning groups retain
all current world coordinates, source normals and faces; the additional packed
payload is 142,668 bytes with no new UVs or images. All four controls default off,
and enabling them preserves the exact browser camera. The actual combined view
was inspected. Main male/female catalogues remain unchanged. Evidence:
`external-sphincter-current-body-packing.json`, `distal-pelvic-floor-context.json`.

### Shared-neighborhood fit comparison and review update

Six fits with zero canal and fixed pelvic-guide intersections were compared
against one common union scope: 87 parts in 63 complete owning groups. The
selected diagnostic improvement uses radial scale 1.6 and half the rotation
between the native ring plane and a transverse canal plane (about 22 degrees).
Its centre, source topology and normal thickness remain unchanged; the rectal
wall and skeleton remain fixed. Compared with the previous local ring, total
triangle pairs fall 771 → 610 and newly contacted neighbours fall six → four.
Both fascia-lata contacts disappear. Both glutei maximi and iliococcygei still
contact the fitted pair; skin, pubo-analis and pararectal-node contacts also
require interpretation. No anatomical or main-model acceptance is claimed.

The local corrected review is updated without adding mesh vertices, triangles,
UVs or texture images. An independent current-body replay confirms 610 pairs.
The Blender review shares the original mesh datablocks. A principal-axis parent
implements nonuniform scale without losing affine shear during object-matrix
decomposition; after scene evaluation, world-coordinate replay errors stay
below 0.2 micrometres. The actual Blender and browser renders were inspected.
Evidence: `external-sphincter-forest-fit-comparison.json`, updated
`external-sphincter-current-body-packing.json`, `external-sphincter-fit-render.json`.

### Exact contact loci and bounded in-plane trial

Current native-source, centre-only translation and fitted-pair contacts are now
reported separately. In the contacted-neighbour scope, native muscles show
125 skin, 20 pararectal-node and 61/61 pubo-analis triangle pairs; neither
iliococcygeus nor gluteus has a native-source intersection. Fitted muscles show
105 skin, 71 nodes, 177/162 pubo-analis, 22/20 iliococcygeus and 26/27 gluteus
pairs. These are surface diagnostics, not accepted or rejected tissue attachments.
The new gluteus/iliococcygeus intersections lie at the posterior upper ring edge.

Six smaller in-plane trials retain lateral scale 1.6 and normal thickness but
use sagittal scales 1.2–1.4. Every trial reintroduced canal intersections, so
none was installed. The current ring remains provisional rather than being
narrowed to improve an isolated packing count.

The local review exposes an optional 'Show exact contact curves' checkbox.
Curves are verified noncoplanar triangle-intersection segments with coplanar
pairs explicit; their endpoints are independently checked against both source
triangle planes and barycentric bounds. Curves load only on first request,
reuse their GPU buffers and obey the current section clipping plane. No camera
fitting or anatomical geometry changes occur when enabling them. TypeScript,
camera-intent tests and the actual checkbox/camera comparison pass. Evidence:
`external-sphincter-contact-locations.json`,
`external-sphincter-planar-proportions.json`; optional curves are local review
assets only. Primary MRI abstract references are recorded in the assembly
inventory; complete external muscle subdivisions remain unrepresented.

### Local posterior edge correction (provisional review geometry)

Whole-ring shrinking failed because it narrowed the canal-facing aperture.
A separate C1 posterior compression instead protects all X/Y coordinates and
all source vertices at/above Z = −66 mm. Its Z derivative remains between 0.5
and 1; this describes the continuous map only. The triangulated source halves
also retain their original topology/components, positive normal orientation,
zero self-intersections and zero pairwise crossings. All canal intersections
remain zero. The selected fit changes 26 vertices per half by at most 4.540 mm;
77 per half remain byte-exact. Normal muscle thickness and shape are locally
changed and remain anatomically unaccepted. This is explicitly a mesh correction,
not another rigid placement or a new exact TA98 subdivision.

The same encompassing dependency-group screen drops 610 → 590 triangle pairs
and four → two newly contacted neighbours: both gluteus contacts disappear;
iliococcygeus remains 30/29 pairs. Skin is 120, pararectal nodes 71 and pubo-analis
178/162 pairs. These contacts still require anatomical interpretation; the main
male and female models are unchanged. The local corrected review receives the
refinement with the same 206 vertices and 404 triangles, no new textures/UVs.

Current-contact curves are regenerated and bound to the exact fitted pair SHA;
the viewer rejects curves from another fit. Blender uses two small local mesh
copies for this nonlinear correction, preserving prior source objects. Evaluated
world coordinates replay below 0.03 micrometres error. Actual Blender/browser
views are inspected. Evidence: `external-sphincter-posterior-edge.json`, updated
`external-sphincter-current-body-packing.json`, `external-sphincter-fit-render.json`.
Linear staging/selection and the posterior applier now guard against silently
undoing or applying this correction twice. The existing export pipeline retains
this source-bound corrected shape.

### Upper posterior edge correction

Measured iliococcygeus contact loci occupy roughly Y 882–884 mm, Z −70 to −68 mm.
Four bounded local Y compressions were tested on the current posterior-corrected
ring, protecting all X/Z, anterior vertices (Z >= −65 mm) and lower vertices.
The selected correction anchors Y at 879.5 mm with a posterior blend and positive
Y derivative >= 0.5. It changes 21 vertices per half by at most 1.998 mm; 82 per
half remain exact. Source topology/components and face orientation are preserved;
self, half-to-half and canal intersections remain zero. This is a local shape
change, not a claim of anatomically accepted muscle dimensions or subdivisions.

In the full common dependency-neighbourhood comparison, pairs fall 590 → 510
and newly contacted neighbours fall two → zero. The independent current-body
replay confirms zero gluteus/iliococcygeus intersections and no newly contacted
part among 57 active parts in 45 complete owning groups. Remaining pairs are
skin 105, pararectal nodes 66, pubo-analis 177/162. The unchanged native source
muscles already contact those structures, but changed intersection counts do not
certify correct attachments or tissue volume. Those interfaces, clinical shape,
whole sphincter subdivisions and rectal/anal/proximal digestive continuity still
prevent acceptance or installation in the main models.

The existing corrected local review now contains this upper-edge refinement;
geometry remains 206 vertices/404 triangles with no new UVs/images. Contact curves
are regenerated and SHA-bound. Actual Blender world coordinates replay below
0.03 micrometres, and the actual Blender/browser render is inspected. The upper
trial reads its archived input after application and the applier rejects applying
it twice. Evidence: `external-sphincter-upper-edge.json`, updated packing, contact
and Blender-render reports under `reports/ta98-male/`. Main catalogues unchanged.
