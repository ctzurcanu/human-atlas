# Structural movement dependencies

Geometry approval and placement approval are separate. An approved organ shape
can be repositioned as a complete assembly, preserving all relative part
coordinates, normals and thickness fields. The approved male stomach exterior
is protected by this policy. Movement plans are review data; they do not move
the camera, alter meshes or change the viewer's explosion layout.

## Authoritative current forests (2026-10-07)

The earlier gastric-only forests below are historical. Geometry planning now uses
`scripts/build-spatial-dependency-forest.py` and
`scripts/plan-structural-movement.py`. The existing
`build-structural-dependencies.py` CLI delegates to this complete builder.

The current male catalogue has **6,578** leaves in **5,106** source-object or
explicit assembly trees under **7** source roots. The female has **2,908** leaves
in **2,747** trees under **6** roots. Every mesh occurs exactly once, including
suppressed source alternatives. The original local reference body catalogues are separate
reference forests: **1,949** male and **1,950** female leaves. Identical part IDs in
different catalogues never mean identical current geometry.

Generated readable trees and machine-readable forests are in
`reports/structural-dependencies/current/`:

* `male-forest.txt` / `male-forest.json`
* `female-forest.txt` / `female-forest.json`
* `native-male-forest.txt` / `native-male-forest.json`
* `native-female-forest.txt` / `native-female-forest.json`

Source roots distinguish local reference body male/female, open anatomy, BodyParts3D,
registered/derived gastric and caecocolic assemblies, upper-limb GLBs and the
separate HRA/Visible Human donor sources. A source root is an explicit calibration
scope, not permission to pull the whole source body when positioning one organ.
All active children of an explicitly selected source root preserve the common
transform when that source is scaled or repositioned. Individual source objects
keep their material fragments and partitions together. Explicit stomach tissue
membership closes over all 33 active cells; the shared caecocolic assembly closes
over all 60 active cells and retains its four suppressed pre-partition layers as
reference leaves. Gastric connectors remain separate trees with shared-interface
refitting obligations; their remote mouths cannot silently follow the stomach.
Anatomical aliases and section cap groupings never establish rigid ownership.

The spatial indexes use **actual indexed vertex coordinates**, including geometry
origin and source offset. Unused points in shared buffers cannot inflate a part's
bounds. Exhaustive closed bounding-box broad-phase graphs contain **351,665**
male and **108,256** female candidate pairs, including touching and contained
bounds. Spatial overlap is a conservative review edge, not proof of attachment,
wall penetration or rigid mechanical coupling. Suppressed references are indexed
but excluded from active movement and packing obligations. Index files total
approximately 1.2 MB compressed; they are offline planning artifacts and add no
viewer meshes, textures, draw calls or startup work.

The exact open male polar/pelvis review is also bound in
`male-polar-scene-forest.json`: seven wall intervals form one native-derived
movement tree; six fixed current-pelvis guides retain their complete owning atlas
trees; the old BodyParts3D rectum is one hidden comparison leaf. Thus fourteen
anatomical geometry instances are explicit, thirteen visible by default. Both
current and original source catalogue hashes, wall/context binaries and manifests
are bound. This prevents a review candidate from being mistaken for installed
current anatomy or an old comparison from being counted as a second active organ.

The current forests retain **16,230 male** and **4,540 female** existing typed
anatomical relationship constraints, in addition to source alternatives and
shared tissue junctions. They remain refit/review edges rather than rigid motion
edges. Neighbor review expands to the complete active owning assemblies, not
only whichever fragment initially intersects the proposed movement. Independent
row-wise all-pairs replay agrees with all four saved spatial indexes; source,
chunk and index hashes plus unique leaf/ownership coverage are verified in
`validation.json`.

Before a current-model movement proposal, run the planner with exact part IDs:

```sh
python3 scripts/build-structural-dependencies.py
# Use a Python environment with NumPy for the planner:
python3 scripts/plan-structural-movement.py --sex male --select FJ2571
python3 scripts/plan-structural-movement.py --sex female --select LOCAL:female:large_intestine_rectum:80
# Original source calibration is explicit and separately scoped:
python3 scripts/plan-structural-movement.py --sex male --catalogue native --whole-source-frame frame:local-reference-body-male
```

The planner hashes the current catalogue, all indexed geometry chunks and the
index before returning a plan; stale evidence rejects movement planning. It
expands mandatory ownership groups and lists spatial neighbor groups, source
alternatives, shared tissue interfaces and existing typed anatomical constraints
separately. Mixed-source neighborhoods are explicitly flagged with the source
frame of every group. The current male rectum belongs to BodyParts3D, while the
sacrum and bladder/prostate sources belong to open anatomy; the proposed native
intestinal assembly is a different local reference body source, not the same rectum.
This distinction must remain visible throughout replacement and packing work.

Twelve tests verify actual complete catalogue coverage, mandatory gastric/colic
membership, model separation, suppressed alternatives, explicit source scaling,
nonrigid aliases/cap groups, shared interfaces, exhaustive broad phase including
containment/touching and exclusion of unused shared-buffer points. Whole-body
mechanical attachment completeness and clinical fit remain unverified; unknown
attachments must not be invented from geometry proximity. No geometry, hierarchy
visibility, camera, ignore rules or public deployment are modified by building
or querying these forests.

## Forest-guided pelvic placement trial (2026-10-07)

`screen-forest-colon-placement.py` now queries both verified ownership forests
before comparing the connected seven-part native male colon with complete pelvic
guide trees. All original source arrays independently replay. The same twelve
source self-contact face pairs remain; no new pairs are introduced. Carrying the
descending colon together with sigmoid and rectum removes the 54 interface
crossings of the earlier independently moved distal trial. The six exact source
joins are preserved. This is an outer carrier trial, not accepted whole-colon
wall construction or installed anatomy.

The first common translation of all three current bladder parts and the complete
prostate clears their 429 colon contacts but introduces 38 contacted neighbors.
It is rejected for installation. A second diagnostic common similarity uses the
median six native/current extent ratios after the existing coarse skeletal fit:
uniform scale 0.6792909472261923, centered on the two current organ centroids,
then common translation (-0.27585, 10.25532, 26.06177) mm. This preserves complete
ownership, source normals and relative anatomy. It is not a clinical size target;
the coarse skeletal fit retains its high held-out residual.

The second trial reviews 207 active neighbor pieces, expanding every group
reached by old/new indexed bounds and typed relationships. Aggregate source
surface-contact counts decrease 2,775 to 675, with zero newly contacted neighbors.
Two existing contacts increase: left inferior rectal vein and prevesical nodes.
The four moved parts remain clear of the candidate colon. Surface counts do not
prove tissue-volume disjointness, and loss of contact at ureters/urethra is an
unresolved interface obligation, not a successful route repair. Consequently
neither trial is installed in either main TA98 model.

Both eight-guide saved contexts independently replay affine positions, original
indices and packed normals, with complete ownership closure. Existing Blender
renders compare original/trial complete groups in lateral, sagittal and whole
colon views, preserving prior scenes without saving the main file. The local
review is `/caecocolic-wall-review.html?assembly=male-native-distal&course=polar&context=forest-sized-pelvis`.
Section controls retain the exact camera. TypeScript, camera-intent,
explosion-plane and explosion-viewport checks pass. Source masters, ignore rules
and public deployment remain unchanged.

Evidence: `reports/ta98-male/forest-colon-placement.json`,
`forest-pelvic-placement-sized-trial.json`, `forest-pelvic-trial-sized-packing.json`,
`forest-pelvic-context-replay.json`, and `forest-pelvic-placement-sized-blender.json`.
Next placement obligations are the complete adjoining urinary/reproductive route
trees, the two increased contacts, source colonic defects and coherent female
pelvic fitting. Rigid proximity edges must not drag remote route endpoints.

## Current placement direction: fixed pelvis, move rectum

The user's subsequent instruction supersedes moving or resizing the bladder and
prostate. Their diagnostic trials above remain uninstalled historical evidence.
The skeleton and current pelvic organs stay fixed; correct the rectal course,
keeping all its finite wall intervals together and refitting/carrying the sigmoid
and adjoining colon without opening shared tissue interfaces.

`stage-forest-rectum-posterior-shear.py` binds the current ownership/scene forests
and applies the same posterior exterior field to all seven connected native
colonic source regions. The first construction moves controls up to 17.9335 mm,
retains the physical distal rim and exact distal/full-colon coordinate equality,
keeps sacral crossings at zero and retains the same twelve whole-colon source
self-contact pairs. It adds no triangles or textures. However, the all-x
silhouette constraint creates a visible angular waist in the upper course, and
all attempted inner-field transports retain invalid finite cells. These trials
are rejected for installation and anatomical acceptance. Actual fixed-pelvis
lateral/sagittal Blender renders record the visible defect in
`reports/ta98-male/fixed-pelvis-rectal-move-blender.json`. The next fit must use
actual medial sections and a smooth complete rectal/sigmoid course; the projected
envelope must not be treated as an anatomical centerline or acceptance target.

Whole-tree urinary source evidence is independently recorded in
`reports/forest-urinary-interfaces.json`: current male bladder has a closed outer
carrier and no established urinary holes; the complete male ureters/urethra have
two physical rims each, but their current bladder connections cross without exact
shared controls. Female native bladder/urethra/ureter joins have exact source
controls, while the right ureter retains nine source self-crossings. These are
explicit later repair obligations, not a reason to move those organs during the
current rectal correction.

## Smooth fixed-pelvis rectal continuation

The rejected all-x envelope is replaced by actual x=0 triangle sections in
`measure-fixed-pelvis-medial-course.py`. There are 196 paired height samples;
midline y-parallel separation is 5.240–58.392 mm (median 23.217 mm), a diagnostic,
not clinical perpendicular fascia thickness. Two smooth posterior course lobes
with quintic endpoint/lateral fades give a maximum 20.0024 mm movement. All pelvic
guides stay fixed. The same field carries every adjoining source colon region;
all eight distal interfaces keep shared source ownership and the physical distal
rim stays exact.

One coarse triangle initially fails the finite-cell gate. Two conforming shared
edge refinements add only six control vertices and twelve distal triangles, with
no images. Source face ancestry, areas, winding, tissue ownership, midpoint
coordinates, complete-colon matching and the unchanged twelve original colon
self-contact ancestors independently replay. All eight sheet topologies/self
checks, twenty-eight distinct-sheet checks and seven closed finite shells pass.
The actual smooth lateral/sagittal Blender renders remove the earlier sharp
silhouette waist but still show unresolved pelvic contacts. No anatomical
acceptance is inferred from the wall checks.

The verified current forest selects 411 active pieces across 334 complete owning
neighbor groups from old/new indexed bounds. Distal outer/luminal contact totals
decrease 13,443 to 7,801; five newly contacted neighbors remain, including anterior
sacrococcygeal ligaments and external iliac nodes. Prostate contacts slightly
increase (263 to 267). This placement remains uninstalled. The full upstream
colon packing, containment, clinical course, original source crossings and
adjoining vessel/fascial interfaces remain unresolved. The lower count does not
constitute acceptance.

Local review: `/caecocolic-wall-review.html?assembly=male-native-distal&course=smooth&context=pelvis`.
It includes the complete fixed bladder and sacral trees, plus a disabled-by-default
old rectum comparison; there is no second active rectum. The shared wall payload
is 235,752 bytes. Sectioning preserves the exact camera. TypeScript, two conforming
subdivision tests and all camera/explosion tests pass. Evidence is bound in
`reports/ta98-male/fixed-pelvis-medial-course.json`,
`fixed-pelvis-smooth-rectal-course-adaptive.json`,
`fixed-pelvis-smooth-rectal-course-replay.json`,
`fixed-pelvis-smooth-rectal-packing.json`, and
`fixed-pelvis-rectal-move-smooth-blender.json`. Both main models and source masters
remain unchanged; female fitting is pending.

## Localized rectal course evidence

The completed smooth-course packing screen is followed by whole contact-pair
localization in `localize-fixed-pelvis-rectal-contacts.py`. Prostate contact affects
both the low rectum and the unclassified distal source. Bladder contacts affect
the rectum; newly contacted sacrococcygeal ligaments lie at 958–965 mm height.
These distinctions prevent treating one total as a single anatomical defect.

A third posterior lobe at the low rectal course and a narrower upper lobe yield
the local `--adaptive --localized` trial. All skeletal, bladder and prostate
guides remain fixed. Rectum/prostate contacts decrease 126 to zero; the distal
source retains 53 prostate pairs. Bladder fundus contacts decrease 264 to 199,
neck contacts 26 to zero, median sacral artery contacts 120 to 32, and both newly
contacted sacrococcygeal ligaments decrease 48 to zero. These are outer-plus-inner
triangle counts and do not prove tissue-volume clearance.

Source/refinement/ownership/shared-interface replay and all eight sheets,
twenty-eight distinct-sheet and seven closed-shell checks pass. The local
candidate adds 27 shared vertices and 54 triangles, no images, and retains the
same twelve original colon crossing ancestors. Actual existing-Blender review
nevertheless rejects the visible sharp return toward the fixed distal mouth and
an angular upper shoulder. No new viewer export or main-model installation is
made. The source mouth remaining fixed was a construction assumption, not the
user's instruction; future fitting must define its pelvic-floor/distal-course
landmarks before moving the complete continuous distal assembly. All fixed
pelvic guides remain protected.

Evidence: `reports/ta98-male/fixed-pelvis-rectal-contact-localization.json`,
`fixed-pelvis-rectal-contact-localization-localized.json`,
`fixed-pelvis-smooth-rectal-course-localized-adaptive.json`,
`fixed-pelvis-smooth-rectal-course-localized-replay.json`, and
`fixed-pelvis-rectal-move-smooth-localized-blender.json`.

## Distal outlet and pelvic-floor continuation

`audit-distal-pelvic-floor-landmarks.py` binds both native and current ownership
forests and records exact welded distal boundaries, nine real transverse sections
per distal source, and complete selected pelvic-floor trees. Male has three
native and nine current floor pieces; female has three in each known native
frame. The current male external anal sphincter bounds are centered near
Z=-60.534 mm, while the current polar source distal rim is near Z=-20.709 mm.
This 39.825 mm source-family discrepancy is a geometric observation, not a
certified aperture landmark or clinical position. The female source differs and
must not inherit a male translation.

The `--adaptive --floor-outlet` plan therefore defines one continuous posterior
field for the whole source colon and all distal wall intervals, allowing the
physical distal rim to move with the low rectum instead of forcing it to return
to its original point. Skeleton, bladder, prostate and floor guides remain fixed.
The plan is recorded before geometry mutation. Shared tissue coordinates remain
exact; the original mouth is transported, not plugged or represented by a new
solid. The 24 shared additional controls and 48 triangles introduce no textures.
All source face/refinement ancestry, fields, ownership and matching full-colon
coordinates independently replay; eight sheets, twenty-eight sheet pairs and
seven closed finite tissue shells pass. The twelve original whole-colon source
crossing ancestors remain unchanged.

Focused outer/inner screening now detects zero prostate contacts in both rectal
and distal regions, zero bladder-neck contacts, and zero contacts at both anterior
sacrococcygeal ligaments. Bladder fundus retains 237 pairs and median sacral artery
36. Complete current-body packing, sphincter/tissue interfaces and containment
remain unverified. Existing Blender sagittal review confirms removal of the sharp
distal return bend but still shows an angular upper rectal/sigmoid shoulder. This
trial is uninstalled and not anatomically accepted. It is progress toward a
continuous course, not completed rectal anatomy.

Evidence: `reports/distal-pelvic-floor-landmarks.json`,
`reports/ta98-male/fixed-floor-outlet-course-plan.json`,
`fixed-pelvis-smooth-rectal-course-floor-outlet-adaptive.json`,
`fixed-pelvis-smooth-rectal-course-floor-outlet-replay.json`,
`fixed-pelvis-rectal-contact-localization-floor-outlet.json`, and
`fixed-pelvis-rectal-move-floor-outlet-blender.json`. No main-model geometry,
camera, hierarchy, ignore rule or public deployment changes.

## Complete outlet-trial packing and exact-section follow-up

The full `--floor-outlet` screen reviews 411 active neighbor parts in 334 owning
trees. Contacts decrease 13,443 to 8,187, but seven previously clear neighbors
are now contacted, including Skin, both external anal sphincters, both
iliococcygeus pieces, intermediate external iliac nodes and pararectal nodes.
The candidate therefore fails placement acceptance despite clearing prostate and
seminal glands. Source-family bounding-box centers are explicitly insufficient
outlet landmarks. Contact localization records skin pairs in both the rectum
and distal source, while sphincter pairs lie in the distal source only. No main
installation is performed.

The upper angular shoulder is also traced to coarse sampling of the smooth
spatial field: an original source edge has 5.328 mm midpoint/chord error. A local
`--curve-fine` trial refines only edges with more than 0.4 mm interpolation error,
along with failed finite-cell facets, under an explicit 1,024-control budget. It
adds 260 shared controls and 520 triangles, no images. Initial finite-cell,
exterior self/sacral and unchanged whole-colon crossing counts pass; complete
independent wall replay and actual render remain pending. This does not repair
the failed outlet packing and is not exported or installed.

`exact-plane-segments.py` now handles on-plane vertices and shared plane edges,
welds duplicate segments and reports coplanar faces separately. Four fixtures
cover vertex/opposite-edge intersections, welded shared edges, coplanar surfaces
and tangent points. Independent regeneration to `fixed-pelvis-medial-course-exact.json`
preserves the original 237 sacral and 115 rectal section segments, all 196 samples
and prior gap statistics; neither source has coplanar faces at this midline.
Historical evidence remains untouched. The precision correction matters for
subsequent skin/sphincter aperture contours rather than changing this source's
existing midline result.

Evidence: `fixed-pelvis-smooth-rectal-packing-floor-outlet.json`,
`fixed-pelvis-rectal-contact-localization-floor-outlet.json`,
`fixed-pelvis-smooth-rectal-course-floor-outlet-fine-adaptive.json`, and
`fixed-pelvis-medial-course-exact.json` under `reports/ta98-male/`.

## Actual distal contour constraints

`audit-distal-passage-contours.py` binds the fixed male Skin (BodyParts3D), both
current external anal sphincter halves (open anatomy), and the complete source-bound
candidate exterior/lumen. Exact transverse planes at 857, 865, 875 and 885 mm
produce complete welded contour graphs, with coincident vertices/edges handled
and coplanar triangles separately reported. Skin has fourteen closed cycles at
each height, including paired shell/body and limb contours. The sphincter halves
have one cycle each at the first three planes and no section at 885 mm. The
candidate outer and inner contours each have one cycle per plane. These are
section facts, not accepted aperture landmarks or a proof of tube patency.

Actual plotted contours confirm candidate/muscle intersection at 865 mm and
candidate crossings of both perineal skin sheets at 875 mm. A group bounding
center or whole-body bounding box cannot constrain these local interfaces. The
next continuous distal fitting must include cross-sectional shape as well as
position, preserving shared tissue geometry and full adjoining ownership rather
than repeating a translation-only trial. Fixed skeleton/bladder/prostate remain
protected. There is no geometry installation or additional accepted TA98 term.

Evidence: `reports/ta98-male/distal-passage-contours.json`, and reproducible SVG/PNG
diagnostics in `outputs/distal-passage-contours.*`. Four exact-plane fixtures pass.
The image is a direct plot of actual source triangles, not an anatomical diagram
or copied reference image; whole uncropped segments remain in the JSON.

## Distal section feasibility and source-frame conflict

`fit-distal-passage-sections.py` reconstructs ordered closed loops from actual
section graphs, checks concave polygon containment, proper segment crossings and
sampled boundary separation. A bounded scale/translation family initially finds
free space by moving the tube anterior to the sphincter; that bypass is rejected
and its report preserved as `distal-passage-section-free-space-only.json`.
The current fit plan additionally requires the tube center to remain within the
combined sphincter contour envelope. This convex envelope is a route-location
constraint, not a certified aperture or a replacement tissue surface.

Under these explicit constraints, there are zero feasible section trials at
857/865/875 mm, and 1,470 at 885 mm where the sphincter has no section. This is
only a failure of the tested centered uniform-section-scale/Z-translation family,
not proof that every possible reconstruction fails. Of the fixed sphincter
contour points, 14/14, 44/50 and 44/44 lie outside the selected current paired
trunk skin interior at those first three heights. These source-frame facts must
be resolved before constructing a purported passage; zero contacts obtained by
bypassing the sphincter are not anatomical progress.

`compare-perineal-skin-source-contours.py` independently loads the complete two
fragment native skin ownership tree. The original source skin remains in its
recorded coordinates, with no registration or substitution. At the three low
planes, 28/28, 44/100 and 34/88 duplicated section-endpoint samples from the
current sphincter lie outside all native contour interiors too. Source families
and sampling definitions are explicit; these counts do not certify clinical
registration or justify moving the fixed skeleton. A coherent perineal frame
and corresponding openings must precede an accepted continuous distal fit.

Five section-feasibility fixtures verify unordered-edge reconstruction, rejection
of noncycle constraints, concave containment rather than bounding-box inclusion,
crossings whose endpoints lie outside the other polygon, and route-envelope
construction. No mesh, runtime tissue/coverage, hierarchy, camera, ignore or
deployment changes occur. Evidence: `reports/ta98-male/distal-passage-section-fit-plan.json`
and `perineal-skin-source-comparison.json`. All failed candidates remain local
and uninstalled; the goal remains incomplete for both sexes.

## Same-source regional coverings and height lineage

The current male catalogue additionally contains open anatomy anal, urogenital and
gluteal regional surfaces. `audit-same-source-perineal-covering.py` expands all
six selected owning trees to fourteen active pieces rather than overlooking
their fragments. They are open regional coverings, not automatically a whole
skin wall or physical anal opening. Each anal half has fourteen source triangles
and one ten-edge boundary. Combined selected covering surfaces retain 164 open
edges across five boundary cycles; boundary roles remain unverified.

Anal-region patch bounds reach 791.602–836.950 mm height, while current external
anal sphincter bounds are 856.990–880.113 mm. This distinction must not be hidden
by treating a region name or common provenance as an established outlet. The
source regional skin and separate BodyParts3D paired whole Skin envelope are
actual different current geometries. Existing-Blender lateral/sagittal renders
inspect both alongside the same fixed sphincter and prostate guides; no automatic
source substitution or copied patch is accepted.

An independent current/master/pre-height replay rules out missed/double scaling
for eight selected fixed guides: positions and faces are exact against the
current protected master, faces exact against the pre-height source, and all
coordinates obey the recorded common 183 cm similarity within 0.0001 mm. This
does not certify source anatomy or tissue packing, but prevents another false
height correction. The mismatch precedes the common transform.

Evidence: `reports/ta98-male/same-source-perineal-covering.json`,
`perineal-height-lineage-replay.json`, `perineal-source-coverings-blender.json`,
and `outputs/perineal-source-coverings-{lateral,sagittal}.png`. Prior Blender
scene/objects are preserved without saving the main file. No atlas geometry,
hierarchy visibility, camera, ignore rule or deployment changes. The next
anatomical plan must establish a single referenced perineal covering/outlet and
fit the continuous distal tube to it, keeping skeletal/internal pelvic guides
fixed and source identities explicit.

## Corrected distal fitting scope and authority branch

The earlier decision to hold every pelvic muscle fixed was a construction
assumption, not a user invariant. The skeleton, bladder and prostate remain
fixed, but the continuous distal passage and surrounding external anal sphincter
require coordinated anatomical fitting. This does not establish a new rigid
attachment by proximity. Complete source ownership remains mandatory, and
puborectalis/perineal remote attachments require explicit refitting instead of
blindly following the outlet.

The authority extract puts external anal sphincter A04.5.04.012 beneath
anal-triangle muscle A09.5.00.003, whose parent is A04.5.05.001; this branch is not
descended from the perineum root chosen in the earlier pelvic inventory. The
new `TA98-DISTAL-SPHINCTER-ASSEMBLY.md` and
`reports/distal-sphincter-assembly-inventory.json` explicitly include all 38
authority descendants of rectum, anal canal, anal-triangle muscle, puborectalis,
perineal body and subcutaneous perineal tissue, with exact parents, current
source parts/groups and missing representations for both models. Independent
recursive SQLite-CTE replay verifies exact unique coverage and literal parent
edges in `distal-sphincter-assembly-inventory-replay.json`. Historical SHA-bound
inventories remain unchanged.

The MRI component study
https://pubmed.ncbi.nlm.nih.gov/16319250/ distinguishes sphincter components using
visible separation and differing origins/insertions in nulliparous women. It is
structural evidence, not a male size template, a source geometry donor, or proof
that every floor muscle should share a rigid transform. Missing sphincter
subdivisions remain missing even when the whole source is present. The current
native distal source remains unclassified until its actual landmarks establish
anal-canal correspondence. This corrected anatomical scope precedes subsequent
geometry construction; no model, camera, hierarchy or deployment changes.

## Historical gastric-only trees


Run `python3 scripts/build-structural-dependencies.py`. The local outputs are
`reports/structural-dependencies/male.json` and `female.json`, bound to the
current runtime SHA256. Every existing mesh appears exactly once, including
suppressed meshes. A suppressed mesh's presence does not enable it.
Each JSON contains a named parent/child `forest`. Readable exports are
`reports/structural-dependencies/male-tree.txt` and `female-tree.txt`.

Current inventory: male **6,301** leaves in **67** groups/candidates and 5,932
standalone leaves; female **2,808** leaves in **3** groups/candidates and 2,767
standalone leaves. Eight regression tests verify exact mesh coverage, named tree
leaves, gastric movement closure, junction/section separation, attachment review
from either endpoint, and coordinated movement without dragging review neighbors.

Each forest has three kinds of information:

* **Whole-organ movement groups:** explicit stomach regional membership;
  selecting one member in the offline planner includes all four regions.
* **Candidate groups:** existing `sectionAssembly` membership, including bone
  and gland partitions. These require attachment review before coupled motion;
  the planner does not automatically expand them.
* **Standalone leaves:** attachments not yet inventoried. They remain explicit
  unknowns. Catalogue concepts are separately retained as reference sets,
  including overlapping aliases and donor collections; they are not rigid edges.

The installed stomach movement group in each model contains all **33** finite
tissue cells recorded by `gastricCompleteWall.partIds`, across the four regions:

```text
Stomach assembly — preserve relative pose
├── Cardia
├── Fundus
├── Body
└── Pyloric part
```

The installed wall tissues share that same transform. The hierarchy below also
includes microscopic features whose distinct representations remain pending:

```text
Stomach wall — installed tissue cells; microscopic coverage incomplete
├── Tunica serosa
├── Tela subserosa
├── Tunica muscularis
│   ├── Stratum longitudinale
│   ├── Stratum circulare
│   │   └── M. sphincter pyloricus — local thickening, shared tissue
│   └── Fibrae obliquae
├── Tela submucosa
└── Tunica mucosa
    ├── Lamina muscularis mucosae
    ├── Plicae gastricae
    ├── Areae gastricae
    ├── Plicae villosae
    ├── Foveolae gastricae
    └── Glandulae gastricae
```

This movement tree does not replace the anatomical inventory or authoritative
TA98 ancestry in `TA98-STOMACH-PARTS.md`. Microscopic detail and boundary/opening
representations do not become independent solids just because they occur here.

## Cross-organ constraints

These connections form a graph outside the part trees. They never silently
propagate a rigid transform into another organ:

| Connection | Required review after movement |
| --- | --- |
| Stomach ↔ oesophagus | Cardial paired wall rim, lumen continuity, local nerve course |
| Stomach ↔ duodenum | Pyloric wall/rim, sphincter, outlet lumen |
| Oesophagus ↔ diaphragm | Hiatus clearance and crural attachment; skeletal constraints |
| Stomach ↔ liver, pancreas, spleen | Actual surface/section packing and attachment endpoints |
| Stomach ↔ omenta, ligaments, vessels, lymphatics, vagal branches | Inventory individual local attachments; refit affected courses |

These are declared review requirements, **not verified mesh junctions**. Gastric
arteries or a vagus nerve do not cause their entire upstream tree to move.
The male `TA98:gastric-tract-wall` cap grouping spans multiple organs; this is
specifically excluded from automatic rigid movement. Female duodenal donor
concept unions do not establish a common source frame with native anatomy.

Installed male cardiac/pyloric connectors explicitly require endpoint refitting
when the stomach moves. They appear in dependency review, outside the rigid
33-cell closure: their remote mouths belong to adjoining conduits. Selecting a
connector also raises whole-stomach review without automatically moving it.
The builder verifies installed wall membership against its assembly record and
rejects a missing tissue cell. Older catalogues without that record retain the
explicit four-region fallback.

## Planning and evidence

An example plan is produced with:

```sh
python3 scripts/build-structural-dependencies.py --sex female --select GASTRIC:female:019:012
python3 tests/structural-dependencies.test.py
```

The planner returns the exact proposed part set and all touched junction,
packing and anchor constraints. It reports unknown attachments and candidate
groups requiring review. Selecting a gastric vessel or vagus candidate also
flags the stomach attachment review, without adding the stomach to its movement
set. `dependencyReviewPartIds` lists affected neighbors separately from
`partIds`; `coupledMoveGroups` lists the actual movement closures. Repeat
`--select` to plan explicitly coordinated selections. It applies no transforms. Before repositioning any
group, compare the supplied local reference geometry and sections, specify
the displacement and attachment refits, replay shared seams/normals, and run
the manual's geometry/packing acceptance checks. Source-relative reference
measurements do not authorize copying proprietary meshes or textures.

`scripts/review-dependency-reference-placement.py` reads the supplied male and
female 2024 GLBs with exact source-hash checks. It has inspected 51 native
reference meshes per sex for stomach, oesophagus, duodenum, liver, pancreas,
spleen and diaphragm. The local `reference-placement.json` retains all 11
alternative candidate bounds per sex separately, with source-relative stomach
offsets. No candidate was silently chosen as a certified whole organ. These
coarse measurements establish context for subsequent actual surface/section
review; they do not establish clearance or a justified displacement by themselves.

Whole-body attachment completeness remains unfinished. The generated forest
inventories every current mesh; it does not claim that every mechanical
dependency has been established or that any anatomical assembly is accepted.

## Male hepatic fit constraints

The current source-bound gastric packing review additionally inventories native hepatic dependants in [TA98-HEPATIC-FIT-DEPENDENCIES.md](TA98-HEPATIC-FIT-DEPENDENCIES.md). The organ and ligament material regions share 100 exact source coordinates; separate hepatic territory surfaces have additional source coincidences. These constraints must survive a coupled local fitting field. Neither proximity nor shared source coordinates certifies an anatomical joint or warrants rigidly moving the entire inflow/outflow tree. Cava parts and lesser omentum remain cross-organ review constraints. No motion group has been promoted to anatomical acceptance.

## Independent mirrored source instances

The current male catalogue retains 215 publisher objects mirrored into a separate left-side instance. Shared `sourceId` describes the publisher object, not a rigid link between original and mirror. Ownership now keys the explicit mirrored instance separately while retaining every material fragment of each instance in one tree. Explicit whole-frame calibration still includes both instances.

The median-root repair scope exposed the error: five left-side parts previously closed over nine parts, pulling right roots and cords along. The regenerated male forest has 6,577 leaves and 5,320 ownership groups; either five-part repair scope now closes over exactly five parts. Female and both native forest owner assignments replay unchanged under the new grouping rule. No geometry, camera or runtime resources changed. Mechanical attachment completeness remains unverified.

Evidence: `server/forest-source-instance.test.py`, `reports/structural-dependencies/current/validation.json`, and the source-bound `reports/ta98-male/median-root-assembly-repair-plan.json`. The 61-term authority inventory and planned cord–root–trunk interfaces are in `TA98-BRACHIAL-PLEXUS-REPAIR.md`; planned anatomical interfaces do not establish accepted mesh joints.
