# Female external anal sphincter: construction inventory

## Separate region inspection (2026-10-07)

The vessel-constrained local review now exposes pars subcutanea (.013), pars
superficialis (.014) and pars profunda (.015) individually, beneath the whole
external anal sphincter (.012). Whole and subdivision presentations are mutually
exclusive. All regions use the same paired control fields and existing ownership;
this adds no anatomical geometry to the atlas.

The review export uses a common diagonal at subdivision interfaces. Exact
coordinate comparison verified 292 paired interface triangles. A shared cap is
omitted while both adjoining regions are visible and restored when its neighbour
is hidden. Smooth field normals are shared across parts; exposed boundary faces
have separate flat normals. This removes the observed alternating triangular
seam artifacts without changing the candidate coordinates or fitting evidence.

Actual lateral section inspected and saved in
`outputs/female-sphincter-three-regions-section.png`. The construction still has a
provisional upper attachment and region boundaries; its native canal source does
not establish complete inner/outer wall anatomy. It remains uninstalled and is
not accepted as completed TA98 coverage. TypeScript check passed. Related native
anal canal, rectum, puborectalis and superior rectal vein remain independently
visible for interface review. Full goal remains incomplete.

## Smoother vessel-bounded rim and Blender inspection (2026-10-07)

A separate `source-guided-smooth-vessel-rim-trial.npz` reduces steep changes in
the additional upper-edge displacement. The limit is 0.35 metres displacement
per metre of carrier-edge length. This is a construction smoothness constraint,
not a biological measurement. It only decreases the existing nonnegative lift,
so all 1,531 native-vein height constraints remain valid. Thirty-seven controls
move; maximum additional lift is 1.554930 mm. Lower parts, shared ownership,
X/Z coordinates and native skin/canal/rectum/puborectalis meshes remain exact.

Independent replay verifies the slope bound, field reconstruction, closed
region shells, additive volumes, zero invalid cells, zero detected self-crossings,
zero superior-rectal-vein crossings and sampled canal exclusion from muscle.
Current-forest screening introduces no new visible contact owners. Embedded
neurovascular interfaces remain unresolved. Field/index/ownership storage stays
36,868 bytes, with no new controls, faces, UV arrays or textures.

The already open Blender 5.2.1 session successfully renders both assembled and
sagittal views. Two new diagnostic scenes preserve the prior scene; the main
file is unsaved. Actual renders show improved upper-edge continuity but retained
local angularity and unfinished attachment morphology. These views are evidence,
not anatomical acceptance. Candidate remains uninstalled. Evidence:
`external-sphincter-smooth-vessel-rim-validation.json` and
`external-sphincter-smooth-vessel-rim-render.json`; renders:
`outputs/female-smooth-sphincter-assembled.png` and
`outputs/female-smooth-sphincter-section.png`.

Supplied local guide morphology can now be inspected separately with
`female-sphincter-interface-review.html?candidate=reference`. Display scale is
diagnostic; its coordinates are not transferred into the atlas. The original
guide defects are retained. Primary reference rechecked through the accessible
[publisher article](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0132226):
its sectional study describes U-shaped deeper slings and a circumferential
subcutaneous component. This differs from the divisions in the
[female MRI study](https://pubmed.ncbi.nlm.nih.gov/16319250/).
Neither account justifies blindly snapping the existing native puborectalis
onto the entire constructed rim. Remaining work must reconcile regional
morphology and perineal/anococcygeal supports with the selected TA98 convention.

Planning only. Original meshes and current catalogues are unchanged. Reference rights and clinical scale remain unverified.

| TA98 | Parent | Latin | Structure | Current female representation | Local reference guide |
|---|---|---|---|---|---|
| A04.5.04.012 | A09.5.00.003 | m. sphincter ani externus | external anal sphincter | missing distinct muscle representation | three-piece whole assembly |
| A04.5.04.013 | A04.5.04.012 | pars subcutanea | subcutaneous part | missing distinct muscle representation | 35b3f477e80919c52a14c98569a7cad7 |
| A04.5.04.014 | A04.5.04.012 | pars superficialis | superficial part | missing distinct muscle representation | decb09b10455f2f141977961e0c7ef1e |
| A04.5.04.015 | A04.5.04.012 | pars profunda | deep part | missing distinct muscle representation | e15e12752af11253f13a5299a338832a |

All four terms describe macroscopic muscle tissue/regions, not openings, wall layers or microscopic features. The whole muscle surrounds the canal; the anal opening is a space and must never receive a solid plug.

Adjoining interfaces to review: native anal-canal exterior and patent lumen, rectum, skin at the anal margin, puborectalis/pelvic floor, perineal body, anococcygeal support, inferior anal nerves and local vessels. Functional supply does not establish mechanical attachment.

Inspect the complete native canal and rectal surfaces before reconstruction. Any new whole-muscle representation must have shared tissue interfaces between its three parts, preserve the lumen, avoid duplicated overlaps and use procedural appearance. Reference vertices, triangles, UVs and images must not be copied.

Detailed current-source topology and exact sections: `reports/ta98-female/external-sphincter-construction-inputs.json`. No construction is accepted by this inventory.

## Native interface and reference findings

The native anus is one annular open surface with 62/11-edge rims and no self
crossings. The native rectum has 24/11-edge rims and no self crossings. Their
11-point rims match exactly; the two source surfaces have no mutual proper
crossings. This establishes a shared source exterior, not finished inner/outer
canal walls, histological coverage or biological lumen acceptance.

The supplied female reference has individually named subcutaneous, superficial
and deep guides. The superficial guide has two components and 11 crossing
pairs, and all guides retain open boundaries. It must not be transferred as
finished geometry. Gross dimensionless bounding proportions guide the original
construction; source vertices, faces, normals, UVs and texture images are not
copied.

## Original local construction stage

The new whole-muscle construction uses 512 carrier controls, two shared
surfaces and 960 carrier triangles. Its finite, hollow assembly has 1,024
surface vertices and 2,048 triangles, with an open central passage rather than
a plug. Carrier faces belong disjointly to three provisional regions. The
field/index/ownership arrays occupy 24,768 raw bytes; this excludes normals and
runtime object overhead. No UV arrays or sampled images are used. Runtime
normal packing and procedural circumferential muscle appearance remain future
integration work.

The native transverse opening is nonhorizontal. The guide begins above its
highest rim point, leaving coverage to the actual skin margin explicitly
unfinished. Flat region interfaces derived from broad reference proportions
are construction boundaries, not accepted anatomical subdivision boundaries.

The first stage has eight skin triangle crossings. Exact local support rays
locate a single outer vertex 0.378 mm below the skin floor. A smooth angular
and axial lift, bounded at 0.8 mm, clears these contacts while preserving the
entire canal-facing muscle surface. Unsupported rays through openings remain
explicit; this diagnostic is not a whole-body containment proof.

Independent replay of the stored Float32 candidate confirms zero skin and
native-canal crossings, zero self crossings and zero invalid tissue cells. The
whole assembly and all three separately closed region shells are genus one.
All sampled native canal points in the construction height range lie outside
the muscle material interval, preserving the void. Original native canal and
rectum meshes are unchanged. Actual assembled and sagittal Blender views were
inspected; all pre-existing scenes and the main file are preserved.

Complete-owner screening retains native inferior/middle rectal vessels,
pudendal nerves and perimuscular rectal plexus contacts for embedded-tissue
interface review. Three contacting HRA uterine/ligament donors are independently
registered hidden reference candidates; their contacts are recorded separately
and must not drive blind changes to native anatomy. Neither nerve/vascular
crossing counts nor hidden donor geometry certifies an anatomical attachment.

The candidate is local and **not installed in the main female atlas**. Remaining
work includes the actual distal margin, sloping subdivision interfaces,
neurovascular/pelvic-floor connections, clinical dimensions and full anatomical
acceptance. This stage adds construction evidence, not completed TA98 coverage.

Evidence under `reports/ta98-female/`: `external-sphincter-original-construction.json`,
`external-sphincter-skin-floor.json`, `external-sphincter-skin-clearance.json`,
`external-sphincter-original-validation.json` and
`external-sphincter-original-render.json`. The original and corrected staged
files are both retained under `.local-models/ta98-review/female-external-sphincter-original-construction/`.

## Curved subdivision plan

Replace only the provisional horizontal region partitions with two smooth, sloping interfaces. Estimate broad angular height trends from area-weighted, second-order harmonic summaries of the three labelled local reference guides; transfer only these dimensionless trends, never reference vertices, faces or images. Clip the already repaired carrier triangles against both interface fields using shared edge controls. Outer and canal-facing surfaces retain their existing piecewise planar shape, distal/proximal rims and local skin clearance. Part ownership is disjoint, and adjoining parts share each tissue boundary with opposite orientation. Replay current native canal/rectum coordinates and skin crossings against the current ownership forest before any acceptance. These guide-based partitions remain provisional until anatomical review; this operation does not resolve the unfinished distal margin or pelvic-floor attachments.

## Curved region evidence

The local carrier now has 658 shared controls and 1,252 disjointly owned triangles. Its complete shell has 2,632 triangles. Field/index/ownership storage is 32,068 bytes, excluding offline replay weights, normals and runtime objects; UVs and images remain zero. All three parts and the whole are closed genus-one shells with zero detected self-crossings. Current native canal, rectum and skin replay without crossings. Distal/proximal rims stay exact; the whole-volume relative change is 1.21 × 10⁻⁸ from Float32 interpolation. The 146 shared boundary edges use one canonical rim diagonal with reverse orientation on the adjoining part; part volumes sum to the whole. Whole outer-rim triangulation is preserved.

Actual assembled and sagittal views were inspected in the existing Blender, preserving all prior scenes and the main file. This remains an uninstalled local construction: the anal margin, clinical proportions and pelvic-floor/neurovascular interfaces remain unfinished. Evidence: `external-sphincter-curved-regions.json`, `external-sphincter-curved-validation.json`, `external-sphincter-curved-render.json`; images `outputs/female-curved-sphincter-assembled.png` and `outputs/female-curved-sphincter-sagittal.png`.

## Distal margin extension plan

Retain the verified curved carrier and all original native organs. Add three low-resolution circumferential rows below its current distal rim, owned solely by the subcutaneous part. The new canal-facing lower guide follows broad angular variation of the native 62-edge anal rim with outward clearance. Its outer guide tapers toward the margin and remains above locally measured skin support where available. Missing vertical support is a diagnostic, not permission to infer body containment: independently test the complete skin, canal and rectum surfaces, every new tissue cell and whole/part shell before retaining the stage. Preserve the two existing curved tissue interfaces and all original carrier vertices/triangles. Clinical muscle-to-skin and pelvic-floor attachments remain separate requirements; no opening receives a solid plug.

The first distal extension is rejected: eight invalid tissue cells, three self-crossing pairs and six skin crossing pairs, with no native canal/rectum/puborectalis crossings. Retain its archive and report. The next bounded posterior junction trial may lift paired controls by at most 0.8 mm within the existing low subcutaneous patch, tapering to zero over 4 mm above the former lower rim. Existing curved muscle-part interfaces and upper controls stay exact. A common-height new annular margin uses measured skin support; replay all source organs and exact full-surface checks again. This local move is a construction trial, not accepted anatomy.

The bounded posterior lift clears self/skin crossings but the common-height lower guide intersects the flaring native canal (86 pairs), and one invalid cell remains. It is also rejected. Recompute canal-facing clearance with exact triangle-ray intersections at each guide's actual height, rather than reusing the narrow opening's radius at a higher station. Reference only the native exterior; absent rays remain flagged. This refines local clearance without moving the canal or inventing an inner wall. Recheck full stored geometry before retaining a candidate.

## Retained distal geometry stage

The current lower extension has 754 shared carrier controls and 1,444 carrier triangles; its closed whole shell has 3,016 triangles. Field/index/ownership storage is 36,868 bytes, excluding packed normals and runtime overhead; no UVs or image assets. Independent Float32 replay verifies zero invalid tissue cells, self-crossings, or crossings against current skin, native anus, rectum and puborectalis. Native canal samples remain outside the muscle material. The source canal/rectum/skin are unchanged. Maximum existing axial movement is 0.800014 mm (Float32 replay); existing radial movement is zero. Controls above the 4 mm lower patch and superficial/deep part interfaces stay exact. The new rows use height-dependent native triangle-ray clearance and a bounded 0.35 mm local radial clearance field at the remaining first-cell patch.

All three parts and the whole remain closed genus one; part volumes sum exactly to the whole within the numerical replay tolerance. Full ownership closure screens 49 current parts and records 13 contacting source owners, with no new visible contact owners. Existing inferior/middle rectal vessels, pudendal nerves and perimuscular rectal plexus require entering/embedded-tissue interface review; independent hidden donor organs remain registration candidates.

Actual assembled and sagittal views were inspected in the existing Blender; the sagittal view includes exact native skin context. Existing scenes and the main file are preserved. This remains a local, uninstalled stage: clinical muscle proportions, anterior/perineal and posterior/anococcygeal attachments, puborectalis/deep-part blending and finished native canal wall/lumen identities remain unaccepted. No claim of complete sphincter anatomy or complete TA98 coverage follows from the geometry pass.

Evidence: `external-sphincter-distal-margin.json`, `external-sphincter-distal-margin-validation.json`, `external-sphincter-distal-forest.json`, `external-sphincter-distal-margin-render.json`; `outputs/female-distal-margin-sphincter-sagittal.png`. Rejected first and lift-only stages are retained separately.


## Current deep interface evidence

Replayed the saved candidate against the current female catalogue after appearance installation. Every previously verified native skin/canal/rectum/puborectalis position/index hash remains exact, and full current intersection/material checks pass. This does not accept the assembly. The attachment inventory now resolves concept-to-mesh mappings as well as direct part fields: A04.5.04.007 is already represented by native Puborectalis, rather than missing as the earlier audit mistakenly reported.

Complete vertex-to-face and interior edge-to-edge testing of the deep outer carrier against native Puborectalis establishes a 2.871045 mm minimum separation (previous whole-muscle vertex/centroid sample: 2.886347 mm). Actual assembled and section views in the already open Blender show the gap; all prior scenes survive and the main file remains unsaved. The guide line is diagnostic, not tissue. No native geometry was moved and no new atlas anatomy installed.

Primary article excerpts support conflicting subdivisions: the nulliparous female MRI study (https://pmc.ncbi.nlm.nih.gov/articles/PMC1479222/) distinguishes winged EAS fibers from puborectalis; the serial-section reconstruction (https://pmc.ncbi.nlm.nih.gov/articles/PMC4549266/) treats the deep EAS portion as puborectal tissue and describes anteriorly open slings. Direct article opens currently hit an access challenge; only indexed primary excerpts were consulted this turn. These do not calibrate dimensions. The next construction must resolve the existing puborectalis/deep-part footprint and anterior perineal support within one native assembly. Do not blindly close the measured gap or add another complete muscle ring over the native sling. Preserve disjoint regional ownership and the existing distal/skin clearance, and repeat all current forest/interface tests after any change.

Related structures remain: native anal canal and rectum, puborectalis and pubococcygeus, perineal body, anococcygeal support, inferior anal/rectal nerves through the pudendal route and local vessels. Missing perineal/anococcygeal tissue remains explicit. Evidence: external-sphincter-deep-interface.json, external-sphincter-deep-interface-render.json; outputs/female-deep-sphincter-interface-section.png.

## Bounded deep-edge construction trial

Test a common axial lift of paired controls in the uppermost 4 mm, restricted smoothly to bilateral native puborectalis proximity near the measured minimum. Maximum trial lift 2.5 mm; X/Z stay exact. Keep subcutaneous/superficial parts, distal margin and part interfaces fixed. This probes native-compartment capacity and leaves explicit separation; it does not invent a connective bridge or establish a clinical muscle subdivision. Use a separate candidate archive. Require finite cells, whole/part self checks, unchanged source hashes, and complete current ownership-forest contact comparison. Any new canal, rectum, skin, floor or other visible-owner contact rejects the trial. Never overwrite the previous verified candidate.


The separate deep-edge capacity trial moves 96 paired controls by at most 2.493143 mm with exact X/Z and fixed lower parts/interfaces. All four tissue shells retain zero detected self-crossings and no invalid cells. Complete current ownership screening introduces no new visible contact owners. The trial remains uninstalled and anatomically unaccepted: it is not proof of a finished upper attachment or the appropriate TA98/puborectalis subdivision. Preserve both source and trial archives while current user-directed female scalp work is completed. Evidence: external-sphincter-deep-edge-lift.json.


## Deep-edge replay and actual section

The current-atlas bounded lift is independently replayed from saved Float32 fields: original faces, region ownership and radial coordinates stay exact; subcutaneous/superficial controls and shared region boundaries remain fixed. All whole/part shells are closed orientable genus-one tissue with additive volumes, zero detected self-crossings and finite cells. Native skin/canal/rectum/puborectalis hashes remain exact; their complete surface checks pass. Native canal and rectal vertices/centroids remain outside the muscle material.

Complete triangle separation falls from 2.871045 to 0.404608 mm at the deep outer surface. Actual assembled and lateral section views in existing Blender show the improved proximity; prior scenes are retained and the main file unsaved. Stored field/index/ownership storage remains 36,868 bytes, with no added controls, faces, UVs or textures. The whole sphincter/attachment is still unaccepted and uninstalled: a remaining gap, proximal morphology, perineal/anococcygeal connective supports, source puborectalis correspondence and clinical dimensions need resolution. A geometry pass is not anatomical completion.

Evidence: external-sphincter-deep-edge-lift-validation.json, external-sphincter-deep-lift-interface.json, external-sphincter-deep-lift-interface-render.json; outputs/female-deep-sphincter-lift-interface-section.png. Related native puborectalis, pubococcygeus, anal canal, rectum, pudendal/inferior-anal route and rectal vessels remain part of the coupled review scope.


## Native deep interface inventory and clipping plan

| TA98 | Parent | Term | Current native sources | Classification |
|---|---|---|---|---|
| A04.5.04.007 | A04.5.04.002 | m. puborectalis | LOCAL:female:puborectalis:236 | macroscopic muscle; no distinct descendants in the current authority |

Clip each existing proximal annular tissue triangle against the projected native puborectalis underside. Enforce the actual height correspondence (native underside above the candidate, within 12 mm) as additional affine clipping constraints. Record exact native face/barycentric ownership and check projected overlap before treating any patch as a possible shared tissue interface. This operation creates a diagnostic surface correspondence only; unsupported spans must remain explicit, not become a fabricated complete attachment ring. Native floor, canal, rectum, skin and vessels remain fixed. No clinical attachment acceptance or main-model geometry change follows from this calculation.

### Independent projected-interface replay, 2026-10-07

The independent barycentric validator initially rejected the map: nearly duplicate arithmetic vertices produced 0.0087823 mm² of duplicate ownership. Diagnostic clipping now removes consecutive arithmetic duplicates within 1e-14 metres; this does not weld or modify any anatomical source vertices. A native-scale regression covers the offending polygon.

The regenerated correspondence contains 355 convex patches covering 178.951482 of 262.180867 mm² (68.254974%). Independent source reconstruction, candidate-triangle inclusion, height/gap constraints, same-cell disjoint ownership and nearest competing native-face checks pass. Residual overlap is about 2.08e-12 mm² and incorrectly nearer ownership about 8.13e-14 mm², both numerical roundoff. Candidate and native catalogue/packed geometry hashes are verified before replay. Five clipping regression tests pass.

This proves only the recorded projected geometry correspondence. It does not establish a clinical insertion footprint, resolve the unsupported 31.745026%, or accept the proposed proximal muscle morphology. No assembly geometry, resource buffers, installed anatomy or camera behavior changes. Blender is listed as running, but its UI connection timed out twice during this continuation; it was not relaunched, reset or saved. The new footprint has not been visually reviewed in Blender.

Evidence: `scripts/validate-female-exact-deep-footprint.py`, `server/projected-clip.test.py`, and local `external-sphincter-exact-deep-footprint-validation.json`. Next inspect the verified support against the actual native underside and define the supported attachment geometry while preserving canal, rectum, skin, lower subdivisions and related floor/vessel constraints.

Actual browser assembly and underside inspection is now available in `/female-sphincter-interface-review.html` using exact current source geometry (4,724 review triangles including diagnostic overlays, 331,560 JSON bytes). Green polygons are projected correspondence only. Native puborectalis has an uneven course; the constructed annular rim is not an established clinical insertion. No source reference tissue is copied into the body catalogue.

A separate source-guided rim trial moves 37 paired upper controls by at most 2 mm, using 11 fully ray-supported radial stations and verified projected height constraints. Lower tissue/ownership and radial coordinates stay fixed; all four shells retain their construction topology/self checks and zero invalid cells. **The trial is rejected for installation:** it creates 14 new triangle contacts with the native superior rectal vein. The current-forest screen, not the isolated muscle checks, rejects the fit. See `external-sphincter-source-guided-deep-rim.json`; candidate is preserved separately as `source-guided-deep-rim-trial.npz`. Its clinical attachments, rendered modified surface, complete lumen and independent replay remain unaccepted/unverified. User eye-color correction interrupted this work; no further female anatomy was installed.
