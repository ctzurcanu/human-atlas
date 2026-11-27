# TA98 mesh coverage implementation

**Full TA98 mesh coverage is unfinished.** This records geometry integrated into the local model, not a completeness certificate.

Recovered 127 TA98 selections and appended 365 source meshes. All 56 phalanges are registered. Another 260 textured region pieces replace 112 original bone/cartilage pieces across all 56 phalanges. Original catalogue order is preserved; superseded pieces remain suppressed recovery records. The repaired gastric geometry is unchanged.

The updated viewer catalogue is `.local-models/ta98-male.json`. The Blender deliverable is `outputs/atlas-ta98-male-review.blend` (6,291 mesh objects; 6,176 active). New unregistered surfaces are hidden by default and shown by explicit selection.

Verification: all 6,291 meshes have finite positions/UVs, valid indices, image textures and valid references. Blender was reopened and checked for counts, TA98 collections, 36 packed textures and preserved gastric coordinates. Source bone surfaces and cross-sections were rendered; the joined exported bone was also inspected in the dark-theme viewer.

| TA98 | Structure | Active mesh pieces | Representation limitation |
| --- | --- | ---: | --- |
| A14.2.01.006 | optic nerve [II] | 2 | bilateral optic nerve surfaces |
| A14.1.09.304 | neocortex | 252 | assembly of named neocortical surface pieces; no invented cortical layers |
| A03.1.08.004 | posterior atlanto-occipital membrane | 1 | named membrane surface |
| A05.1.03.004 | incisor tooth | 8 | available permanent incisor teeth |
| A05.1.03.005 | canine tooth | 4 | available permanent canine teeth |
| A05.1.03.006 | premolar tooth | 8 | available permanent premolar teeth |
| A05.1.03.007 | molar tooth | 8 | available permanent molar teeth; third molars absent |
| A14.1.09.572 | white substance of telencephalon | 2 | bilateral cerebral white matter surfaces |
| A14.1.06.012 | brachium of inferior colliculus | 2 | bilateral inferior collicular brachia |
| A14.1.06.013 | brachium of superior colliculus | 2 | bilateral superior collicular brachia |
| A12.2.07.066 | artery of postcentral sulcus | 2 | bilateral named arterial surfaces |
| A05.6.04.002 | terminal ileum | 11 | distal ileal tissue pieces; exclude the separate ileocecal junction |
| A02.3.01.005 | costal cartilage | 14 | available first through seventh costal cartilages; remaining cartilages not supplied by this group |
| A12.3.07.011 | ascending lumbar vein | 3 | ascending lumbar veins; same anatomy occurs in two TA98 routes |
| A12.3.09.004 | ascending lumbar vein | 3 | same ascending lumbar vein assembly, reused without duplicated geometry |
| A12.2.12.030 | right branch | 8 | named right hepatic arterial tree |
| A12.2.12.035 | left branch | 7 | named left hepatic arterial tree |
| A12.3.12.002 | right branch | 9 | named right portal venous tree |
| A12.3.12.005 | left branch | 8 | named left portal venous tree |
| A12.3.12.003 | anterior branch | 2 | anterior branch of the right portal vein |
| A12.3.12.004 | posterior branch | 3 | posterior branch of the right portal vein |
| A12.3.12.011 | medial branches | 2 | medial branches of the left portal vein |
| A12.3.12.009 | lateral branches | 2 | lateral branches of the left portal vein |
| A12.2.09.043 | anterior branch | 2 | anterior ulnar recurrent arterial branch |
| A12.2.09.044 | posterior branch | 2 | posterior ulnar recurrent arterial branch |
| A15.2.07.014 | check ligament of lateral rectus muscle | 2 | named lateral rectus check ligaments |
| A12.2.07.039 | intermediomedial frontal branch | 2 | intermediomedial callosomarginal arterial branches |
| A12.2.07.040 | posteromedial frontal branch | 2 | posteromedial callosomarginal arterial branches |
| A12.2.07.057 | middle temporal branch | 2 | middle temporal branches of the middle cerebral artery |
| A12.2.01.110 | middle lobar artery | 8 | available branches of the right middle lobar artery |
| A12.2.01.102 | superior lobar arteries | 24 | right superior lobar arterial tree |
| A12.2.01.203 | superior lobar arteries | 15 | left superior lobar arterial tree |
| A12.2.01.113 | inferior lobar arteries | 22 | right inferior lobar arterial tree |
| A12.2.01.214 | inferior lobar arteries | 22 | left inferior lobar arterial tree |
| A12.2.01.115 | basal part | 13 | right basal segmental arterial tree |
| A12.2.01.216 | basal part | 15 | left basal segmental arterial tree |
| A12.2.12.033 | anterior segmental artery | 2 | anterior segmental hepatic arteries |
| A12.2.12.034 | posterior segmental artery | 3 | posterior segmental hepatic arteries |
| A12.2.12.037 | medial segmental artery | 2 | medial segmental hepatic arteries |
| A12.2.12.038 | lateral segmental artery | 3 | lateral segmental hepatic arteries |
| A12.2.12.041 | pancreatic branches | 4 | available pancreatic branches of the splenic artery |
| A12.2.12.079 | superior segmental artery | 2 | superior segmental renal arteries |
| A12.2.12.082 | inferior segmental artery | 2 | inferior segmental renal arteries |
| A12.2.12.084 | posterior segmental artery | 3 | posterior segmental renal arteries |
| A12.2.01.103 | apical segmental artery | 8 | right apical segmental arterial tree |
| A12.2.01.107 | posterior segmental artery | 6 | right posterior segmental arterial tree |
| A12.2.01.104 | anterior segmental artery | 9 | right anterior segmental arterial tree |
| A12.2.01.114 | superior segmental artery | 9 | right superior segmental arterial tree |
| A12.2.01.116 | anterior basal segmental artery | 3 | right anterior basal segmental arterial tree |
| A12.2.01.117 | lateral basal segmental artery | 4 | right lateral basal segmental arterial tree |
| A12.2.01.119 | posterior basal segmental artery | 6 | right posterior basal segmental arterial tree |
| A14.1.06.303 | accessory nuclei of oculomotor nerve | 2 | bilateral named accessory oculomotor nuclei |
| A14.1.07.203 | declive [VI] | 1 | named cerebellar declive surface |
| A03.6.10.508 | dorsal intercuneiform ligament | 2 | bilateral dorsal intercuneiform ligament surfaces |
| A12.3.09.009 | renal veins | 2 | bilateral renal vein surfaces |
| A14.2.01.015 | motor root | 2 | bilateral trigeminal motor roots |
| A14.2.01.013 | sensory root | 2 | bilateral trigeminal sensory roots |
| A03.2.04.002 | alar ligaments | 1 | source alar ligament assembly; side completeness requires review |
| A12.3.10.021 | inferior rectal veins | 2 | bilateral inferior rectal veins |
| A12.3.10.018 | middle rectal veins | 2 | bilateral middle rectal veins |
| A12.2.15.006 | lateral sacral arteries | 2 | bilateral lateral sacral arteries |
| A12.3.05.008 | meningeal veins | 2 | source meningeal veins |
| A12.3.05.038 | deep temporal veins | 2 | bilateral deep temporal veins |
| A12.3.07.025 | posterior spinal veins | 1 | posterior spinal venous surface |
| A12.3.07.024 | anterior spinal veins | 1 | anterior spinal venous surface |
| A12.2.07.048 | anterolateral central arteries | 2 | bilateral lenticulostriate arterial surfaces |
| A04.2.04.001 | infrahyoid muscles | 8 | source infrahyoid muscle assembly |
| A04.2.03.001 | suprahyoid muscles | 11 | source suprahyoid muscle assembly |
| A04.1.01.001 | extra-ocular muscles | 14 | source extraocular muscle assembly |
| A12.2.12.058 | ileal arteries | 1 | available ileal arterial surface |
| A12.3.12.023 | pancreaticoduodenal veins | 3 | source pancreaticoduodenal venous assembly |
| A12.3.12.020 | ileal veins | 1 | available ileal venous surface |
| A12.3.10.007 | obturator veins | 2 | bilateral obturator veins |
| A02.4.08.001 | carpal bones | 16 | bilateral assembly of all sixteen carpal bones |
| A02.5.09.001 | tarsal bones | 14 | bilateral assembly of all fourteen tarsal bones |
| A12.2.07.042 | paracentral branches | 6 | paracentral callosomarginal branches |
| A12.2.07.044 | precuneal branches | 2 | precuneal pericallosal branches |
| A12.2.07.061 | superior terminal branches | 22 | source superior terminal middle cerebral arterial branches |
| A12.2.07.055 | inferior terminal branches | 8 | source inferior terminal middle cerebral arterial branches |
| A12.2.08.021 | pontine arteries | 4 | bilateral pontine arteries |
| A12.2.08.022 | medial branches | 2 | medial pontine arterial branches |
| A12.2.08.023 | lateral branches | 2 | lateral pontine arterial branches |
| A12.2.07.095 | anterior temporal branches | 2 | anterior temporal lateral occipital arterial branches |
| A14.2.01.008 | superior branch | 2 | bilateral superior oculomotor nerve branches |
| A14.2.01.009 | inferior branch | 2 | bilateral inferior oculomotor nerve branches |
| A14.2.01.027 | long ciliary nerves | 2 | bilateral long ciliary nerves |
| A14.3.02.005 | short ciliary nerves | 2 | bilateral short ciliary nerves |
| A06.1.01.006 | nasal cartilages | 5 | source nasal cartilage assembly |
| A12.3.09.002 | inferior phrenic veins | 6 | inferior phrenic venous assembly |
| A12.3.01.012 | anterior veins of right ventricle | 2 | source anterior cardiac veins |
| A12.2.09.031 | superficial palmar branch | 2 | bilateral superficial palmar branches of radial arteries |
| A12.2.09.055 | deep palmar branch | 2 | bilateral deep palmar branches of ulnar arteries |
| A12.2.16.018 | saphenous branch | 2 | bilateral saphenous branches of descending genicular arteries |
| A12.2.16.028 | ascending branch | 2 | bilateral ascending branches of lateral circumflex femoral arteries |
| A12.2.16.030 | transverse branch | 2 | bilateral transverse branches of lateral circumflex femoral arteries |
| A14.1.09.258 | body | 2 | bilateral source bodies of the fornix |
| A02.5.17.001 | metatarsals [I - V] | 10 | all ten metatarsal bones, reused as one bilateral assembly |
| A02.2.02.001 | cervical vertebrae [C I - C VII] | 7 | seven cervical vertebrae, retained as individually selectable components |
| A02.3.01.001 | ribs [I - XII] | 24 | all twenty-four source rib bones |
| A02.3.01.002 | true ribs [I - VII] | 14 | bilateral ribs I through VII |
| A02.3.01.003 | false ribs [VIII - XII] | 10 | bilateral ribs VIII through XII |
| A04.3.02.213 | rotatores cervicis | 2 | bilateral cervical rotator muscle components |
| A04.3.02.211 | rotatores lumborum | 2 | bilateral lumbar rotator muscle components |
| A04.3.02.301 | interspinales | 2 | source interspinalis muscle assembly |
| A04.3.02.401 | intertransversarii | 4 | source intertransversarius muscle assembly |
| A02.4.10.002 | proximal phalanx | 50 | all 10 named proximal phalanges of the hand, including both sides |
| A02.4.10.003 | middle phalanx | 40 | all 8 named middle phalanges of the hand, including both sides |
| A02.4.10.004 | distal phalanx | 40 | all 10 named distal phalanges of the hand, including both sides |
| A02.4.10.001 | phalanges | 130 | all 28 source phalanges of both hands; original bone and articular cartilage surfaces retained |
| A02.5.18.002 | proximal phalanx | 50 | all 10 named proximal phalanges of the foot, including both sides |
| A02.5.18.003 | medial phalanx | 40 | all 8 named middle phalanges of the foot, including both sides |
| A02.5.18.004 | distal phalanx | 40 | all 10 named distal phalanges of the foot, including both sides |
| A02.5.18.001 | phalanges | 130 | all 28 source phalanges of both foots; original bone and articular cartilage surfaces retained |
| A14.1.08.106 | stria medullaris of thalamus | 2 | bilateral named striae medullares of the thalami |
| A06.1.01.003 | dorsum of nose | 2 | named dorsum-of-nose surface pieces; a surface feature, not an extra solid nose |
| A05.8.01.033 | hepatic segmentation | 8 | original eight named hepatic segment tissue meshes, reused without copying geometry |
| A05.8.01.037 | left liver | 4 | TA98 left hepatic part assembled from its II, III, IV and posterior I descendants; this is not the traditional left lobe |
| A05.8.01.038 | left lateral division | 2 | left lateral division: original segment II and III tissue |
| A05.8.01.041 | left medial division | 1 | left medial division: original segment IV tissue |
| A05.8.01.043 | posterior liver | 1 | posterior hepatic part: original segment I tissue |
| A05.8.01.045 | right liver | 4 | right functional hepatic part: original segments V, VI, VII and VIII; not the traditional right lobe |
| A05.8.01.046 | right medial division | 2 | right medial division: original segment V and VIII tissue |
| A05.8.01.049 | right lateral division | 2 | right lateral division: original segment VI and VII tissue |
| A05.8.01.063 | anterior branch | 2 | anterior branches of the right hepatic duct: named superior and inferior source tributaries |
| A05.8.01.064 | posterior branch | 3 | posterior branches of the right hepatic duct: named superior and inferior source tributaries |
| A05.8.01.066 | lateral branch | 3 | lateral branches of the left hepatic duct: named superior and inferior source tributaries |
| A05.8.01.067 | medial branch | 3 | medial branches of the left hepatic duct: named superior and inferior source tributaries |

## Installed phalanx regions

- 36 proximal and middle phalanges now have separate base, shaft and head surfaces, using shared triangle-edge intersections. No overlapping faces or internal caps were added.
- All 56 original whole bones passed connected, genus-zero surface topology and non-adjacent triangle intersection checks. All 168 measured cross-sections of the subdivided bones were closed.
- All 56 exported assemblies passed exact shared coordinate/normal checks, duplicate/degenerate face checks and UV/material validation.
- Collar interfaces remain **provisional**, following the source articular patches. Base and shaft selectors now include all 56 source instances, with **provisional anatomical boundaries**. Head/trochlear selectors retain their proximal/middle scope.
- The two trochlear selectors expose only the distal articular cartilage surface; they do not establish separately modeled bony trochleae.
- No marrow, cortical thickness, growth plates or microscopic structures were fabricated.
- Existing saved selection/hidden URLs expand suppressed original pieces into their replacement surfaces.

Proof: `outputs/ta98-phalanx-region-review.png` (regions and measured cross-sections), `outputs/ta98-phalanx-joined-view.png` (live viewer), and `reports/ta98-male/phalanx-export-validation.json` (actual exported buffers).

## Distal phalanx additions

20 distal bones now have independently selectable base and shaft surfaces, plus 14 provisional tuberosity surfaces. Six source tips remain explicitly unresolved: both fifth fingers and both third/fourth toes. The source surfaces are preserved; no absent terminal anatomy was fabricated. The 80 new pieces replace 40 original pieces using shared coordinates, normals and UVs.

Rendered verification includes every distal bone from two angles and 60 closed measured cross-sections. Proof files: `outputs/ta98-distal-phalanx-region-review.png`, `outputs/ta98-distal-phalanx-region-review-side.png`, and `outputs/ta98-distal-phalanx-joined-view.png`.

Blender image maps now multiply the original tissue colors and source vertex colors, matching the viewer material behavior instead of replacing those colors with neutral texture pixels.

## Rebuilding

`npm run build:ta98-male-review` now includes recovery, texturing, reviewed phalanx partitions and hierarchy generation. `npm run build:ta98-phalanges` validates existing partitions rather than subdividing them again. Python needs NumPy and SciPy; set `ATLAS_GEOMETRY_PYTHON` if necessary. Blender is installed locally; `ATLAS_BLENDER` can select its executable.

Flat-color source materials have neutral illustrative image textures and UVs. These textures preserve source colors and are **not histological reconstructions**. Existing source image maps remain intact.

## Liver assemblies and additional recovered selections

Fourteen additional TA98 selectors now expose actual named meshes: eight hepatic segment assemblies, four right/left biliary branch groups, bilateral striae medullares and the dorsum-of-nose surface. Only two new source pieces were necessary; the other selections reuse existing tissue meshes. No extra geometry was copied for a second hierarchy path.

The liver assemblies use parenchymal segment meshes. Functional hepatic divisions are kept distinct from traditional surface lobes and their source vessel/duct groups. All 69 liver terms, parents, mesh sources and missing representations were recorded in `reports/ta98-male/liver-source-assembly-plan.json` before integration. Missing capsule, serosa, lobules and microscopic structures remain unresolved. Source tissue boundaries and duct registration are not certified complete.

Recovered concepts and new region concepts now expose their own TA98 identifiers in search, WebMCP and the inspector. The inspector displays the supplied coverage limitation instead of assigning the first child mesh's identifiers to an assembly.

## Remaining work

Recovering source assemblies does not establish independently modeled wall layers, microscopic structures, subregions or complete arterial trees. Partial molar and costal-cartilage groups are explicitly retained as subsets. New source surfaces require alignment and anatomical review.

The remaining structures need source geometry or anatomically justified modeling, followed by rendered and cross-section inspection. Generic containing surfaces and unrelated vessels must not be relabeled as missing tissue. Local review assets have not been deployed.

## Pancreas regional meshes

Five textured source-surface regions are installed locally: head, uncinate process, neck, body and tail. Head selection includes its uncinate child. All source selectors expand to the replacement assembly; the old whole surface is suppressed. Original material and interpolated UVs are retained.

All 19 TA98 pancreas terms were inventoried before subdivision in `reports/ta98-male/pancreas-parts-plan.json`. Interfaces remain provisional: proximal SMV landmarks, the clinical left-aortic-border body/tail convention, and the inferior medial source hook. This does not certify individual anatomical boundaries. Named body surfaces and borders, pancreatic notch, omental eminence, duct sphincter, accessory pancreas and pancreatic islets still require dedicated geometry.

The installed assembly has 6,318 triangles and 3,161 unique vertices: one closed connected genus-zero surface, no duplicated faces, exact shared normals and positions, preserved source-face areas and UVs. Blender detected zero nonadjacent triangle intersections; five measured section contours are closed. Front/back region renders and actual section contours are in `outputs/ta98-pancreas-*`.

`npm run build:ta98-pancreas` validates an existing installation or prepares a fresh candidate. After a fresh build, inspect its renders before `python3 scripts/run-ta98-pancreas.py --install-reviewed`. Full TA98 coverage remains unfinished.

| Pancreas TA98 code | Installed textured selection |
|---|---|
| A05.9.01.002 | Head: head surface plus uncinate child |
| A05.9.01.003 | Uncinate process surface |
| A05.9.01.005 | Neck surface |
| A05.9.01.006 | Body surface |
| A05.9.01.014 | Tail surface |

These five additions are provisional regions, separate from the 127 source-recovered selectors.

## Cardiac source additions

An individually traced, textured anterolateral papillary head (BodyParts3D FJ2418) is installed as a **partial** anterior left papillary representation, A12.1.04.007. It is registered against four corresponding primary papillary pieces; maximum fitted bounding landmark residual is 0.284 mm. This fitting error does not certify anatomical attachments. Source FMA7264 also contains FJ2429, whose terminology conflicts with the primary inferior left papillary naming; no whole-muscle coverage is inferred from the group.

The source contained a redundant 13-vertex/22-face internal closed shell. All its vertices tested inside the main component along three independent parity ray directions. Removing it preserves the visible outer surface. The final mesh has 754 triangles, one closed genus-zero component, no detected self-intersections, no duplicate faces, and three closed measured section contours.

Five TA98 selections added: ventricles (A12.1.00.012), atria (.017), papillary muscles (.022, partial), anterior left papillary (.04.007, head only), posterior left papillary (.04.008, primary inferior piece with provisional correspondence). Whole heart includes the new component. Other cardiac terms remain unresolved as recorded in the complete heart inventory, `reports/ta98-male/heart-source-parts-plan.json`. Full TA98 coverage is unfinished.

Use `npm run test:ta98-cardiac` for exported binary and selector verification. `npm run build:ta98-cardiac` prepares candidates after a fresh build; review front/back and measured-section renders before installation with `python3 scripts/run-ta98-cardiac-selections.py --install-reviewed`.

## Complete thyroid replacement and regions

The primary thyroid was an incomplete unilateral open surface: 84 vertices, 125 triangles and 41 boundary edges. It is now suppressed and replaced by the complete evaluated local reference thyroid, partitioned into right/left lobes, isthmus and a visible pyramidal lobe. No spatial transform was applied to the reference source. The original whole-thyroid selection and all original part selectors expand to the four replacement regions. Local reference redistribution remains restricted pending license review and deployment authorization.

All nine thyroid TA98 terms were inventoried before modification in `reports/ta98-male/thyroid-parts-plan.json`. The three new term selections are A11.3.00.002 (bilateral lobes), .003 (isthmus), and .004 (pyramidal lobe). Lateral boundaries and pyramidal base follow measured upper silhouette valleys of the source, sampled at 0.1 mm; tissue interfaces remain provisional. The outer surface does not establish accessory glands, fibrous capsule, stroma or microscopic lobules.

The complete assembled export contains 3,124 triangles and 1,564 unique vertices, one connected closed genus-zero surface, exact shared normals/positions, preserved source-face areas and continuous UVs. Blender detected zero nonadjacent intersections; five measured sections are closed and their contours rendered for review. Both front/back rendered geometry and measured section images were inspected before installation.

`npm run test:ta98-thyroid` validates installed binary geometry and selectors. After a reset, `npm run build:ta98-thyroid` stages the whole source and region render gates; install a visually reviewed candidate with `python3 scripts/run-ta98-thyroid.py --install-reviewed`. Full TA98 geometry coverage remains unfinished.

## Evaluated parathyroid gland replacement

The four primary superior/inferior parathyroid glands each had 9 vertices and 8 triangles. Their stable primary IDs now reference complete evaluated local neck-source geometry: 66 vertices and 128 triangles per gland, 512 triangles total. Reference coordinates and packed normals are preserved exactly; each gland has UVs and the existing image texture. No registration or shared tissue partition was introduced.

All four TA98 parathyroid terms and their source selections were inventoried before replacement in `reports/ta98-male/parathyroid-parts-plan.json`, against the authoritative chapter 11. The parent A11.4.00.001 now selects the four superior/inferior glands rather than the unrelated coarse supplemental assembly. A11.4.00.004 accessory glands remain missing; macroscopic gland surfaces do not establish microscopic layers.

Each replacement is closed, connected and genus zero, with zero detected self-intersections and three closed measured sections. Posterior and oblique renders were inspected against the complete thyroid source, followed by rendered mid-gland section contours. Source placement remains that of the local reference. Review assets and redistribution remain local. Use `npm run test:ta98-parathyroid` to verify exact installed coordinates, indices, normals, topology, UVs, textures and parent selection. Full TA98 coverage remains unfinished.

## Additional named vascular and lymphatic source selections

Seven missing selections were recovered from explicitly named local reference sources: hand intercapitular veins (A12.3.08.024), foot intercapitular veins (.11.019), lumbar veins (.07.012 and .09.003), arteries of the pterygoid canal (A12.2.05.082 and .06.006), and retropharyngeal nodes (A13.3.00.025). The original primary geometry is reused for five distinct assemblies; repeated terminology routes do not duplicate meshes. Only the two hand intercapitular vein meshes were appended, with original source coordinates and illustrative texture UVs. Catalogue total: 6,293 parts. Local source licensing restrictions still apply.

Both new hand-vein sources have three connected tubular components, 60 open edges, no nonmanifold edges, and **36 self-intersections per side** near a junction. Their renders were inspected in context and isolation; they are explicitly marked as defective review candidates, hidden by default. Their adjoining venous interfaces and repair are unfinished. The staged exact self-union experiment removes detected intersections but preserves only 48 of the original 60 boundary vertices, so it has **not been installed or accepted**. See `reports/ta98-male/intercapitular-source-review.json` and the complete hand-vein inventory/interface plan in `intercapitular-repair-plan.json`. No vessel wall layers are claimed, and the added candidates do not prove complete anatomical coverage.

### Provisional original hand-vein wall selection

The hand source contains seven **edge-connected face components**, rather than three clean wall objects: three open tube surfaces and four closed unlabelled node primitives. Exact-coordinate welding joins some nodes to the tubes at individual vertices, producing misleading vertex-component counts and invalid vertex topology. One node overlaps the tube and accounts for the detected intersections. The exact Boolean experiments are rejected: both lose an original end ring and introduce nonmanifold edges.

Male TA98 now uses the three original open tube surfaces per side, retaining every selected source face and its positions, normals and texture coordinates. The four closed nodes are omitted **provisionally**: their anatomical identity remains unresolved; they are not inferred to be venous valves or wall layers. All original source records and binary files remain available in the backup/metadata. This is a component subset, not a certified complete vein assembly.

Each installed wall selection has 1,610 triangles, three components with Euler characteristics 0, 0 and 1, five unchanged 12-vertex end rings, zero nonmanifold edges, no duplicate or degenerate faces, and zero detected self-intersections. Three measured section planes per side produce closed contours. Isolated wall, surrounding vein context and measured-section renders were inspected before installation. Exact end-ring coordinates and retained normals/UVs are verified by `npm run test:ta98-intercapitular`. Adjoining venous continuity remains unresolved, so the selections remain local review geometry hidden by default. Full TA98 coverage is unfinished.

## Individually selectable lung segment surface patches

Twenty-five lung-related TA98 source-surface selections now reuse the existing textured primary lung geometry, including 18 individual segment patches, bilateral superior/inferior lobes, right middle lobe, left lingula, both lung assemblies and the segment aggregate. Clicking a numbered patch selects its segment concept rather than the complete source lobe. No geometry, material, UVs, normals, coordinates or original binary offsets were changed. Full lung terminology (50 terms) and source/interface definitions were inventoried before mapping in `reports/ta98-male/lung-source-selections-plan.json`.

These are **component subsets**, not certified whole segment volumes. The source patches numbered 1–9 are interpreted in their specific lobe context. Both inferior lobes use an ambiguous `Lung-1` source label for a posterior basal candidate; those two S X mappings remain provisional. The left medial basal S VII has no separate source patch and remains missing. Lung-base source patches are unallocated hilar surfaces; they are not asserted to be anatomical bases, wall layers or microscopic lung lobules.

Original rendered front/back surfaces and the isolated right middle lobe were inspected. Lobe source integrity remains unfinished: all assemblies have open edges, the inferior right lobe has nonmanifold edges and seven detected intersections; no lung repair completion is claimed. Existing lung geometry is reused, while coarser local alternatives remain reversible review sources. `npm run test:ta98-lungs` checks 25 exact selections, segment click identities, missing/provisional flags and preservation of every original textured binary field. Full TA98 coverage remains unfinished.


## Six spinal study alternatives — 2026-10-01

Added six individually selectable textured reference surfaces: lateral vestibulospinal tract, tectospinal tract, rubrospinal tract, anterior/posterior spinocerebellar tracts, and gracile fasciculus. Their TA98 entries now link to local study alternatives. Each contains two continuous, open-ended source surfaces. Exact source coordinates, indices and normals were retained; existing primary 5 mm study slabs remain unchanged. Alternatives are hidden by default and revealed for explicit TA98 inspection.

The full spinal term inventory, assemblies, interfaces and entity distinctions are in `reports/ta98-male/spinal-tract-source-plan.json`. The six surfaces have no detected nonadjacent triangle intersections or nonmanifold edges; each component has Euler characteristic zero and legitimate end-ring boundaries. Measured sections at y=1.20, 1.30 and 1.40 m have closed contours. Longitudinal montage and section renders were inspected. Generated UVs use existing illustrative tissue maps; they do not show axons or histology.

**Coverage remains partial.** These inherited fitted extrusions have repeated cross-sectional profiles. Anatomical extent, origins, terminations, decussations, connectivity and variation by spinal level are unverified. No attachment, wall closure or new anatomical partition was invented. They cannot establish anatomically complete pathways. Unverified redistribution rights keep these additions local.

`npm run test:ta98-spinal-tracts` checks source fidelity, source preservation, texture mapping, partial-status selections and hidden defaults. The full male build/test commands include this step.


## Accessory tendon sheath and bursa source surfaces — 2026-10-01

Added 12 textured named reference source envelopes for six TA98 selections: flexor pollicis longus, flexor carpi radialis, extensor pollicis longus and extensor carpi ulnaris tendon sheaths; finger fibrous sheaths; and gluteus medius trochanteric bursae. Both sides are represented. Eight finger-source components are carried by two meshes. These selections are explicitly partial: the surfaces do not certify synovial lining, complete wall thickness, cavity continuity, tendon fit or registration. The fibrous-sheath selection does not yet include the thumbs.

The original thumb sheath meshes each have 14 detected nonadjacent triangle intersections. They were reviewed and excluded from installation, with identifiers and reasons preserved. Twelve imported surfaces have no detected triangle self-intersections, open/nonmanifold edges, or open measured section contours. Extensor pollicis longus source shells each have Euler characteristic zero (a tunnel); the other single source envelopes have Euler characteristic two. None is labeled as an anatomically complete wall repair. Source positions, normals and indices are retained exactly. All existing male geometry remains unchanged; new local alternatives are hidden by default.

Inventory of all 53 applicable hand/lower-limb accessory terms, parentage and source correspondence: `reports/ta98-male/accessory-sheath-source-plan.json`. Source diagnostics and sections: `accessory-sheath-source-review.json`. Palmar, dorsal and isolated bursa renders were inspected. The mid-bursa plane initially coincided with source vertices; noncoincident section planes close. No extra solid was substituted for microscopic lining or a cavity.

`npm run test:ta98-accessory-sheaths` verifies source fidelity, partial TA98 links, textures, preserved existing geometry and defect exclusions. Added to the full male build/test commands. Local assets remain undeployed pending redistribution rights.


## Thumb fibrous sheath source repair — 2026-10-01

Added both previously excluded thumb envelopes after a complete source-surface repair. Original source shells each had 14 detected nonadjacent triangle intersections. Exact self-intersection resolution removed the overlapping fold; a local six-iteration relaxation (factor 0.2, 4 mm radius) removed a tiny residual inner protrusion identified by dense sections. Repair parameters and the full accessory structure inventory were recorded before modification in `reports/ta98-male/thumb-sheath-repair-plan.json`.

Each installed thumb envelope has 695 vertices and 1,390 triangles, one connected closed component, Euler characteristic zero (one retained through passage), no open/nonmanifold edges, no duplicate/degenerate triangles, and no detected nonadjacent triangle intersections. All 44 measured contours close: 43 slices show separate inner/outer contours, and the terminal slice crosses the original oblique end lip as one contour. Front/back surfaces and measured sections were inspected. The original and repaired exceptional sections are compared in `outputs/ta98-thumb-sheath-exceptional-sections.png`.

Sampled vertex-to-surface distance is at most 0.18 mm from repaired vertices to source, and 0.91 mm from removed/intersecting source vertices to the repaired surface. These are sampled distances, not a certified whole-surface Hausdorff bound. No new anatomical attachment, wall layer or microscopic lining is claimed. Anatomical registration, source wall thickness and tendon fit remain provisional; TA98 A04.8.04.015 stays `component-subset`.

The fibrous-sheath selection now includes both finger envelopes and both repaired thumb envelopes. Every existing binary position, normal, index and UV remains unchanged; existing finger limitation text was updated to reflect inclusion of the thumbs. Original defect exclusions remain as provenance. New thumb pieces have illustrative fiber textures, individual side concepts, and hidden defaults. Local-only redistribution restrictions remain.

`npm run test:ta98-thumb-sheaths` verifies reviewed candidate hashes, installed binary fidelity, topology, retained passage sections, textures and selections. Included in the full build/test commands.

## Accessory gaps: extensor digiti minimi sheath and subfascial prepatellar bursa

Four previously absent reference meshes are installed in Male TA98 (two sides for each structure):

- A04.8.04.013: `REF:upper-limb-left:218:0`, `REF:upper-limb-right:218:0` (1,852 triangles each).
- A04.8.05.015: `REF:lower-limb-left:403:0`, `REF:lower-limb-right:396:0` (992 triangles each).

Correspondences checked against the IFAA TA98 muscle chapter. The reference spelling “Subfacial” is retained in source concept identities and corrected to “Subfascial” in displayed mesh names. The full 100-term accessory inventory, parents, current representations, missing entries, assembly, and interface limitations were recorded before installation in `reports/ta98-male/accessory-gap-plan.json`.

Each mesh retains exact source coordinates, triangle indices and packed normals. UVs and an illustrative image finish are added. No existing geometry is modified. Each source has one closed connected component, Euler characteristic 2, zero nonmanifold edges and zero detected nonadjacent triangle intersections. Three section planes per side form closed contours. Front/back rendered geometry and all 12 sections were inspected in `outputs/ta98-accessory-gap-montage.png` and `outputs/ta98-accessory-gap-sections.png`.

**Partial coverage:** The hand source is C-shaped in section, with a longitudinal slit and joined inner/outer lips. It does not establish a full enclosed tendon sheath lumen. The prepatellar bursa is a flattened envelope without a resolved inner cavity, wall layers or synovial lining. Neither has verified registration against Male TA98 tendons/fascia. These limitations remain in mesh/concept metadata and TA98 selections; this is not a claim of anatomically complete repair. Both bilateral assemblies use `component-subset` and are hidden by default. Original detailed geometry is preserved. Review assets remain local because redistribution rights are unverified.

Reproduce: `python3 scripts/ta98-accessory-gap-sources.py prepare`, then Blender `--background --python scripts/review-ta98-sheath-sources.py -- --gap`. Inspect actual renders and sections, record the review findings, then use the `install` and `validate` modes. Installation refuses missing or mismatched reviewed preparation. `npm run test:ta98-accessory-gaps` and the full review test include source fidelity, UV/texture, selection and preservation checks. Male TA98 now contains **6,317 mesh records**; this count does not imply complete TA98 coverage.

## Explicit repeated TA98 routes: six formerly unresolved entries

Six entries now select the existing textured primary surfaces:

- A14.2.01.084 and A14.2.01.118: both chorda tympani nerves, shared identities across lingual and facial/intermediate routes.
- A05.1.01.104 and A05.2.01.003: the same named soft-palate source patch in oral and pharyngeal routes.
- A13.3.00.038: the same bilateral retropharyngeal-node meshes already represented under A13.3.00.025.
- A03.1.03.003: the gingival source surface already represented under A05.1.01.108.

These are explicit reviewed correspondences, not a global rule that repeated names imply identical structures. IFAA chapters 03, 05, 13 and 14 define the terminology; the nervous-system chapter explicitly identifies the repeated chorda tympani entries. The plan records affected terms, descendants, parents, source meshes, missing details and interface limitations before modifying selection metadata.

No mesh was added, replaced, transformed or partitioned. Geometry, normals, UVs, image maps, binary chunks and default visibility remain unchanged. The new `TA98ROUTE:<code>` concepts refer to the original mesh identities; all six entries remain `component-subset`, not anatomically complete. The model still contains 6,317 mesh records.

Actual dark front/back renders and component cross-sections were inspected. Both nerve tubes are closed genus zero. Each retropharyngeal-node mesh has two closed components, totaling four node envelopes; all twelve component sections close. Gingiva has three closed components, including a detached 34-triangle component of unresolved identity. Soft palate remains a visibly incomplete open patch with 46 boundary edges and three open sections. No self-intersections or nonmanifold edges were detected in these sources. No repair of the open palate, detached gingival component, nerve connections, or microscopic anatomy is claimed. These findings are included in the term/concept limitations and review report.

Artifacts: `reports/ta98-male/repeated-route-plan.json`, `repeated-route-review.json`, `repeated-route-validation.json`; `outputs/ta98-accessory-routes-montage.png` and `ta98-accessory-routes-sections.png`.

Reproduce review with `python3 scripts/ta98-repeated-route-selections.py prepare`, Blender `--background --python scripts/review-ta98-sheath-sources.py -- --routes`, then `python3 scripts/render-ta98-selection-review.py --report reports/ta98-male/repeated-route-review.json --prefix ta98-accessory-routes-`. Recording visual findings is required before installation. `npm run test:ta98-repeated-routes` validates the source-preserving selections and runs within the full review suite. All assets remain local.

## Soft palate: replace the incomplete coarse patch with the intact source envelope

The prior soft-palate source review found 46 boundary edges and three open sections in the 216-triangle primary patch. A denser original reference surface (`REF:head:719:0`) was available. Before replacement, the plan inventoried all seven applicable entries under both soft-palate routes, existing sources, unresolved descendants, entity distinctions, assembly and interfaces.

The denser whole source has 1,884 welded vertices and 3,764 triangles, one connected component, Euler characteristic 2, zero open/nonmanifold edges, zero detected nonadjacent self-intersections, and no topological through-tunnel. Actual dark front/back renders and three measured sections were inspected; all section contours close. This improves the original surface integrity without manufacturing caps or subdividing tissue. It does not prove anatomical wall thickness or interface fit.

`ZA:Soft palate` now contains the exact original dense source positions, indices and packed normals, with fresh illustrative planar UVs and the existing image material. Its stable identifier, concept links and default visibility are retained. No duplicate palate is appended and no other geometry is modified. The model remains 6,317 mesh records. Both A05.1.01.104 and A05.2.01.003 still use `component-subset`: complete descendant structures, muscle/aponeurosis separation, microscopic layers and adjoining hard-palate/pharyngeal interfaces remain missing or unverified. Their current limitations describe the closed replacement, superseding the earlier open-patch findings without deleting that history.

Files: `reports/ta98-male/palate-source-plan.json`, `palate-source-review.json`, `palate-source-validation.json`; `.local-models/ta98-review/palate-source.bin`; dark rendered geometry and sections in `outputs/ta98-palate-source-montage.png` and `ta98-palate-source-sections.png`.

Reproduce staging with `python3 scripts/prepare-ta98-palate-source.py`; review with Blender `--background --python scripts/review-ta98-sheath-sources.py -- --prepared .local-models/ta98-review/palate-source-prepared.json --report reports/ta98-male/palate-source-review.json --output-prefix ta98-palate-source-`; generate inspectable section panels with `python3 scripts/render-ta98-selection-review.py --report reports/ta98-male/palate-source-review.json --prefix ta98-palate-source-`. Record rendered review findings before running `install-ta98-palate-source.py`. The dedicated validation checks exact source geometry/normals, closed genus-zero topology, duplicate/degenerate faces, texture UVs, stable route selections and preservation of all other geometry/materials/visibility.

Earlier preservation validators now permit only this explicitly recorded successor: exact original-part metadata, exact installed-part metadata and SHA-256 of its binary must match. Arbitrary changes to existing anatomy still fail. The full review suite includes dedicated replacement validation. All review assets stay local.

## Gingiva: dense original envelope and verified tooth-context inspection

`ZA:Gingiva` now uses the exact denser original reference surface `REF:head:228:0`: **24,596 triangles replace 4,168**. The old stable primary identifier, existing image material and selection identities remain. No extra whole gingiva or overlapping copied tissue is appended. The model remains 6,317 mesh records.

Before replacement, the checked-in plan recorded five applicable entries under the oral/periodontal gingiva routes, their parents, existing source meshes and unresolved representations. Whole gingiva, individual papillae, margin, sulcus, periodontal attachment and microscopic layers remain distinct. No coarse envelope is assigned as all those descendants.

The original source contains three closed connected components: 12,424 and 12,008 triangle jaw envelopes, plus one detached **164-triangle component of unresolved anatomical identity**. Each has Euler characteristic 2. Zero open/nonmanifold edges and zero nonadjacent self-intersections were detected. Nine cross-section planes through the individual component height ranges produce closed contours. Whole front/back surfaces, all sections, front/oblique/back context with all 28 primary tooth envelopes, and the isolated tiny component were actually rendered and inspected. No coordinate transform was used. The tooth context appears coherent, but exact periodontal surface attachment/clearance is not certified. The tiny component is retained with an explicit unresolved label; it is not reassigned to an invented papilla or removed without evidence.

Fresh illustrative UVs use the existing image map. Final positions, triangle indices and packed normals equal the dense source exactly. A03.1.03.003 and A05.1.01.108 select this same local geometry and remain `component-subset`. Individual papillae, microscopic layers and periodontal interfaces still need anatomical evidence and geometry work. This is a source-detail improvement, not completion of gingival TA98 descendants.

Artifacts: `reports/ta98-male/gingiva-source-plan.json`, `gingiva-source-review.json`, `gingiva-context-review.json`, `gingiva-source-validation.json`; `outputs/ta98-gingiva-source-montage.png`, `ta98-gingiva-source-sections.png`, `ta98-gingiva-context-front.png`, `ta98-gingiva-context-oblique.png`, `ta98-gingiva-context-back.png`, and `ta98-gingiva-small-component.png`.

Reproduce staging with `prepare-ta98-gingiva-source.py`; review the original surface using Blender `review-ta98-sheath-sources.py -- --prepared .local-models/ta98-review/gingiva-source-prepared.json --report reports/ta98-male/gingiva-source-review.json --output-prefix ta98-gingiva-source-`. `review-ta98-gingiva-context.py` uses all primary `· Teeth` envelopes and adjusts the near clipping plane for the tiny-component inspection. Review findings are required before `install-ta98-gingiva-source.py`. The dedicated validator checks exact geometry/normals, three closed genus-zero components, duplicate/degenerate faces, sections, UVs, textures, stable routes and preservation of other anatomy. Full review tests include it. Earlier preservation checks allow only both precisely recorded source replacements, with matching original/installed part metadata and binary SHA-256. All assets remain local.

## Ligament source review: six partial selections, two explicit exclusions

Reviewed fifteen existing local meshes against the TA98 articular chapter and recorded the full articular terminology inventory before selection changes. Six entries now have explicit partial selections referencing twelve existing textured meshes:

- A03.0.00.007: bilateral forearm interosseous-membrane examples; the general term is not fully represented by these examples.
- A03.2.01.004: bilateral intertransverse source groups (24 coarse closed block-like components per side; anatomical level completeness unverified).
- A03.2.08.004 and A03.2.08.005: bilateral anterior and lateral sacrococcygeal source surfaces, with open attachment ends and unverified shared bone interfaces.
- A03.6.08.015: bilateral coarse closed patellar ligament sources; attachment fit remains unverified.
- A03.6.10.006: bilateral anterior tibiotalar source surfaces, not the entire medial ankle ligament; open ends and attachment fit remain unresolved.

No mesh was added, transformed, replaced or partitioned. Surface geometry, textures, UVs, chunks and default visibility are preserved. Each new concept remains `component-subset`. Actual dark front/back surface renders and first/middle/last measured component sections were inspected; closure was computed at three levels for every connected component. No anatomical completeness or repair of open source ends is asserted.

Two sources require further work and are explicitly excluded from automatic term correspondence:

1. A03.2.04.006: the source labeled transverse ligament of atlas renders as a cruciform-looking assembly with longitudinal extensions. TA98 distinguishes the cruciate ligament, longitudinal bands and transverse ligament. Assigning this whole source to the transverse component would be unsupported. Full part inventory is present; dissection and component correspondence remain pending.
2. A03.6.07.002: right/left zona orbicularis surfaces have 10/6 detected self-intersections. Their whole surfaces need repair and renewed section/geometry review. Closed edges and genus-one topology alone do not establish sound surfaces.

Prior unreviewed selections for these two entries are retained as provenance in `ligamentSourceSelections.excludedPriorSelections` and removed from current audit selections. Source concepts carry `ta98ExcludedTerms`; the hierarchy builder propagates the exclusion by exact source concept identity across source/combined models so name matching cannot silently restore rejected correspondences. Geometry remains accessible as an unverified local source.

Artifacts: `reports/ta98-male/ligament-selection-plan.json`, `ligament-selection-review.json`, `ligament-selection-validation.json`; `outputs/ta98-ligament-source-montage.png`, `ta98-ligament-source-sections.png`, and `ta98-ligament-section-sample.png`.

Stage with `python3 scripts/ta98-ligament-selections.py prepare`; review with Blender `review-ta98-sheath-sources.py -- --prepared .local-models/ta98-review/ligament-selection-prepared.json --report reports/ta98-male/ligament-selection-review.json --output-prefix ta98-ligament-source-`. Render inspection panels with `render-ta98-selection-review.py`. Record visual findings before installation. `npm run test:ta98-ligament-selections` validates the explicit source references, preserved geometry/materials/visibility and unresolved exclusions; it is included in the full review suite. Review assets remain local and the full TA98 goal remains incomplete.

## Zona orbicularis source rings — 2026-10-01

Installed two repaired, textured whole source rings for **A03.6.07.002** in Male TA98, retaining source identifiers and default visibility. The hip-joint terminology, parents, source matches and missing representations are recorded in `reports/ta98-male/zona-repair-plan.json`. No capsule layer, attachment, or microscopic structure was inferred from this envelope.

The first exact self-Boolean attempt produced nonmanifold branching edges despite reporting zero triangle intersections; it was rejected. Independent local relaxation repaired the left fold but caused additional branching on the right. The source rings have bijective reflected vertex correspondence and identical triangle connectivity. The accepted right repair therefore uses only homologous left fold patches, reflected with the measured source x offset of -0.00064094 m; the remaining right surface is retained. The left fold uses ten local smoothing iterations with factor 0.2 in the original intersection neighborhoods (4 mm radius). Right patch transfer is restricted to 5 mm neighborhoods of the original intersection points. This correction is provisional anatomy, not proof of ligament-to-capsule or femoral-neck fit.

Both installed surfaces: **718 vertices, 1,436 triangles, one connected component, Euler characteristic zero, no open/nonmanifold edges, no detected nonadjacent triangle intersections**. Each has 67 closed measured sections. Front, back, lateral and passage renders and representative fold-region sections were inspected. The ring passage is retained; there are no added caps. Sampled vertex-to-surface distances are at most 0.819 mm from candidate to source and 0.747 mm from source to candidate; these are not certified Hausdorff distances.

Stable part IDs were replaced rather than duplicated. Existing image materials were retained with newly generated UVs and repaired normals. The repaired selection is provided by `TA98ZONA:A03.6.07.002`; the original source concept exclusions remain to prevent reintroducing the unmodified defective meshes from other models. The rejected cruciform/transverse source correspondence also remains excluded. Original zona defect provenance remains in ligament metadata. Dedicated validation checks exact installed geometry/normals against the inspected candidate, packed binary hash, edge topology, textures, hierarchy selection, unchanged visibility and other mesh preservation.

Review images: `outputs/ta98-zona-repair-montage.png`, `outputs/ta98-zona-repair-sections.png`. Candidate and Blender assets remain local; redistribution rights are unverified. Full TA98 coverage remains incomplete.

## Head muscle source cleanup and selections — 2026-10-01

Male TA98 now has explicit reviewed selections for eight additional terms: alar part of nasalis (A04.1.03.011), orbicularis oculi (.013), depressor supercilii (.019), auricularis anterior (.020), superior (.021), posterior (.022), levator labii superioris alaeque nasi (.032), and lateral pterygoid (A04.1.04.006). These reuse **20 existing textured mesh records**. Orbicularis oculi uses the dense orbital and palpebral components rather than the coarse whole source. Lateral pterygoid selects both superior and inferior heads on both sides. These selections remain partial anatomical coverage: ciliary bundle, deep/lacrimal eyelid part, fascia, microscopic architecture, insertion footprints and tissue registration are not established.

The six original ear-muscle envelopes contained collapsed faces or duplicated zero-thickness fins. Whole-envelope cleanup removes all copies of those fins and collapsed faces; it preserves every remaining surface coordinate. Anterior ear muscles now have 48 triangles each, superior 66 left/64 right, posterior 38 each. Each is one closed genus-zero component with no detected triangle intersections, duplicate/degenerate faces or nonmanifold edges.

Both inferior pterygoid heads originally had Euler characteristic one despite closed edges: one vertex had two separate link cycles, indicating a surface pinch. The repair separates those incident fans by a provisional **0.05 mm** shift toward their respective neighboring fan centroids. All other coordinates remain exact. Each repaired head has 326 triangles and one closed genus-zero component. This is a local surface correction, not a certified anatomical fissure or attachment reconstruction.

Eight whole surfaces were repaired before any proposed partition. Full head-muscle TA98 inventory, parents, source matches and missing representations are recorded in `reports/ta98-male/head-muscle-plan.json`. Whole bilateral renders, individual repaired-source closeups and all 60 measured sections were inspected. All 20 selected surfaces have closed measured sections and no detected nonadjacent triangle intersections or nonmanifold edges. The eight replacements additionally pass single-cycle vertex-link validation. Original materials, stable identifiers, locations and default visibility remain; repaired meshes receive new UVs and normals. Unmodified meshes and previous repairs are preserved by explicit binary-backed replacement records.

Scripts: `repair-ta98-head-muscles.py`, `install-ta98-head-muscle-repairs.py`, `ta98-head-muscle-selections.py`, `validate-ta98-head-muscle-repairs.py`, `run-ta98-head-muscles.py`. Review outputs include `ta98-head-repair-closeups-0.png`, `ta98-head-repair-closeups-1.png`, and `ta98-head-muscle-repaired-sections.png`. All source assets stay local pending license verification. Full TA98 textured geometry coverage remains incomplete.

## Limb muscle assemblies and local surface corrections — 2026-10-01

Seven additional whole-muscle TA98 selections now use dense existing source heads: pronator teres (A04.6.02.025), flexor carpi ulnaris (.030), flexor digitorum superficialis (.033), extensor carpi ulnaris (.045), flexor pollicis brevis (.055), flexor hallucis brevis (A04.7.02.057), and adductor hallucis (.060). The selections reference **48 existing textured tissue pieces**, forming **28 bilateral source heads**. The 21 applicable parent/head TA98 terms, original source matches, missing representations and planned interfaces are recorded in `reports/ta98-male/limb-muscle-plan.json`.

The review assembles each source head in its exact original coordinates, including both muscle and tendon partitions where supplied. All 20 two-piece heads had identical paired muscle/tendon boundary coordinates and normals. Original face ownership is preserved in the exported partitions; no duplicated overlapping whole muscle, attachment marker or invented microscopic layer was added.

Original surfaces had detected local triangle intersections in ten heads and pinched vertices in two ulnar extensor-carpi-ulnaris heads. Whole-head local corrections precede partition re-export. Fold correction uses ten distance-weighted relaxation iterations (factor 0.15), restricted to a 2 mm neighborhood of the detected crossings, or 8 mm for the lateral flexor-hallucis-brevis folds. Coincident vertex fans in the ulnar extensor heads are separated by a provisional 0.05 mm shift toward their own incident fan centroid. All 22 affected tissue partitions retain their source face counts and identities, use normals computed on the complete corrected head, and preserve exact shared seam coordinates/normals when re-exported.

All 28 inspected source heads now have no detected nonadjacent triangle intersections or nonmanifold edges, with single-cycle vertex links. **116 measured sections** close, including additional cuts immediately around corrected crossings and pinches. Whole front/back assemblies, individual corrected-head closeups and focused measured sections were inspected. Dedicated validation proves that installed tissue pieces reassemble exactly into the inspected whole-head candidate, preserves tissue ownership/counts, verifies every shared muscle/tendon seam and normal, checks texture UVs and binary hashes, and preserves all unrelated geometry and default visibility.

**Existing source tunnels remain unresolved.** Several original heads have nonzero genus. Their anatomy has not been established; local fold/pinch correction does not constitute a complete organ repair or justify filling those tunnels. Attachment footprints, inter-head union, fascia and microscopic architecture remain unverified. All seven selections are marked partial (`component-subset`), not complete coverage. Source models remain local pending redistribution-rights verification.

Scripts: `repair-ta98-limb-muscles.py`, `install-ta98-limb-muscles.py`, `validate-ta98-limb-muscles.py`, `render-ta98-limb-muscle-repairs.py`. Reports: `limb-muscle-interfaces.json`, `limb-muscle-repaired-review.json`, `limb-muscle-validation.json`. Whole renders, individual closeups and focused sections are in `outputs/ta98-limb-*`. Full TA98 textured mesh coverage remains incomplete.


## Axial muscles: provisional installed geometry

Ten existing source records were corrected and installed with stable IDs and image textures; nine TA98 assemblies use 28 unique textured records. Assemblies cover supplied digastric bellies/tendons, anterior cervical intertransversarii, iliocostalis, multifidus, semispinalis capitis, thoracic rotatores, puborectalis, and psoas minor. These are partial representations, not complete anatomical coverage.

Four multifidus cervical/lumbar surfaces received local fold/pinch corrections. Six thoracic multifidus, semispinalis and rotatores surfaces required provisional 0.15 mm voxel envelopes and simplification. Duplicate zero-thickness simplification fins were removed. All 28 installed candidates match the inspected geometry exactly; their coordinate-welded edges have two incident faces, vertex links are single cycles, and no duplicate faces, degenerate triangles or detected nonadjacent intersections remain. All 337 recorded cross-sections close. Source-to-envelope sampled vertex deviation is below 1.1 mm; this is not a Hausdorff bound or anatomical certification.

Actual whole assemblies, closeups and representative measured sections were inspected in `outputs/axial-whole-page-*.png`, `outputs/axial-close-page-*.png` and `outputs/axial-section-page-*.png`. Jagged source attachment ends, source tunnels, disconnected fragments, subvoxel detail and continuous anatomical joins remain unresolved. The six reconstructed envelopes are explicitly provisional. Existing detached puborectalis fragments are retained without invented identities. A vague posterior cervical intertransversarii source is excluded from lateral-only correspondence because medial/lateral membership is unverified.

Evidence: `reports/ta98-male/axial-muscle-plan.json`, `axial-muscle-repaired-review.json`, `axial-muscle-candidate-validation.json`, and `axial-muscle-validation.json`. Recheck with `npm run test:ta98-axial-muscles` using the bundled NumPy Python runtime. Local model only; redistribution rights remain unverified.

Integration verification: the full model review suite passed; hierarchy selectors and descendant sets passed for 7,474 terms across seven models. The rebuilt Blender file reopened with 6,317 objects and 38 packed textures; all ten axial replacements retained their exact exported coordinates. See `reports/ta98-male/axial-muscle-integration-validation.json`. Full TA98 coverage remains incomplete.


## Masseter whole and regional assemblies

The missing whole-masseter TA98 selector now uses all eight existing textured muscle/tendon pieces. Both superficial and deep regional selectors were explicitly reviewed and installed. Four coordinate-welded whole source envelopes are closed genus-one surfaces, have no detected nonadjacent intersections or pinched vertices, and contain no duplicated or degenerate faces. All 12 measured sections close. Muscle/tendon seams have identical coordinates and quantized normals; no caps or duplicated tissue faces were introduced. Whole region renders and all recorded sections were inspected. Original source tunnels, attachment ends and internal architecture remain unresolved; no complete anatomical repair is claimed. Geometry, materials and default visibility are preserved exactly. Evidence: `reports/ta98-male/masseter-plan.json`, `masseter-review.json`, and `masseter-validation.json`.


## Major muscle assemblies and local geometry corrections

Six previously unresolved whole-muscle TA98 selectors now use 60 existing textured tissue pieces: biceps brachii, triceps brachii, iliopsoas, quadriceps femoris, biceps femoris and gastrocnemius. Thirty complete source-head envelopes were inspected. Eight heads had local surface intersections corrected in 2 mm neighborhoods; both medial gastrocnemius heads had four zero-thickness duplicate faces removed each. The latter artifacts had escaped edge-valence and triangle-intersection checks but failed the dedicated duplicate/vertex-link validator. All ten corrections were installed as twenty muscle/tendon partitions with stable source IDs and image maps.

All thirty whole heads now pass coordinate-welded edge-valence, single-cycle vertex-link, duplicate/degenerate face and detected nonadjacent intersection checks. All 114 remaining recorded sections close. Every source coordinate outside the recorded repair neighborhoods/artifact exemptions is retained exactly. Tissue pieces reassemble to the inspected whole heads; all shared seam coordinates and quantized normals match. Face ownership and counts are unchanged except the eight explicitly recorded zero-thickness artifact faces. The original source tunnels and attachment ends remain anatomically unresolved; no complete organ repair or microscopic representation is claimed. Rectus femoris straight/reflected subhead boundaries and common tendon interfaces remain unresolved.

Rendered whole assemblies, all repaired-head closeups and every measured section were inspected in `outputs/major-whole-page-*.png`, `outputs/major-close-page-*.png` and `outputs/major-section-page-*.png`. Evidence is recorded in `reports/ta98-male/major-muscle-plan.json`, `major-muscle-repaired-review.json`, `major-muscle-candidate-validation.json`, `major-muscle-partition-validation.json` and `major-muscle-validation.json`. Local only; redistribution rights remain unverified. Full TA98 mesh coverage remains incomplete.

Integration checks passed: full model review suite, all 7,474 hierarchy selectors/descendant sets across seven models, and Blender reopening with 6,317 objects and 38 packed textures. All twenty major-muscle replacement coordinates survived export exactly. Evidence: `reports/ta98-male/major-muscle-integration-validation.json`.


## Regional lumbar/cervical muscles and rib elevators

Six previously missing TA98 selectors now use twelve existing textured source records: whole lateral lumbar intertransversarii, its dorsal and ventral parts, medial lumbar intertransversarii, cervical multifidus, and levatores costarum. The cervical multifidus correspondence uses the authoritative TA98 cervicis/colli synonym; it reuses the already corrected source geometry. The dorsal/ventral source pieces are reused for their parent group without duplicating tissue geometry. Alternate lateral lumbar source FJ1547/FJ1547M and ambiguous posterior cervical sources are excluded from these selections.

Four rib-elevator records had duplicate or collapsed zero-thickness faces removed: four faces from each long source and two from each short source. Removing the short-source artifacts resolved the invalid vertex link observed despite edge closure. No surviving surface coordinates were moved and no faces or caps were added. Dedicated validation reconstructs the expected surviving coordinate-face multiset from the originals and compares it exactly to the installed candidates.

All twelve selected source surfaces pass coordinate-welded edge-valence, single-cycle vertex-link, duplicate/degenerate-face and detected nonadjacent-intersection checks. All 252 recorded sections close. Whole group renders, four cleaned source closeups and representative sections from every selected source were inspected. Original source tunnels, small fragments, vertebral/rib level completeness, attachment footprints and microscopic architecture remain unresolved; no complete anatomical repair is claimed. Textures, other geometry and default visibility are preserved. Evidence: `reports/ta98-male/regional-muscle-plan.json`, `regional-muscle-repaired-review.json`, `regional-muscle-candidate-validation.json` and `regional-muscle-validation.json`. All assets remain local pending redistribution rights. Full TA98 coverage remains incomplete.

Integration: the full model review suite and 7,474 term selectors/descendant sets across seven models passed. Blender reopened with 6,317 mesh objects and 38 packed textures; the four regional replacements retained exact exported coordinates. See `reports/ta98-male/regional-muscle-integration-validation.json`.


## Scalp muscles and epicranial aponeurosis

The previously missing epicranius and occipitofrontalis selections now reuse the supplied bilateral frontal bellies, occipital bellies, temporoparietalis and epicranial aponeurosis. All six applicable TA98 entries are explicitly reviewed; their twelve textured source records are reused without duplicate geometry. The full descendant inventory, source correspondence, tissue distinctions and planned interfaces are recorded in `reports/ta98-male/scalp-muscle-plan.json`.

Four whole source envelopes had local nonadjacent triangle intersections: one in each occipital muscle and three in each aponeurosis half. Ten distance-weighted relaxation iterations, factor 0.15 and radius 2 mm, corrected these crossings before re-export. Source faces and their tissue ownership remain unchanged. All coordinates outside the recorded repair neighborhoods remain exact. No caps or new anatomical subdivisions were invented.

All eight source envelopes pass coordinate-welded edge-valence, single-cycle vertex-link, duplicate/degenerate-face and detected nonadjacent-intersection checks. All forty measured sections close and were inspected, along with front/back whole assemblies and corrected-source closeups. Every muscle/tendon partition pair has exact shared boundary coordinates and quantized normals. Installed pieces reassemble exactly to the inspected whole candidates; image maps, finite UVs, unrelated geometry and default visibility are preserved.

**Anatomical completion remains unresolved.** Source genus is retained, including unidentified tunnels. Closed section contours alone do not establish anatomical correctness. The original jagged attachment ends, bilateral aponeurosis midline continuity and true muscle/aponeurosis joins require further verification. These supplied outer surfaces do not establish scalp layers or microscopic architecture. Six selections remain partial; no complete whole-organ repair is claimed.

Evidence: `scalp-muscle-repaired-review.json`, `scalp-muscle-candidate-validation.json`, `scalp-muscle-partition-validation.json`, `scalp-muscle-validation.json` and `scalp-muscle-integration-validation.json` in `reports/ta98-male/`. Recheck with `npm run test:ta98-scalp-muscles` using the bundled NumPy Python runtime. The full review suite, all 7,474 term selectors across seven models and Blender reopening passed. The review file has 6,317 mesh objects and 38 packed textures; four replacement coordinates survived export exactly. All assets remain local pending redistribution-rights verification. Full TA98 textured mesh coverage remains incomplete.


## Omohyoid superior/inferior bellies and transitional tendon regions

Six actual textured surface regions now replace the two undivided omohyoid source surfaces: bilateral superior belly (A04.2.04.004), inferior belly (.005), and intermediate-tendon examples (A04.0.00.045). The whole muscle remains the continuous assembly of all three regions on each side. Existing digastric intermediate-tendon examples are retained. The original whole mesh records are suppressed and all concept/selection references are remapped, avoiding overlapping visible whole meshes.

Before subdivision, all applicable authoritative TA98 descendants, parent relationships, source matches, missing boundary representations and interfaces were inventoried in `reports/ta98-male/omohyoid-plan.json`. Both original whole sources were inspected front/back and in six measured sections. Each is a single closed genus-zero surface without detected nonadjacent triangle intersections. No surface repair was required.

The narrow source transition is bounded by provisional planes at y=1.429123 and y=1.441123 m. These are observed geometric approximations, not verified histological myotendinous boundaries. Each source triangle is conformingly clipped using shared barycentric coordinates. Its area is conserved independently to within 1e-13 square metres before float32 storage. No new caps, wall layers or duplicate faces are added. All original vertices remain and texture coordinates/normals are interpolated consistently at new seam vertices. The reconstructed whole surfaces remain closed genus zero with valid single-cycle vertex links and no detected intersections. Every shared interface coordinate, quantized normal and UV is identical between adjacent regions.

Both colored assemblies and all fourteen measured sections, including cuts immediately around each interface, were inspected. Installed binary meshes match the reviewed candidates exactly. Original materials and default visibility are retained; unrelated geometry is preserved through explicit original/installed source records and exact recorded region additions. The superior and inferior belly terms are marked `reconstructed-region-boundary-provisional`; no anatomical completion claim is made. Fascial sling, precise attachment footprints, verified tendon tissue identity and microscopic architecture remain unresolved.

Evidence: `omohyoid-review.json`, `omohyoid-partitioned-review.json`, `omohyoid-validation.json` in `reports/ta98-male/`. Scripts: `prepare-ta98-omohyoid.py`, `partition-ta98-omohyoid.py`, `render-ta98-omohyoid-regions.py`, `install-ta98-omohyoid.py`, `validate-ta98-omohyoid.py`. Source terminology: https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/04%20TA98%20EN.htm . Anatomical study consulted for the two-belly/tendon relationship and variation: https://pmc.ncbi.nlm.nih.gov/articles/PMC2664130/ . Exact boundaries are our provisional inference from the source surface, not supplied by these references. Local model only; full TA98 textured mesh coverage remains incomplete.

Integration verification passed: the full model review suite, all 7,474 hierarchy selectors/descendant sets across seven models, and Blender reopening with 6,323 objects and 38 packed textures. All six region coordinates survived export exactly. See `reports/ta98-male/omohyoid-integration-validation.json`.


## Diaphragm: exact central tendon extraction

All sixteen TA98 diaphragm entries and their parentage, source matches, missing tissue representations and planned adjoining interfaces were inventoried before geometry changes in `reports/ta98-male/diaphragm-plan.json`. The whole muscle/tendon assembly has 9,592 triangles, one closed component, Euler characteristic -2 (genus two), no detected nonadjacent triangle intersections, no duplicate/degenerate faces and valid single-cycle vertex links. Its three measured sections close. Whole front/back, superior/inferior and posterior-oblique renders were inspected. Two enclosed source passages and the posterior arch are retained without filling or caps; their exact registration with vessel/oesophageal walls remains unverified.

The named tendon source contains two disconnected components. Treating all its faces as central tendon would incorrectly include the inferior crural attachment component. Exact graph-component extraction creates a central tendon mesh of 820 triangles and a separate inferior tendon component of 369 triangles. No coordinates, normals, UVs or source face ownership change. The central component supplies the previously missing A04.4.02.013 selection; it is marked partial. The inferior component remains in the whole diaphragm assembly and is not asserted to represent any particular arcuate ligament.

Both open tissue partitions were inspected independently from front/back and in six measured sections; open contours at tissue interfaces are expected and explicitly recorded. Their coordinate-face multiset equals the original tendon exactly. Reassembly with the unchanged muscle reproduces the complete original whole surface. All 163 shared muscle/tendon coordinates and quantized normals match exactly. Original image materials and UVs are retained. The original combined tendon record is suppressed and all concept/selection references are remapped, so no overlapping visible source remains. All unrelated geometry and default visibility are preserved.

The costal/sternal/lumbar muscle regions, right/left crural boundaries, medial/lateral/median arcuate ligaments and phrenico-oesophageal ligament remain unresolved. Existing supplied tissue surfaces do not establish these ligaments or histology. Openings and anatomical spaces are retained as spaces, not separate solids. Full TA98 coverage is incomplete; local assets remain undeployed pending redistribution-rights verification.

Evidence: `diaphragm-review.json`, `diaphragm-tendon-components-review.json` and `diaphragm-tendon-validation.json` in `reports/ta98-male/`. Scripts: `prepare-ta98-diaphragm.py`, `partition-ta98-diaphragm-tendon.py`, `render-ta98-diaphragm-views.py`, `render-ta98-diaphragm-components.py`, `install-ta98-diaphragm-tendon.py`, `validate-ta98-diaphragm-tendon.py`. Run `npm run test:ta98-diaphragm-tendon` with the bundled NumPy Python runtime.

Integration verification passed: full model suite, all 7,474 hierarchy selectors and descendant sets across seven models, and Blender reopening with 6,325 mesh objects and 38 packed textures. Both tendon components retained exact exported coordinates. See `reports/ta98-male/diaphragm-tendon-integration-validation.json`.


## Diaphragm: five textured muscle regions

The sixteen-entry diaphragm inventory and planned whole assembly/interfaces preceded this partition. Five actual image-textured patches now provide provisional selections for lumbar part (A04.4.02.002), right crus (.003), left crus (.004), costal part (.008), and sternal part (.009). All 8,403 original source muscle triangles belong to exactly one patch. Positions, triangle winding, quantized normals, UVs and material remain exact; no surface subdivision, smoothing, caps or additional tissue was introduced. The original whole muscle record is suppressed and its concept references remapped to these regions. The lumbar selection includes both crural muscle patches, the remaining lumbar patch and the existing inferior tendon component.

Anatomical boundaries are provisional geometric inferences along existing source edges, recorded explicitly in the plan and selection metadata. Counts are right crus 867, left crus 574, remaining lumbar 2,110, costal 4,798 and sternal 54 triangles. The lumbar assignment has two large components plus a one-triangle fragment; the costal assignment has one large component and 7-, 26-, and 16-triangle fragments. These fragments are unresolved partition assignments, not identified anatomical structures. The source-edge boundaries remain jagged. Individual crural tendon allocation, exact attachments, microscopic architecture, arcuate ligaments and phrenico-oesophageal ligament remain unresolved.

Superior, inferior and posterior oblique colored assemblies, each isolated region front/back and all thirty measured region sections were inspected. Open patch contours represent shared tissue interfaces. Reassembly with both installed tendon patches reproduces the original one-component, closed 9,592-triangle diaphragm exactly, including its two source passages. All 163 tendon/muscle seam coordinates and quantized normals agree. Passage registration with neighbouring organs remains unverified. No anatomical completeness claim follows from the topology checks.

Evidence: `reports/ta98-male/diaphragm-muscle-candidate-validation.json`, `diaphragm-muscle-regions-review.json`, and `diaphragm-muscle-validation.json`. Run `npm run test:ta98-diaphragm-muscle` with the bundled NumPy Python runtime. Original geometry, materials and default visibility are preserved. Full TA98 mesh coverage remains incomplete; this model stays local.

Integration checks passed: full model suite, 7,474 hierarchy selectors/descendant sets across seven models, and Blender reopening with 6,330 mesh objects and 38 packed textures. All five installed region coordinates survived export exactly. See `reports/ta98-male/diaphragm-muscle-integration-validation.json`.

## Pelvic diaphragm and levator ani: partial textured assemblies

Before source selection, all applicable pelvic diaphragm descendants and adjoining external anal sphincter descendants were inventoried with their authoritative TA98 parentage, source matches and unresolved representations in `reports/ta98-male/pelvic-muscles-plan.json`. Two formerly empty aggregate terms now select actual existing textured meshes: levator ani (A04.5.04.002) uses bilateral pubococcygeus, pubo-analis, puborectalis and iliococcygeus; pelvic diaphragm (A04.5.04.001) additionally uses bilateral coccygeus. These are partial assemblies, not newly manufactured tissue. Source origin annotation patches and the external anal sphincter are excluded from muscle tissue selections. Urogenital hiatus remains a space rather than an invented solid.

All ten source envelopes were inspected front/back, together in superior/inferior/oblique colored views, and in every one of 129 measured cross-sections. Their existing triangle coordinates/winding, textures, UVs and metadata are unchanged. Each passes closed edge, single-cycle vertex-link, duplicate/degenerate face and measured-section checks, with no detected nonadjacent self-intersections. No exact coordinate-identical faces occur twice within either assembly selection. These checks do not establish anatomical tissue joins or eliminate intermesh volume overlaps.

Visible gaps and irregular source envelopes remain. Pubococcygeus and pubo-analis each retain genus-three topology; iliococcygeus retains genus-two topology. The puborectalis sources retain genus-eleven main envelopes with 11 left and 22 right components; coccygeus retains a main envelope plus an eight-triangle component on each side. These tunnels and detached components have no verified anatomical identity. Cross-source registration, puboperinealis/puboprostaticus boundaries, anococcygeal body/tendon/raphe, fascial layers and precise attachment footprints remain unresolved. No continuous watertight whole pelvic floor or complete TA98 coverage claim is made.

Evidence: `reports/ta98-male/pelvic-muscles-plan.json`, `pelvic-muscles-review.json`, `pelvic-muscles-validation.json`; scripts `prepare-ta98-pelvic-muscles.py`, `render-ta98-pelvic-assembly.py`, `install-ta98-pelvic-muscles.py`, `validate-ta98-pelvic-muscles.py`. Run `npm run test:ta98-pelvic-muscles` with the bundled NumPy Python runtime. These assets stay local. Terminology reference: https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/04%20TA98%20EN.htm .

Integration passed: full model review suite, all 7,474 term selectors and descendant sets across seven models, and Blender reopening with 6,330 objects and 38 packed textures. See `reports/ta98-male/pelvic-muscles-integration-validation.json`.

## Masticatory muscles, spinotransversales and triceps surae

The full authoritative TA98 descendant inventories and planned assemblies/interfaces preceded changes in `reports/ta98-male/muscle-groups-plan.json`. Three previously empty group terms now select existing actual textured source geometry: masticatory muscles (A04.1.04.001), spinotransversales (A04.3.02.101), and triceps surae (A04.7.02.043). The first uses bilateral superficial/deep masseter, temporalis, superior/inferior lateral pterygoid and medial pterygoid; the second splenius capitis/cervicis; the third medial/lateral gastrocnemius and soleus. Plantaris, source origin/insertion annotation patches and fascial surfaces are excluded. Group selections reuse 36 source mesh records, avoiding newly duplicated visible geometry.

Whole-surface inspection found one small self-intersection in each splenius capitis and one two-fan vertex pinch in each medial pterygoid. Four complete source envelopes were corrected before installation: splenius relaxation is confined to recorded 2 mm neighborhoods; pterygoid fans receive an explicitly provisional 0.05 mm separation. Six existing source partition records are re-exported. Original muscle/tendon face ownership is preserved, no faces/caps added or removed, shared positions/quantized normals agree exactly, and source coordinates outside the recorded correction neighborhoods are retained. Original material image textures remain; corrected tissue pieces use illustrative planar UVs, not a histological texture claim.

All 22 complete envelopes now pass duplicate/degenerate face, closed edge, single-cycle vertex-link and nonadjacent self-intersection checks. All 74 measured sections, including cuts around correction locations, close and were inspected, together with individual muscle views, repair closeups and complete group views. Existing nonzero-genus source tunnels are retained with their anatomical meaning unverified. Bone contacts, muscle-to-muscle overlaps, tendon tissue allocation, complete fascia/internal architecture and exact attachments remain unresolved. This is local geometry and selection progress, not full anatomical certification or complete TA98 coverage.

Evidence: `reports/ta98-male/muscle-groups-plan.json`, `muscle-groups-repaired-review.json`, `muscle-groups-candidate-validation.json`, `muscle-groups-partition-validation.json`, and `muscle-groups-validation.json`. Scripts: `prepare-ta98-muscle-groups.py`, `repair-ta98-muscle-groups.py`, `render-ta98-muscle-group-repairs.py`, `install-ta98-muscle-groups.py`, and associated candidate/partition/installed validators. Run `npm run test:ta98-muscle-groups` with the bundled NumPy Python runtime. Source terminology: https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/04%20TA98%20EN.htm .

Integration passed: full model suite, all 7,474 hierarchy selectors/descendant sets across seven models, and Blender reopening with 6,330 objects and 38 packed textures. All six corrected source partitions retained exact exported coordinates. See `reports/ta98-male/muscle-groups-integration-validation.json`.


## Oral muscles: two partial groups, three individual selections, two surface corrections

The authoritative descendant inventory for tongue and soft-palate/fauces muscle groups, whole tongue and soft palate precedes geometry changes in `reports/ta98-male/oral-muscles-plan.json`. Five selections now use twelve existing textured sources: tongue muscle group A04.1.05.001 (bilateral genioglossus, hyoglossus and styloglossus), soft-palate/fauces group A04.1.06.001 (bilateral levator veli palatini, tensor veli palatini and palatopharyngeus), and individual styloglossus A05.1.04.105, levator A05.2.01.102 and tensor A05.2.01.103. These partial selections reuse geometry; no duplicated tissue was created.

Both palatopharyngeus sources contained a vertex joining two disjoint surface fans. Each was separated with a provisional 0.05 mm local displacement. No source faces were added or removed; all other source coordinates and original face ownership remain exact. Original image materials remain, with illustrative planar UVs on the two corrected meshes. Twelve envelopes pass closed edge, single-cycle vertex-link, duplicate/degenerate face and measured self-intersection checks. All forty measured sections close and were inspected, including cuts adjacent to both corrections. Complete group views, individual source views and correction closeups were inspected. Coarse source contours, retained nonzero-genus tunnels and cross-source attachments remain unverified; topology checks do not certify anatomy or absence of intermesh overlaps.

FJ2762 labels both whole uvula and uvular muscle. Its outer envelope is retained as a whole-uvula source but explicitly excluded from musculus uvulae A05.2.01.104 throughout hierarchy correspondence. Intrinsic tongue muscle volumes, independently established uvular muscle, palatoglossus, palatine aponeurosis and finer tissue boundaries remain missing. No outer tongue or palate envelope is relabeled as those tissues.

Evidence: `oral-muscles-plan.json`, `oral-muscles-repaired-review.json`, `oral-muscles-candidate-validation.json`, `oral-muscles-partition-validation.json` and `oral-muscles-validation.json` under `reports/ta98-male/`. Run `npm run test:ta98-oral-muscles` with the bundled NumPy Python runtime. This is local geometry and selection progress; complete TA98 mesh coverage remains unfinished.

Integration passed: full model suite, all 7,474 hierarchy selectors and descendant sets across seven models, and Blender reopening with 6,330 objects and 38 packed textures. Both corrected source coordinates survive export exactly. See `reports/ta98-male/oral-muscles-integration-validation.json`.


## Pharyngeal tissue sources: two raphes and salpingopharyngeus

The initial 44-term pharyngeal descendant inventory, source correspondences, tissue distinctions and adjoining interfaces precede work in `reports/ta98-male/pharyngeal-sources-plan.json`. Three reviewed selections reuse five named records: pharyngeal raphe A05.3.01.101 (FJ2749), bilateral pterygomandibular raphe A05.3.01.102 (FJ2744/FJ2756), and bilateral salpingopharyngeus A05.3.01.115 (FJ2745/FJ2757). Raphes are not relabeled muscular tissue; salpingopharyngeal mucosal folds are not represented by the muscle envelope. Independent constrictor subparts, fascia/wall layers and contacts remain missing or unverified.

FJ2749 originally contains three duplicated-face pairs, eleven nonmanifold edges, one measured nonadjacent self-intersection and five connected components. Duplicate removal/local correction and a subsequent exact Boolean both failed strict topology checks and were not installed. The installed provisional candidate is reconstructed directly from the original source at 0.05 mm voxel resolution, followed by 0.05 decimation ratio: 56,058 triangles, one closed genus-zero component. Original vertex-sampled distance to the candidate reaches 0.86581 mm; candidate vertex-sampled distance to the original reaches 0.01775 mm. These are bidirectional vertex samples, not a Hausdorff certificate. Original thin/internal fragments and subvoxel contacts can disappear or merge. Their anatomical identity is unresolved. Therefore this source is marked reconstructed-region-boundary-provisional, with no claim of complete anatomical repair.

The visible main contour remains source-shaped and rodlike; lower-end irregular cross-section contours remain. All five candidate envelopes pass closed-edge, single-cycle vertex-link, duplicate/degenerate face and measured self-intersection checks. All nineteen measured sections and the front/back views were inspected. Four source records retain their original geometry, UVs and image materials. The raphe retains its image material with illustrative planar UVs; this is not histological fiber evidence. No tissue attachment union or complete pharyngeal assembly is asserted.

Evidence: `pharyngeal-sources-plan.json`, `pharyngeal-sources-review.json`, `pharyngeal-sources-repaired-review.json`, `pharyngeal-sources-candidate-validation.json`, `pharyngeal-sources-partition-validation.json` and `pharyngeal-sources-validation.json` under `reports/ta98-male/`. Run `npm run test:ta98-pharyngeal-sources` with the bundled NumPy Python runtime. All local model redistribution restrictions remain in effect. Full TA98 geometry coverage remains unfinished.

Inventory correction: six tonsil-related terms cross-listed under the lymphoid hierarchy were added after initial geometry staging, bringing the whole-organ scope to fifty terms. This sequencing gap is recorded in the plan; no further geometry modifications preceded that supplement. Future preparation includes all A05.3 terms and their descendants.

Integration passed: full model suite, 7,474 hierarchy selectors across seven models, and Blender reopening with 6,330 objects and 38 packed textures. The repaired raphe coordinates survive export exactly. See `reports/ta98-male/pharyngeal-sources-integration-validation.json`.


## Duodenum: four textured regions with exact shared pyloric surface

All applicable duodenal descendants and common small-intestinal wall terms were recorded before geometry work in `reports/ta98-male/duodenum-parts-plan.json` (35 terms). Superior A05.6.02.002, descending A05.6.02.005, inferior A05.6.02.007 and ascending A05.6.02.008 regions now select original duodenal wall faces. Regional assignments follow observed source turns and provisional spatial zones: superior y>1.116 m; descending remaining x<-0.027 m/y>1.065 m; ascending remaining x>=0.004 m/y>1.066 m; inferior remaining lower curvature. These boundaries are geometric inferences, not certified flexure locations. Ampulla, papillae, gland volumes, suspensory muscle and independent histological layers remain unresolved. Flexures are not invented solid tissues.

Before partition, the complete shared gastric tract was inspected front/back and in measured sections. Its 412,470 triangles form one closed genus-zero component with valid single-cycle vertex links and no detected nonadjacent intersections. Four initial apparent cardia intersections were traced to adjacent triangles with duplicated attribute vertices at exact identical coordinates. The surface checker now recognizes geometric adjacency; a Blender regression retains a real nonadjacent crossing while excluding an attribute seam.

The original open duodenal ownership patch has two branched boundary vertices at the existing pyloric seam, despite the valid complete tract. Their one-triangle fans are retained as two supplementary superior-region patches. Six new mesh records cover four TA98 regions without changing, discarding or duplicating any of the original 94,070 faces. The source object is suppressed, and every concept/audit selection referencing it uses its six replacements. Every original coordinate, quantized normal, UV and material image is preserved; no new cap or wall is generated. All six patches have valid cycle/path vertex links and paired regional boundaries. The visible gastric junction validator uses the new patches, preserving the existing pyloric ring exactly.

Whole colored front/back/oblique views, individual regional views and all thirty measured patch sections were inspected. Open arcs terminate on regional/source interfaces; isolated regional edges remain visibly jagged along original source triangles. The original illustrative approximately 1 mm nested wall and distal closure remain provisional. No complete anatomical duodenum, jejunal continuity, histological layer or attachment claim follows from this partition.

Evidence: `duodenum-source-review.json`, `duodenum-source-vertex-links.json`, `duodenum-regions-review.json`, `duodenum-candidate-validation.json` and `duodenum-validation.json` under `reports/ta98-male/`. Run `npm run test:ta98-duodenum` with the bundled NumPy Python runtime. Full TA98 mesh coverage remains incomplete; all model assets remain local.

Integration passed: full model suite, 7,474 hierarchy selectors across seven models, and Blender reopening with 6,336 objects and 38 packed textures. All six patch coordinates survive export exactly. The visible tract has 264 exact pyloric shared edges with identical packed normals. See `reports/ta98-male/duodenum-integration-validation.json`.


## Gallbladder: four textured exterior regions with conforming shared seams

The 29-term biliary inventory, authoritative TA98 parentage, source matches, tissue distinctions and planned interfaces were recorded before geometry work in `reports/ta98-male/gallbladder-parts-plan.json`. Four exterior regions now have separate textured mesh selectors: fundus A05.8.02.002, infundibulum A05.8.02.003, body A05.8.02.004 and neck A05.8.02.005. These are provisional morphological surface zones, not histologically established tissue volumes. Independent wall layers, folds, glands, sphincters and an internal lumen remain unresolved; no exterior envelope is relabeled as those structures.

Whole-source review found one collinear triangle. Two neighboring triangles were retriangulated over exactly the same source surface before partitioning. All 363 original positions, packed normals and UVs remain bit exact; covered area is unchanged. The repaired 722-triangle exterior is one closed genus-zero component, with valid cycle vertex links, no duplicated or degenerate faces and no detected nonadjacent intersections. Its measured sections close. These checks establish the tested exterior topology, not a hollow organ lumen.

The first source-edge regional candidate had jagged boundaries and was superseded before installation. The installed regions use three shared planes along the source long axis, subdividing triangles and interpolating source UVs/normals. They contain 159 fundus, 489 body, 215 infundibulum and 95 neck triangles (958 total). Every original source face surface is covered once within float32 tolerance; maximum source-plane deviation is 5.38137e-8 m and maximum per-source-face accumulated area error is 1.70870e-10 square meters. Original vertex attributes remain exact. New boundary coordinates, packed normals and UVs agree exactly across patches. No overlapping caps are added. The original whole source is suppressed, and its concept/audit selections use the four replacements.

Whole colored front/back/oblique views, individual regional views and all twelve measured patch sections were inspected. Isolated regions expose open tissue interfaces; assembled exterior remains continuous. Source low-polygon contours and the capped neck persist. Nearest original source vertex distances are 7.30308 mm between gallbladder and cystic duct, 5.05915 mm between cystic and bile ducts, and 0.31757 mm between gallbladder and bile duct. These are vertex samples, not exact surface-distance certificates. Duct registration, capped tips and lumen continuity remain unresolved. This partition does not complete the biliary assembly.

Evidence: `gallbladder-source-repaired-review.json`, `gallbladder-source-validation.json`, `gallbladder-regions-review.json`, `gallbladder-candidate-validation.json` and `gallbladder-validation.json` under `reports/ta98-male/`. Run `npm run test:ta98-gallbladder` with the bundled NumPy Python runtime. Full TA98 mesh coverage remains incomplete. All review meshes remain local under the redistribution restriction.

Integration passed: full model suite, all 7,474 hierarchy selectors across seven models, and Blender reopening with 6,340 objects and 38 packed textures. Four regional patch coordinates and the corrected whole-source coordinates survive export exactly. See `reports/ta98-male/gallbladder-integration-validation.json`.


## Biliary trunk and tributaries: two continuous exterior regions, ten reviewed selectors

The full 29-term biliary inventory and named source matches precede geometry changes in `reports/ta98-male/biliary-sources-plan.json`. Seventeen sources were inspected separately, grouped by TA98 term, and together in colored front/back/oblique assembly views. All 72 measured source sections close. Sixteen sources comprise 23 closed genus-zero components; the alternate ZA:Bile duct is an open genus-zero tube with 24 terminal boundary edges. All sources pass duplicate/degenerate face and cycle/path vertex-link checks, with no detected nonadjacent self-intersections. Several hepatic tributary records contain disconnected capped components; their tissue identity/contact continuity remain unresolved. Closed topology does not establish internal duct lumens.

The source FJ3079, named “Common hepatic duct,” extends below the cystic contact while the separate ZA:Bile duct follows a displaced alternative route. Anatomical terminology places the common hepatic/bile division at the cystic junction (reference: https://pubmed.ncbi.nlm.nih.gov/29266455/). The new allocation is an explicit inference from source morphology/contact, not patient-specific anatomical certification. The nearest sampled cystic/trunk vertices are 0.63014 mm apart. A transverse source-Y plane at 1.1619912385940552 m partitions the existing FJ3079 exterior into proximal common hepatic A05.8.01.061 (232 triangles) and distal bile A05.8.02.013 (448 triangles). Those 680 triangles cover the 626 original face surfaces once within float32 tolerance, with shared seam coordinates, packed normals and UVs. All original source vertices and texture attributes remain exact. Maximum source-plane deviation is 1.37413e-9 m. No caps, internal wall or connection bridge were created.

Both candidates were inspected front/back, together from three views, and in all six measured sections. The original whole trunk and displaced alternative bile route remain preserved but suppressed. Existing concepts and audit selections use the visible replacement meshes, including proximal-only selection for the common hepatic concept. Eight further reviewed selections reuse original sources for right/left hepatic ducts, anterior/posterior and lateral/medial tributaries, left caudate duct and cystic duct. Their geometry and textures remain unchanged. Two trunk terms retain reconstructed-region-boundary-provisional status; the other eight are component-subset selections.

The visible assembly now has one trunk exterior, but capped duct ends, cystic-to-trunk attachment, gallbladder registration, disconnected tributary fragments, internal lumen continuity and duodenal termination remain unresolved. Right caudate and interlobular duct sources have not been identified. Independent mucosa, folds, glands, sphincters and ampulla are missing. No generic exterior is substituted for those tissues. Full TA98 mesh coverage remains incomplete; these assets stay local.

Evidence: `biliary-sources-review.json`, `biliary-source-topology.json`, `biliary-trunk-review.json`, `biliary-trunk-candidate-validation.json` and `biliary-trunk-validation.json` under `reports/ta98-male/`. Run `npm run test:ta98-biliary` with the bundled NumPy Python runtime.

Integration passed: full model suite, all 7,474 hierarchy selectors across seven models, and Blender reopening with 6,342 objects and 38 packed textures. Both new trunk patches retain exact exported coordinates. Three views of all 21 installed biliary/gallbladder records were inspected with their exported materials; existing visible branch/contact gaps remain unresolved. See `reports/ta98-male/biliary-integration-validation.json`.


## Bladder regional subdivision — 2026-10-01

Installed three textured regional patches: body (A08.3.01.003, 672 triangles), fundus (.004, 647), and neck (.005, 101). The original 506-face bladder envelope is retained as a suppressed source record; its surface is covered once by 1,420 conforming triangles. All original positions, packed normals and UVs are unchanged. New shared coordinates and texture attributes match exactly; rounding deviates from source face planes by at most 28.72 nm. Area conservation is exact before float32 rounding and remains inside its analytical rounding bound afterward.

The complete 72-term bladder/ureter/male urethral/attachment inventory is in `reports/ta98-male/bladder-parts-plan.json`. Source and regional front/back views, whole colored front/back/oblique/bottom views and measured sections were inspected. Both source inner and outer contours and the existing lower opening remain. Body/fundus have separate inner and outer surface components; those are not independently claimed histological layers. The source neck rim connects the whole envelope. Open regional edges are paired tissue interfaces, without duplicated faces or added caps.

Boundaries are explicitly provisional: the neck uses a localized rectangular spatial footprint around the source lower opening. Source wall thickness, detailed anatomy, and ureter/urethral contacts are not certified. Nearest source vertex samples are 1.533 mm left ureter, 1.615 mm right ureter and 0.189 mm urethra; these are not surface contact distances. Uvula, trigonal and detrusor muscle layers, mucosa/submucosa, crests, ligaments and microscopic detail remain unresolved. Openings/apex are not substituted with solids. This is additional provisional regional geometry, not a completed bladder reconstruction.

Local review only. Reproduce installed checks with `npm run test:ta98-bladder`; integration evidence is recorded in `reports/ta98-male/bladder-integration-validation.json`.


## Ureter abdominal/pelvic regions — 2026-10-01

Installed four textured surface patches for A08.2.01.002 and A08.2.01.003, bilaterally. Whole left/right ureters remain connected uncapped tubes with exactly their original two endpoint boundary loops. Shared regional seams reuse coordinates, packed normals and UVs; no source surface is removed, duplicated or capped. Original 6,048 triangles become 6,152 triangles through conforming clipping only.

Inventory `reports/ta98-male/ureter-parts-plan.json` covers all seven authoritative ureter terms. The Medical University of Białystok urinary anatomy lecture, pages 33–34, was downloaded locally and its diagrams inspected. It identifies the iliac crossing separating abdominal and pelvic courses. The original common-iliac distal boundary-ring centroids establish a model height of approximately 0.982879 m for each side; the planar boundary remains provisional. Whole context and colored source/region assemblies, individual regions and twelve transverse sections were inspected. Terminal renal/bladder tube shapes remain unchanged.

Intramural course, adventitia, muscle coat and mucosa remain missing; a single source tube does not establish these layers. Bladder-side crossings revealed by the packing checker remain unresolved and are not certified urinary connections. Source endpoints are intentional uncapped boundaries, not approved final anatomy. Region proof is in `ureter-candidate-validation.json`; installed attribute/preservation checks are in `ureter-validation.json`. Model and runtime remain local review assets.


## Male urethral regional subdivision — 2026-10-01

The complete 31-entry male urethral inventory, including all named regions and descendant structures, precedes geometry work in `reports/ta98-male/urethra-parts-plan.json`. The Medical University of Białystok sixth abdominal lecture, page 61, was visually inspected; its diagram distinguishes preprostatic, prostatic, intermediate and spongy segments. Original source urethra/prostate/spongiosum surfaces and all nine measured source sections were inspected. The source urethra is a single uncapped surface with two original twelve-edge endpoint loops, not a separately verified wall/lumen model.

The source tessellation supplies 121 twelve-point rings; the script verifies every face connects consecutive rings. Centerline points cross the original prostate envelope between rings 92/93 and the spongiosum envelope between 82/83. Direct ray/triangle hits establish points at (approximately) y=0.811633,z=-0.00809660 m and y=0.807903,z=0.00656851 m. Shared Z-plane regional boundaries pass through these measured source landmarks but remain anatomically provisional over the complete tube circumference.

Installed A09.4.02.004 prostatic urethra (692 triangles), A09.4.02.017 intermediate part (292), and A09.4.02.021 spongy urethra (2,000). All 2,880 original source face regions are covered exactly once by 2,984 conforming triangles. Original vertices/normals/UVs and both endpoint edges are unchanged; new seam coordinates/normals/UVs are shared exactly. Three regional front/back views, nine sections and the colored whole assembly were inspected. The male whole-urethra selector now assembles the three regions, with the incomplete anatomy caveat retained.

No new intramural/preprostatic segment, prostatic proximal/distal boundary, urethral crest, seminal colliculus, utricle, glands/ducts, sphincter or wall layers are asserted. The source upper end intersects the neck/prostatic context; the intramural transition is unresolved. Neither the prostate nor the spongiosum has an independently established urethral tissue/lumen boundary. Packing and bladder/prostatic/urethral contact crossings remain unresolved; open endpoints are not certified urinary junctions. See `urethra-candidate-validation.json`, `urethra-validation.json`, `urethra-landmarks.json`, and `urethra-integration-validation.json`.

## Oesophagus: whole-wall repair before regional construction

The complete 17-term oesophageal inventory is recorded in `reports/ta98-male/oesophagus-parts-plan.json`, including wall layers, glands, constrictions and supporting muscle/tendon terms. No new regional coverage is claimed yet.

Five rendered landmark views and eleven measured transverse sections of the existing source were inspected. The current reconstructed oesophagus is free of measured self-intersections and shares 397 exact cardia boundary edges, but it has 212 triangle crossings with the cricoid cartilage and 733 with the central tendon. The original native oesophagus also has tissue crossings, so simply replacing it with the native mesh is insufficient. Crossing counts are triangle pair counts, not anatomical overlap volume or severity.

A local proximal repair candidate is constructed at the measured inferior cricoid border (Y = 1.4568747474 m). It retains the cardia boundary exactly and the source per-corner UV layer, removes the upper extension, and reconstructs only an annular wall rim. The lumen disk remains open. The candidate has 98,853 triangles, zero measured self-intersections, zero cricoid intersections, valid vertex links, no duplicate/degenerate/nonmanifold faces, one connected component and the two original cardia boundary loops. Three rendered proximal views and three transverse sections showing two closed wall contours were inspected.

Candidate files remain local: `.local-models/ta98-review/oesophagus-proximal-candidate.blend` and `.json`. This candidate is **not installed**: the oesophagus/diaphragm tissue crossings, correct right-crural hiatus, pharyngeal tissue connection, gastric assembly proof and regional/histological representations remain unresolved. Nearest centerline approach to the central tendon was rejected as a hiatus landmark; its minimum fell at the sampled endpoint. The right-crural nearest-point minimum alone likewise cannot establish the anatomical boundary. Source landmark measurements and this rejection are retained for reproducibility.

Scripts: `prepare-ta98-oesophagus.py`, `review-ta98-oesophagus-landmarks.py`, `measure-ta98-oesophagus-landmarks.py`, `check-ta98-oesophagus-packing.py`, `compare-ta98-native-oesophagus.py`, `repair-ta98-oesophagus-proximal.py`, `render-ta98-oesophagus-proximal.py`. Reports include source/candidate hashes, exact interface checks, topology, crossing examples and explicitly incomplete anatomy. Existing installed geometry and the current coverage numerator remain unchanged.

### Coupled passage, generic wall and crural repair candidates

The complete source diaphragm is a closed genus-two assembly. An oblique ray projection at slope 0.5 locates two enclosed paths; the small patient-left path is used provisionally as the oesophageal passage, not as a certified anatomical rim. The source cardiac tube's measured section centroids guide a common field applied to every oesophageal, gastric and duodenal region. This retains face ownership and per-corner UVs and fixes all shared coordinates consistently. It clears the central tendon but originally compressed the generic wall to about 0.45 mm and retained 472 right-crural and 406 left-crural triangle crossing pairs.

`restore-ta98-passage-wall.py` now produces a locally verified wall candidate. A nearest-point offset collapsed two triangles; a subsequent unrestricted normal move folded the cardia seam. Both were rejected. The successful correction moves 6,466 canonical inner vertices along interpolated exterior normals projected into XZ, preserving Y and the fitted exterior exactly. Maximum additional displacement is 1.114 mm. The whole 408,500-triangle assembly has one closed genus-zero component, valid vertex links, no open/nonmanifold edges, duplicate/degenerate faces, winding defects or detected self-intersections. `validate-ta98-restored-wall.py` independently confirms shared coordinates, unchanged exterior/Y, face indices, corner UVs and materials. Every changed inner vertex measures approximately 1.0000–1.0002 mm from the exterior. Six inspected transverse sections have 0.9405–1.0014 mm minimum 2D contour gaps. These are generic wall measurements, not histological layers or a continuous minimum-thickness proof.

Candidate: `.local-models/ta98-review/oesophagus-wall-thickness-candidate.json` and `.blend`. Four textured whole/context views and the six sections were inspected. Review rendering computes common normals from the entire assembly, preventing artificial UV-split vertex shading seams. Visible source surface stepping remains, particularly on the stomach and duodenum; corrected normals alone do not repair it.

`measure-ta98-gastric-penetration.py` tested 162,964 gastric vertices against the complete diaphragm. Before crural fitting, 2,950 vertices were inside its tissue, reaching 6.179 mm from its nearest boundary. The short exits were positive X and anterior Z; a downward ray can run along a crus for over 100 mm and is unsuitable as a simple repair direction.

The archived common 6.5 mm X / 6.5 mm Z field candidate preserves the oesophageal seam at Y≥1.223, shared coordinates, regional Y coordinates, faces and corner UVs. It has no detected self-intersections or crossings with the complete diaphragm, cricoid or T10. Vertex containment finds zero confirmed inside points. All 695 ambiguous parity-ray cases were independently resolved as outside by float64 solid-angle winding. Six inspected section gaps are 1.0002–1.1104 mm. Source-based field directions and boundaries are provisional; no anatomical certification follows.

Archived candidate: `.local-models/ta98-review/gastric-crural-fit-6.5-candidate.json` and `.blend`; reports `gastric-crural-fit-6.5-{validation,sections,containment,containment-winding,neighbors}.json`. Four context views and six sections were inspected. The junction-side context view is largely occluded by the diaphragm and cannot prove the junction's shape.

**This candidate is not accepted or installed.** The neighbor comparison finds liver surface crossing pairs increasing from 1,854 to 1,919. Pair counts are not overlap volume, but the changed contact requires correction and review. Other tested crossings are inherited and unchanged: pancreas 615, right kidney 553, transverse colon 1,136, ascending colon 142, jejunum 163, primary heart source parts 2,546, and inferior right lung source parts 460. Tested spleen, left kidney, suprarenals, descending colon and other lung-lobe assemblies have no detected surface crossings. This does not establish complete tissue volume disjointness or cover vessels, nerves, ligaments and missing neighboring anatomy.

Pure-X and reduced-anterior alternatives retain right-crural crossings. Tighter neck support can fold the discretely triangulated wall even when the analytical field has positive Jacobian; these attempts are rejected by the unchanged self-intersection guard. Diagnostic reports retain parameters and crossing coordinates. Further work must jointly resolve liver/crural packing, remove visible surface stepping, establish all adjoining tissue interfaces and replace provisional anatomical boundaries. The inherited distal duodenal closure does not prove an intestinal lumen connection.

Additional scripts: `locate-ta98-diaphragm-passages.py`, `measure-ta98-cardia-tube.py`, `fit-ta98-oesophagus-passage.py`, `restore-ta98-passage-wall.py`, `validate-ta98-restored-wall.py`, `validate-ta98-passage-fit-sections.py`, `render-ta98-passage-fit.py`, `measure-ta98-gastric-penetration.py`, `fit-ta98-gastric-crura.py`, `resolve-ta98-diaphragm-containment.py`, and `check-ta98-gastric-neighbors.py`. No new TA98 regional coverage or main-model installation is claimed in this batch.

Reference text consulted: [NCI SEER stomach anatomy](https://training.seer.cancer.gov/ugi/anatomy/stomach.html), describing the stomach below the diaphragm and predominantly left of midline, with cardia surrounding the oesophageal junction and fundus above/left of that opening. This establishes general relationships, not model-specific dimensions. The diagram link was opened; its pixels were not inspected, so diagram-based shape validation remains unproven.

### Bounded whole-wall surface regularization

`smooth-ta98-gastric-wall.py` stages a smoother version of the archived 6.5 mm crural-fit candidate. It uses exact-coordinate whole-assembly adjacency and twelve bounded Taubin cycles. Oesophageal vertices, the existing distal closure, near-diaphragm vertices and a two-ring collar around all existing liver-crossing triangles are fixed. An initial pass produced two local folds and was rejected. The successful pass additionally preserves the four source triangles implicated in those folds and their two-ring neighborhoods. These protected patches remain unregularized; their anatomical feature classification is unestablished. Their source edges range from 0.508 to 6.565 mm, exposing uneven tessellation that requires a subsequent geometry-quality repair, not merely more aggressive smoothing.

Candidate: `.local-models/ta98-review/gastric-smoothed-wall-candidate.json` and `.blend`. It moves 128,114 canonical vertices by at most 0.350060 mm including float32 storage rounding, adds no triangles, and retains the 408,500 source triangles and per-corner UVs. `validate-ta98-gastric-smoothing-fidelity.py` independently proves original indices/materials/source ownership, every shared coordinate, unchanged oesophagus, protected patch geometry and distal closure. Its analytical storage rounding bound is under 0.060 micrometres. Whole topology remains one closed genus-zero component without edge/link, winding, duplicate or degenerate-face defects; no self or complete-diaphragm crossings are detected. The 1,919 existing liver-crossing pairs remain; this is not liver packing acceptance.

Four textured context views were inspected. Surface stepping is visibly reduced but persists in protected regions and parts of the underside/duodenum. The side view remains occluded by the diaphragm. The six upper junction sections remain closed with approximately 1.0002–1.1104 mm minimum 2D gaps. `review-ta98-gastric-wall-sections.py` additionally compares twenty complete assembly sections from Y=1.10817 to 1.22717 m, showing every closed contour and its source owners, nesting and minimum 2D separation. All twenty contour graphs were inspected and remain closed with unchanged contour counts. Candidate minimum nested-contour gaps range from approximately 0.6017 to 1.0335 mm in those sections; source minimum gaps range from 0.6166 to 1.0344 mm. These section diagnostics expose variable generic body-wall thickness; they do not establish a uniform 1 mm wall, clinical thickness, histology or full lumen/junction acceptance.

Reports: `gastric-wall-smoothing-validation.json`, `gastric-wall-smoothing-fidelity.json`, `gastric-smoothed-wall-sections.json`, `gastric-smoothed-body-sections.json`, and `gastric-smoothing-protected-source-patches.json`. Rejected-pass evidence remains in `gastric-wall-smoothing-trial-12-validation.json`. The main Male TA98 catalogue is unchanged, no additional terms are counted as completed, and this local candidate remains uninstalled pending wall/triangulation quality and joint liver/crural packing corrections.


## On-demand runtime integration

The lossless local runtime is now previewable with `model=local-ta98&delivery=1`. The viewer loads only selection/visible-layer/section chunks, with three concurrent requests, cancellation of obsolete work and reuse of completed geometry. Hidden review geometry is packed separately. The isolated textured stomach was visually inspected in dark theme: four required chunks, four loaded, zero pending. Full attribute byte comparison and gzip round trips passed. No anatomy was changed, no TA98 term added, no review asset deployed. Detailed delivery budgets are in `reports/ta98-male/runtime-on-demand.json`. Whole-body rendering, frame-time measurements, progressive overview LOD and memory eviction remain open.


## Detailed male placement constraints

The required construction order is now explicit: primary detailed-male skeleton/fascia → valid compartment/available-space boundaries → adjacent structures → whole-organ fit/repair → full TA98 partition. `prepare-ta98-placement-reference.py` stages original geometry with source offsets applied once and independently compares exact world vertices/triangles against the master. Gastric reference set: 633 native source guide parts inventoried, 372 staged nearby, all 372 frame-identical, 641,472 source triangles. No skeleton/fascia was deformed and no candidate placement has been accepted. `check-ta98-placement-reference.py` and placement overlay rendering provide the next geometric constraints. Open sheets, ribs with intentional gaps and missing compartment walls cannot establish an exact closed available volume by themselves.

## Native oesophageal profile repair

The skeletal/fascial frame exposed an oversized circular oesophageal reconstruction. A staged radial fit now follows the native detailed-male outer wall/lumen contours: 42 section samples, 43,273 moved vertices, max 4.817 mm. The final thoracic-end variant preserves all faces/UVs/materials and lower junction coordinates, passes whole topology/self-crossing checks, and clears all tested bones. The intermediate thoracic variant that crossed C6 was rejected. Adjacent heart/airway/aortic/pulmonary contacts remain; the candidate is unaccepted and uninstalled. Full evidence and reproducible scripts are recorded in `TA98-OESOPHAGUS-PLACEMENT.md`. No new distinct TA98 terms are counted.

## Adjacent oesophageal/atrial packing continuation

The first bounded course trial passed topology/self-intersection and all skeletal-guide gates but was rejected after actual sections exposed passage into the atrial chamber. The fitter now includes measured posterior-atrium sidedness and refuses to stage the failed +/-10 mm search. Four sampled levels have no independently bone-clear/posterior-centroid grid option. The primary atrial source and its original UVs were staged for subsequent whole-source/interface repair; all 113 heart terms and five atrial descendants are recorded. Two open source boundary cycles remain unidentified for anatomical purposes. No original source was deformed, no failed trial installed, and no copied mesh counted as new term coverage. See `TA98-OESOPHAGUS-PLACEMENT.md` and `reports/ta98-male/left-atrium-packing-plan.json`.


### Atrial local triangulation repair

Removed three near-collinear source triangles using internal diagonal flips, preserving all 2,607 original vertices, UVs, materials and opening coordinates. The candidate retains 5,105 faces and independently passes self-crossing/topology checks; matching textured views were inspected. No additional TA98 term is counted, and no installed model change is claimed. Atrial wall interfaces and oesophageal/cardiac packing remain unresolved. See `TA98-OESOPHAGUS-PLACEMENT.md` and the candidate-hashed atrial sliver reports. Rejected posterior-column deformation remains uninstalled.

## Patellar source surface dissection

Male TA98 now includes independently selectable bilateral base, apex, articular surface and anterior surface (`A02.5.05.002`–`.005`). These are four **provisional surface-feature representations**, not four newly constructed solid tissues. The five-term whole patellar inventory and adjoining interfaces were recorded before partitioning. Every original bone/cartilage triangle is assigned exactly once; remaining border faces retain explicit unclassified ownership. Original coordinates, packed normals, UVs and image materials are bit-exact, and the active triangle count is unchanged.

The original and partitioned assemblies each retain one closed genus-zero component per side, 276 triangles each, with no topology defects or detected self-crossings. Textured anterior/posterior/side views and 24 actual source cross-sections were inspected. Regional boundaries remain provisional; unchanged source placement is not a packing certificate. No wall layer or microscopic term is inferred from an exterior patch.

Patellar installation catalogue SHA-256: `f6c27f984a9bf6f4a17053f05b254e36756437f87294d2f214bfa69455b22b43`. The export validator passes exact face/attribute/material preservation; the rebuilt lossless runtime passes byte and gzip round-trip checks; all 7,474 hierarchy selectors pass for all eight models. Evidence: `reports/ta98-male/patella-{parts-plan,source-review,source-sections,regions-review,subdivision,export-validation}.json`. Candidate construction/preservation reports describe the pre-install stage. These four surface features are reported separately rather than added to a material-term completion percentage. Full anatomical acceptance remains incomplete.

## Hyoid regional dissection

Male TA98 now includes the body (`A02.1.16.002`), lesser horn (`.003`) and greater horn (`.004`) as three provisional material regions, with bilateral horns giving five textured mesh parts. The entire four-term hyoid inventory and body/horn, ligament, membrane and muscle interfaces were recorded first. Actual source anterior/posterior/superior/side views were compared with labelled Gray plate 186, and eighteen real axial/sagittal/coronal triangle sections were inspected before cutting. The original hyoid has 1,324 triangles and one closed genus-zero component with no topology defects or detected self-crossings.

The first facet-aligned partition was rejected for coarse jagged regional boundaries. Conforming cuts use one canonical edge/plane intersection cache and shared coordinate, normal and UV arrays. Provisional planar junctions are derived from the native protrusions and body shoulder; numerical shifts of at most 30 micrometres avoid sliver cells without moving source vertices. These boundaries are not histological or subject-specific anatomical acceptance. The refined whole assembly has 944 vertices and 1,884 triangles, one closed genus-zero component, no topology defects and no detected self-crossings at working or runtime float32 precision. No artificial caps or duplicated exterior tissue were added. Independent per-source-triangle barycentric checks prove full coverage without positive-area overlap (maximum overlap `1.11e-16` in source barycentric coordinates), original surface distance below `5.59e-17` m, original/interpolated UVs and normals, and unchanged source vertices. Actual four-view textured candidate renders were inspected.

Installed catalogue SHA-256: `ef0aca2a70c92ea8a219509ba5599271d6e044a88f0cdc8f2479d05fbb19e4dd`. Exported buffers match the reviewed candidate at their packed precision; unrelated parts, materials and hidden-state metadata remain identical. The cuts add 560 triangles and approximately 8.3 KB compressed geometry. Rebuilt runtime byte checks and all eight models' 7,474 hierarchy selectors pass. Regional exterior patches share real seams but do not independently reconstruct complete internal region volumes; native whole-body packing remains unaccepted.

Evidence: `reports/ta98-male/hyoid-{parts-plan,source-review,source-sections,regions-conforming-construction,regions-conforming-preservation,regions-conforming-review,runtime-surface,export-validation}.json`; local Blender review: `outputs/ta98-hyoid-regions-conforming.blend`. Scripts with the same hyoid names preserve each construction and verification step. `update-ta98-created-coverage.py` now checks every declared constructed/dissected mesh against active catalogue records and authoritative entity types. Current provisional installed material coverage: **44 distinct terms / 288 mesh records, 0.7374% of the existing 5,967 eligible-term classification**; separately **4 immaterial surface-feature terms / 8 mesh records**. This is geometry presence, not verified anatomical completion. Copied source anatomy and correspondence-only aliases remain excluded.

## Bilateral fibular regions and corrected malleolus correspondence

The full 17-term fibular descendant inventory is recorded in `reports/ta98-male/fibula-parts-plan.json`. It distinguishes material bone regions, surface concavities/articular facets, and boundary curves. Native bilateral bone/cartilage, adjacent tibia/talus, relevant ligaments/nerve and the original malleolar skin landmarks were staged in one measured coordinate frame. Three whole source views and 24 actual triangle cross-sections were inspected against the labelled Gray plate 258. The first side-camera setup cropped the long bone; its orientation was corrected and the complete side view rechecked. Multiple closed sagittal contours reflect the curved source bone and are retained rather than mistaken for holes.

Four provisional material regions are now installed bilaterally: head (`A02.5.07.002`), neck (`.005`), shaft (`.006`) and lateral malleolus (`.014`). They comprise 12 bone/cartilage mesh records. The reusable `ta98_source_partition.py` makes conforming cuts with a common coordinate-edge cache across material seams, preserves all original vertices and interpolates original normals/UVs. Per-source-face barycentric checks prove complete coverage without positive-area overlap, area error below `1.67e-19` square metres, UV error below `4.93e-16`, and surface distance below `6.36e-17` metres. Each whole side retains one closed genus-zero component, with zero topology defects or detected self-crossings at working and runtime precision. No false caps are added; internal regional tissue volumes and anatomical boundary acceptance remain unresolved. Native packing is not certified by preserving source placement.

The original `.014` correspondence incorrectly selected skin objects named “lateral malleolus.” Those objects remain skin with unchanged geometry, while explicit source exclusions remove the false bone correspondence. The local TA98 hierarchy now selects the dissected bony malleolus. Two articular-facet selections (`.003`, `.015`) follow the actual original proximal/distal cartilage footprints; these are recovered correspondences and are **not counted as newly constructed geometry**. Ten additional fibular feature terms remain explicit: `.004`, `.007`–`.013`, `.016`, `.017`. A regional cut does not pretend to establish the full inventory's crests, borders, fossa or groove.

Installed catalogue SHA-256: `57e460fb47272f8d3b6dedc8af5fd8c37e3a96c1fc12c73e6b80b30791c3d108`. Actual packed buffers match the reviewed candidate and surrounding part/material metadata is unchanged. The cuts add 400 triangles total and about 6.4 KB compressed geometry. Lossless runtime verification and all eight models' 7,474 hierarchy selectors pass after rebuilding; the false skin-to-bone mapping is independently checked absent from the local bone selector.

Current provisional created/dissected material coverage: **48 distinct terms / 300 mesh records, 0.8044% of the existing 5,967 eligible-term classification**. Separately: **4 constructed surface-feature terms / 8 mesh records**. Copied source meshes, recovered correspondence aliases and repairs of existing whole anatomy remain excluded. Fully verified anatomical completion remains unproven. Evidence: `fibula-{parts-plan,source-review,source-sections,regions-construction,regions-preservation,regions-review,runtime-surface,export-validation}.json`; Blender artifact `outputs/ta98-fibula-regions.blend`.

## Native cardiac valve context (2026-10-03)

Named cardiac relationships now include native valve regions, adjoining chamber envelopes and available outflow vessels. Nine male valve/leaflet records and all four native female valves resolve their available context. Four native male papillary muscles also resolve their ventricular chamber. Reciprocal adjacency adds 44 directed links in the male catalogue and 16 in the female catalogue. The female source's combined ventricular envelope is explicitly described as combined; it does not establish separate left/right ventricular geometry. Independently registered HRA donors, suppressed parts, opposite chambers and nervous-system ventricles are excluded.

The references are NIH NHLBI's [valve and blood-flow description](https://www.nhlbi.nih.gov/health/heart/blood-flow) and [OpenStax heart anatomy](https://openstax.org/books/anatomy-and-physiology-2e/pages/19-1-heart-anatomy). Whole-valve descriptions distinguish leaflets from the complete valve. Relationships are anatomical context, not certified annular contact, source wall attachment or lumen continuity. No geometry, texture, source placement, triangulation or anatomical acceptance status changed. The graph is computed once per catalogue and cached; it adds no meshes or textures.

`npm run test:cardiac-relations` checks both actual catalogues, reciprocal context, precise descriptions, exclusion of donors/wrong frames/suppressed parts/wrong chambers, and absence of a semilunar papillary apparatus. Full reciprocal checks pass for 6,273 active male parts / 16,088 links and 2,774 female parts / 3,996 links. A live female valve relationship-expansion check preserves position, target, up, quaternion and projection offset exactly within a stable viewport. Evidence: `reports/ta98-completeness/cardiac-relations.json`, `cardiac-camera-review.json`, and `outputs/female-valve-native-relations.jpg`. This closes a relationship gap, not the remaining full-body geometry coverage gap.
