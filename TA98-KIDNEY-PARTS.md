# Kidneys and collecting systems: priority 9

Whole-assembly plan and complete current TA98 source inventory before renal geometry changes.

## Completion work: one canonical kidney, then transfer

Use the installed male-left assembly as the master. The new registry inventory
at `.local-models/geometry-workbench/renal-completion/male-authority-inventory.json`
binds all 61 descendants to the local authoritative terminology; the two pelvic
variant classifications remain recorded separately below. Preserve every current
source and all skeleton geometry. Complete and inspect the master before transfer.

The first repair closes the four disconnected collecting perforations into one
coherent hilar capsule termination; the hilum remains an opening. Reconstruct the
full reference-guided cortical union and renal columns against the same native
inner capsule, eleven medullary units and complete collecting envelope. Inspect
the remaining sinus compartment rather than relabeling an arbitrary leftover as
cortex or fat. All interfaces must be reciprocal in one encoded bank.

Then subdivide the real medullary and collecting tissues into their named regions
and coats, retain open papillary/urine passages, and resolve vascular courses outside
the collecting lumen. Gross landmarks are selectable boundaries/surface regions;
microscopic organization needs explicit efficient detail geometry. Use Columbia's
[renal histology](https://histologylab.ctl.columbia.edu/lab12/kidney/) and the
[University of New England histology collection](https://digitalhistology.org/organs-systems/urinary/tubules-and-collecting-system/tubules-and-collecting-system-18/)
for the cortex/ray/labyrinth and medullary regional organization. These describe
structure; they do not establish patient-specific regional measurements.

Transfer the reviewed complete assembly with the recorded recipient size/pose and
free-ureter course. Recheck each recipient's fixed skeleton, suprarenal, vascular,
collecting and surrounding-organ interfaces; publish locally through the guarded
installer only after actual encoded surface/section review. No public deployment.

## Current local renal installation

The collecting-layer revision is installed in **both ordinary local TA98 models**, with one assembly per side. The male-left donor was finished and inspected first, then transferred at the recorded recipient sizes and placements. Current evidence is under `.local-models/geometry-workbench/renal-completion/layer-delivery/`: `male-installed.json`, `female-installed-v2.json`, both reviewed geometry jobs, and `served-verification.json`. The earlier construction remains recoverable through the installation rollback manifests.

| TA98 | Current installed representation in each kidney | Remaining work |
|---|---|---|
| A08.1.01.013 | One finite fibrous capsule with a single rounded hilar termination | Reciprocal cortical charts |
| A08.1.01.015 | Finite partial peripheral cortical field; numerical slivers cleaned within a 0.1 micrometre source-vertex simplification tolerance | Complete cortical extent, columns and fine organization |
| A08.1.01.020 | Gross medullary assembly selecting all eleven pyramidal bulk cells and eleven papillae | Outer/inner zones and finer subdivisions |
| A08.1.01.027 | Eleven finite papillary cells | Duct openings and microscopic surface organization |
| A08.1.01.029 | Eleven finite pyramidal bulk cells, with separate papillary controls | Individual whole-pyramid grouping and finer anatomy |
| A08.1.05.001 | Three finite calyceal/pelvic coats and three communicating UPJ coats | Individual calyceal regional controls and complete vascular fit |
| A08.1.05.009–011 | Separately selectable adventitia, muscular layer and mucosa of the calyceal/pelvic and UPJ walls | Patient-specific histological thickness and finer organization |
| A08.2.01.001 | Three continuous finite ureter coats along the retained recipient course | Abdominal/pelvic regional controls and intramural bladder junction |
| A08.2.01.005–007 | Separately selectable ureteric adventitia, muscular layer and mucosa | Finer histological organization |

Each side has **33 finite selectable pieces**, sharing one packed position/normal bank: capsule, peripheral cortex, eleven pyramidal bulk cells, eleven papillae and nine collecting-coat cells. Each coat has finite thickness, reciprocal shared interfaces and real annular regional joins. Computational urine disks are not emitted. Three-step coat peeling and individual TA98 selection are enabled. The older kidneys and unpartitioned walls remain suppressed and recoverable. Central leftovers and coarse source partitions are not relabelled as completed renal columns or microscopic tissue.

The female left/right uniform scales relative to the male-left donor are **0.735639 / 0.724906**. The male-right donor scale is **0.985672**, with a further common 6 mm inferior correction returning smoothly to the retained distal ureter course. Skeleton, suprarenals, native blood vessels and all unrelated anatomy remain literal. Fixed bone/suprarenal/whole-liver surface checks introduce no newly crossed previously disjoint neighbour; this is not full packing acceptance.

Actual donor exterior and finite sections, followed by recipient finite sections, were inspected in the existing Blender session. All **132 encoded tissue cells** pass closed topology, expected component/genus and exhaustive self-crossing checks; both reviewed packing jobs also pass occupied-volume separation. The female-left and male-right free-ureter coats are regenerated inside the transferred physical carriers, retaining both source annuli, rather than independently deforming interior faces. Two internal male-right controls required a bounded Float32 rounding correction of at most **0.466 nm**; external boundaries and literal source controls stayed fixed. The local server returns the installed catalogues and all four exact reviewed banks, with matching gzip payloads. No new UV arrays or images are used. Camera-intent, explosion-plane, explosion-viewport, hierarchy/coats and TypeScript checks pass. These are construction findings, **not complete renal anatomical acceptance**: full cortical/vascular organization, individual calyceal controls, microscopic terms, intramural bladder interfaces and full packing acceptance remain unfinished. Browser/GPU performance is not measured.

## Source inventory before reconstruction

The complete parent/classification/source table below records the pre-reconstruction inventory. Its “Unbuilt” cells describe that historical source snapshot; the installed additions and outstanding distinctions are listed above.

| TA98 | Term | Parent | Representation | Male sources | Female sources |
|---|---|---|---|---|---|
| A08.1.01.001 | kidney | A08.0.00.000 | whole organ | `ZA:Kidney.l`<br>`ZA:Kidney.r` | `LOCAL:female:l_kidney:189`<br>`LOCAL:female:r_kidney:192` |
| A08.1.01.002 | lateral border | A08.1.01.001 | boundary curve | Unbuilt | Unbuilt |
| A08.1.01.003 | medial border | A08.1.01.001 | boundary curve | Unbuilt | Unbuilt |
| A08.1.01.004 | hilum of kidney | A08.1.01.003 | passage or cavity landmark, not a solid plug | Unbuilt | `HRA:VH_F_hilum_of_kidney_L` (hidden reference)<br>`HRA:VH_F_hilum_of_kidney_R` (hidden reference) |
| A08.1.01.005 | renal sinus | A08.1.01.003 | passage or cavity landmark, not a solid plug | Unbuilt | Unbuilt |
| A08.1.01.006 | anterior surface | A08.1.01.001 | surface or positional landmark | Unbuilt | Unbuilt |
| A08.1.01.007 | posterior surface | A08.1.01.001 | surface or positional landmark | Unbuilt | Unbuilt |
| A08.1.01.008 | superior pole | A08.1.01.001 | surface or positional landmark | Unbuilt | Unbuilt |
| A08.1.01.009 | inferior pole | A08.1.01.001 | surface or positional landmark | Unbuilt | Unbuilt |
| A08.1.01.010 | renal fascia | A08.1.01.001 | distinct tissue envelope or wall layer; outer shell alone insufficient | Unbuilt | Unbuilt |
| A08.1.01.011 | paranephric fat | A08.1.01.010 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.012 | perinephric fat | A08.1.01.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.013 | fibrous capsule | A08.1.01.001 | distinct tissue envelope or wall layer; outer shell alone insufficient | Unbuilt | `HRA:VH_F_kidney_capsule_L` (hidden reference)<br>`HRA:VH_F_kidney_capsule_R` (hidden reference) |
| A08.1.01.014 | kidney lobes | A08.1.01.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.015 | renal cortex | A08.1.01.001 | tissue region or regional assembly | Unbuilt | `HRA:VH_F_outer_cortex_of_kidney_L` (hidden reference)<br>`HRA:VH_F_outer_cortex_of_kidney_R` (hidden reference) |
| A08.1.01.016 | cortical labyrinth | A08.1.01.015 | internal fine tissue organization; generic surface insufficient | Unbuilt | Unbuilt |
| A08.1.01.017 | cortex corticis | A08.1.01.016 | internal fine tissue organization; generic surface insufficient | Unbuilt | Unbuilt |
| A08.1.01.018 | medullary rays | A08.1.01.015 | internal fine tissue organization; generic surface insufficient | Unbuilt | Unbuilt |
| A08.1.01.019 | renal columns | A08.1.01.015 | tissue region or regional assembly | Unbuilt | `HRA:VH_F_renal_column_L` (hidden reference)<br>`HRA:VH_F_renal_column_R` (hidden reference) |
| A08.1.01.020 | renal medulla | A08.1.01.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.021 | outer zone | A08.1.01.020 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.022 | outer stripe | A08.1.01.021 | internal fine tissue organization; generic surface insufficient | Unbuilt | Unbuilt |
| A08.1.01.023 | inner stripe | A08.1.01.021 | internal fine tissue organization; generic surface insufficient | Unbuilt | Unbuilt |
| A08.1.01.024 | vascular bundles | A08.1.01.023 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.01.025 | interbundle region | A08.1.01.023 | internal fine tissue organization; generic surface insufficient | Unbuilt | Unbuilt |
| A08.1.01.026 | inner zone | A08.1.01.020 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.027 | renal papilla | A08.1.01.026 | tissue region or regional assembly | Unbuilt | `HRA:VH_F_renal_papilla_L_a` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_b` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_c` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_d` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_e` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_f` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_g` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_h` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_i` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_j` (hidden reference)<br>`HRA:VH_F_renal_papilla_L_k` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_a` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_b` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_c` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_d` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_e` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_f` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_g` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_h` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_i` (hidden reference)<br>`HRA:VH_F_renal_papilla_R_j` (hidden reference) |
| A08.1.01.028 | renal crest | A08.1.01.026 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.01.029 | renal pyramids | A08.1.01.001 | tissue region or regional assembly | Unbuilt | `HRA:VH_F_renal_pyramid_L_a` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_b` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_c` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_d` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_e` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_f` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_g` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_h` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_i` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_j` (hidden reference)<br>`HRA:VH_F_renal_pyramid_L_k` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_a` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_b` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_c` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_d` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_e` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_f` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_g` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_h` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_i` (hidden reference)<br>`HRA:VH_F_renal_pyramid_R_j` (hidden reference) |
| A08.1.01.030 | cribriform area | A08.1.01.027 | surface or positional landmark | Unbuilt | Unbuilt |
| A08.1.01.031 | openings of papillary ducts | A08.1.01.030 | passage or cavity landmark, not a solid plug | Unbuilt | Unbuilt |
| A08.1.02.001 | renal segments | A08.1.01.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.02.002 | superior segment | A08.1.02.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.02.003 | anterior superior segment | A08.1.02.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.02.004 | anterior inferior segment | A08.1.02.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.02.005 | inferior segment | A08.1.02.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.02.006 | posterior segment | A08.1.02.001 | tissue region or regional assembly | Unbuilt | Unbuilt |
| A08.1.03.001 | intrarenal arteries | A08.1.01.001 | vascular tissue network; lumen and junction review required | `ZA:Intrarenal arteries of left kidney`<br>`ZA:Intrarenal arteries of right kidney` | Unbuilt |
| A08.1.03.002 | interlobar arteries | A08.1.03.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.003 | arcuate arteries | A08.1.03.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.004 | cortical radiate arteries | A08.1.03.003 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.005 | afferent glomerular arteriole | A08.1.01.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.006 | efferent glomerular arteriole | A08.1.01.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.007 | perforating radiate arteries | A08.1.03.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.008 | vasa recta | A08.1.03.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.03.009 | capsular branches | A08.1.03.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.04.001 | intrarenal veins | A08.1.01.001 | vascular tissue network; lumen and junction review required | `ZA:Intrarenal veins of left kidney`<br>`ZA:Intrarenal veins of right kidney` | Unbuilt |
| A08.1.04.002 | interlobar veins | A08.1.04.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.04.003 | arcuate veins | A08.1.04.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.04.004 | cortical radiate veins | A08.1.04.003 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.04.005 | venulae rectae | A08.1.04.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.04.006 | stellate veins | A08.1.04.001 | vascular tissue network; lumen and junction review required | Unbuilt | Unbuilt |
| A08.1.05.001 | renal pelvis | A08.1.01.001 | collecting conduit with walls and connected lumen | `ZA:Renal pelvis.l`<br>`ZA:Renal pelvis.r` | `HRA:VH_F_renal_pelvis_R` (hidden reference)<br>`HRA:VH_F_renal_pelvis_L` (hidden reference) |
| A08.1.05.002 | branching type | chapter-08 | variant classification, not separate geometry | Unbuilt | Unbuilt |
| A08.1.05.003 | major calices | A08.1.05.001 | tissue region or regional assembly | Unbuilt | `HRA:VH_F_major_calyx_R_a` (hidden reference)<br>`HRA:VH_F_major_calyx_R_b` (hidden reference)<br>`HRA:VH_F_major_calyx_R_c` (hidden reference)<br>`HRA:VH_F_major_calyx_L_a` (hidden reference)<br>`HRA:VH_F_major_calyx_L_b` (hidden reference)<br>`HRA:VH_F_major_calyx_L_c` (hidden reference)<br>`HRA:VH_F_major_calyx_L_d` (hidden reference) |
| A08.1.05.004 | superior calyx | A08.1.05.003 | collecting conduit with walls and connected lumen | Unbuilt | Unbuilt |
| A08.1.05.005 | middle calyx | A08.1.05.003 | collecting conduit with walls and connected lumen | Unbuilt | Unbuilt |
| A08.1.05.006 | inferior calyx | A08.1.05.003 | collecting conduit with walls and connected lumen | Unbuilt | Unbuilt |
| A08.1.05.007 | minor calices | A08.1.05.001 | tissue region or regional assembly | Unbuilt | `HRA:VH_F_minor_calyx_R_a` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_b` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_c` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_d` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_e` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_f` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_g` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_h` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_i` (hidden reference)<br>`HRA:VH_F_minor_calyx_R_j` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_a` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_b` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_c` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_d` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_e` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_f` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_g` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_h` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_i` (hidden reference)<br>`HRA:VH_F_minor_calyx_L_j` (hidden reference) |
| A08.1.05.008 | ampullary type | chapter-08 | variant classification, not separate geometry | Unbuilt | Unbuilt |
| A08.1.05.009 | adventitia | A08.1.05.001 | distinct tissue envelope or wall layer; outer shell alone insufficient | Unbuilt | Unbuilt |
| A08.1.05.010 | muscular layer | A08.1.05.001 | distinct tissue envelope or wall layer; outer shell alone insufficient | Unbuilt | Unbuilt |
| A08.1.05.011 | mucosa | A08.1.05.001 | distinct tissue envelope or wall layer; outer shell alone insufficient | Unbuilt | Unbuilt |

## Assembly plan

Each native kidney in its fixed current body frame. Inspect complete source surface, sinus and real collecting source before deriving tissues.
Finite fibrous capsule with a true hilar termination; cortex with columns between separately modeled medullary pyramids. A single eroded inner shell is not renal medulla. Calyces/pelvis/ureter require finite walls and connected urine space. Microscopic nephron organization needs a separate efficient detail representation.
Keep skeleton, suprarenal assemblies, stomach/pancreas and native renal vessels fixed. Retain literal proven source interfaces. Vessels are allowed to enter renal tissue; they must not enter collecting urine space. Hilum/sinus are openings/cavity regions, never blocking solids. Fat and renal fascia remain distinct tissues.
Male first. Transfer the coherent repaired assembly with explicit size/pose fit and fixed recipient hilar/ureter landmarks; never expose a second overlaid donor kidney.
Current inventory is not anatomical acceptance. Hidden HRA sources remain local and are not enabled or copied as native finished anatomy.

Current geometry, hashes and source diagnostics: `reports/kidneys-current-v249.json`. The 2026-10-03 reports are historical evidence and do not establish current renal acceptance.

## Historical first finite trial

The current male whole surfaces each have zero tested self-crossings. The female source left/right retain 122/43 crossings, and the right retains 12 open edges. Male source external, posterior, collecting and two real section views were inspected in the existing Blender window.

The left candidate has a finite capsule, unpartitioned remaining renal tissue, and adventitia/muscular/mucosal collecting-wall cells. Actual packed topology, shared interfaces, outward normals and six finite tissue section/peel renders were reviewed. It is uninstalled; cortex, pyramids, papillae, calyceal subdivisions and complete vascular/ureteral interfaces are unfinished. The source calyceal ends remain coarse and do not establish papillary openings.

The bilateral v251 trial was stopped at the right collecting outlet: a computational cap intersects kidney tissue. Do not remove that guard or render the cap. Establish the actual native ureter continuation. The successful left construction is retained separately as v252 with packed replay v253 and full-frame render v255. See `reports/renal-first-candidate-review-v256.json`. Both installed atlases stay unchanged.

## Historical connected finite construction

The later v257 carrier joins both native male pelvis outlets to their ureters. The right collar is rounded with source-section area retained; both physiological urine passages remain open. The v264 construction contains 22 finite cells: capsule, explicitly unpartitioned renal tissue, three pelvis coats and three coats for each of the abdominal and pelvic ureter regions on each side. The native regional plane and constructed thicknesses remain provisional. This construction does not establish cortex, pyramids, papillae, microscopic anatomy or an intramural bladder junction.

The actual packed banks replay with reciprocal positions and normals, closed tissue cells, no detected cell self-crossings, one urine space per side and no stored disk across the pelvis–ureter opening. The two banks total 672,166 compressed bytes and add no images. Ten exported surface, peel and finite-section views were rendered and inspected in the existing Blender window (v266–267); the prior scene was preserved.

Installation and anatomical acceptance remain pending. Native blood vessels enter the collecting urine space. Holding every internal source-cylinder end fixed cannot repair all of these conflicts: v268 finds eight left venous and four right arterial rim vertices already inside the closed urine classifier. The [renal venous endocast study](https://pubmed.ncbi.nlm.nih.gov/2231877/) was read at abstract level: venous branches pass around both sides of the collecting tree, so moving every vein into one anterior plane is unsupported. A coupled vascular/collecting correction must preserve proven external hilar attachments and keep blood outside the innermost collecting surface. Both current atlas catalogues remain unchanged by this construction.
