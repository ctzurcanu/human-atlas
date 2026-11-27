# Suprarenal glands: priority 8

Current source inventory before geometry changes. Each entry remains an implementation obligation; a label or outer carrier does not establish internal tissue.

| TA98 | Term | Parent | Representation | Male meshes | Female meshes |
|---|---|---|---|---|---|
| A11.5.00.001 | suprarenal gland | A11.0.00.000 | whole organ | `ZA:Suprarenal gland.l`<br>`ZA:Suprarenal gland.r` | `LOCAL:female:l_adrenal_gland:190`<br>`LOCAL:female:r_adrenal_gland:191` |
| A11.5.00.002 | anterior surface | A11.5.00.001 | surface region | Unbuilt | Unbuilt |
| A11.5.00.003 | posterior surface | A11.5.00.001 | surface region | Unbuilt | Unbuilt |
| A11.5.00.004 | renal surface | A11.5.00.001 | surface region | Unbuilt | Unbuilt |
| A11.5.00.005 | superior border | A11.5.00.001 | boundary curve | Unbuilt | Unbuilt |
| A11.5.00.006 | medial border | A11.5.00.001 | boundary curve | Unbuilt | Unbuilt |
| A11.5.00.007 | cortex | A11.5.00.001 | tissue volume | Installed finite left/right cells; provisional histological contour | Installed finite left/right cells; provisional histological contour |
| A11.5.00.008 | medulla | A11.5.00.001 | tissue volume | Installed finite left/right cells; provisional central core | Installed finite left/right cells; provisional central core |
| A11.5.00.009 | hilum | A11.5.00.001 | vascular passage, not a solid | Unbuilt | Unbuilt |
| A11.5.00.010 | central vein | A11.5.00.001 | finite vessel wall and lumen | Unbuilt | Unbuilt |
| A11.5.00.011 | accessory suprarenal glands | A11.5.00.001 | variant organ; only if evidenced | Unbuilt | Unbuilt |

## Assembly and interfaces

Each native gland in its current skeleton, kidney, diaphragm and vascular frame. Inspect and correct the complete outer surface and renal pole relationship before tissue partitions.
Capsule is an additional tissue control; cortex and medulla need distinct finite volumes and a common reciprocal interface. Named surfaces partition an existing exterior; borders are curves. Hilum remains an opening. A central vein needs an adrenal source and a communicating finite lumen; no hepatic-lobule source.
Adrenal and kidney are separate organs. Surrounding vessels may enter the gland but need separate wall/blood space. No shared adrenal/renal faces, random accessory copies, skeleton movement or suppression of unrelated sources. Superseded whole-gland sources stay retained and hidden behind their explicit replacement tissues.
Female first, then size and placement in male; retain distinct left/right shapes. Provisional boundaries and unfinished vascular interfaces remain explicit.

Current coordinates, fixed neighboring meshes and exact source diagnostics are in `reports/suprarenal-current-v216.json`. No geometry or acceptance change in this stage.

## First clearance trial

The actual sixteen native surface/section views show the adrenal renal-facing surface crossing the renal pole. Test a bounded superior translation of the whole female gland, first, without rotating or reshaping it. Search 0–6 mm, retaining the first crossing-free interval only as a candidate. Then check every staged neighboring surface and exact stored topology/normal orientation. Any new nonphysiological neighbor intersection rejects installation. Existing vessel entry locations are unverified and require explicit refitting; do not claim that translation establishes a central vein or hilar continuity. No kidney, diaphragm, accepted organ or skeleton moves. Capsule/cortex/medulla construction follows whole-surface and packing review; provisional dimensions remain labeled as construction values.


All sixteen current packed adrenal/renal views were inspected in existing Blender with the original scene preserved. The first bounded female superior-translation trial did not remove all renal surface crossings within 6 mm on either side; it is uninstalled. Do not extend this into an unbounded displacement search. Next construct a local renal-facing shape correction with fixed external landmarks and complete neighboring checks, then derive the capsule/cortex/medulla interfaces. See `reports/female-suprarenal-clearance-v218.json`. OpenStax Figure 17.17 and its histology inset were actually inspected in the browser; no pixels, UVs or reference mesh were copied.

## Local renal-facing reconstruction

Retain the original whole-gland pose and every exterior triangle outside the renal contact patch. Because the female renal source has unverified internal topology, use its convex occupied envelope only as a conservative clearance constraint, not as a replacement kidney or an assertion of its actual tissue volume. Offset its supporting planes by 0.3 mm as an explicit construction clearance. Intersect the original healthy gland solid with the exterior of that constraint. Reject a disconnected result, new tunnels, self-crossings or removal of more than 25% of native volume; that bound limits a proposed local correction and is not a physiological measurement. Verify the actual Float32 candidate, all fixed neighbors and retained source face lineage, and inspect complete renal-facing contours and cuts before deriving the capsule/cortex/medulla. Source renal vessels and adrenal entries remain explicit refit obligations. This does not authorize moving the kidney or replacing its anatomy with a convex hull.

## Finite tissue construction

The corrected whole surfaces are healthy, single closed genus-zero boundaries. Both actual renal-facing cuts and complete contours were inspected. The proposed first tissue version uses one exterior plus two shared internal interfaces: a connective-tissue capsule, cortex and medulla. Derive the internal interfaces by morphological erosion with a small eight-segment geodesic kernel, not by scaling the gland about an arbitrary point in its concavity. The 60 µm capsule kernel and approximately 20% medullary volume target are explicitly illustrative construction values, not measured histology or anatomical acceptance. The cortex/medulla boundary is a provisional central-core model pending detailed histological shape and central-vein coupling. Every shared interface has reciprocal triangle ownership and opposite normal signs in one buffer bank; no new UV or image assets. Check stored cells, nesting, intersections, volumes and all shared attributes, then review actual peeling and sections. Remaining source contacts with diaphragm, aorta, posterior muscles, connective tissue and spine are recorded in v219; the renal reconstruction does not certify their anatomy. Do not infer a finished central vein, hilum, named surface or accessory gland from these three tissues.

## Local first versions installed in both TA98 models

Female v226 and male v238 replace the displayed native whole-gland carriers with six finite cells per body: capsule, cortex and medulla on each side. Both retain the original source organs hidden for review. One shared bank per body stores reciprocal interfaces and normal signs; female 203,009 raw / 89,371 gzip bytes, male 1,002,953 raw / 441,679 gzip bytes, with no new images. The female construction and controls were reused on the male native contours, retaining their separate left/right outlines and body size rather than translating another body's organs without fit checks.

All skeleton, kidney and other-organ coordinates remain literal. Renal crossing pairs decreased female 76/54 to 0/0 and male 248/202 to 0/0. The conservative renal clearance removes female 4.014%/1.095% and male 21.815%/22.272% of the original occupied gland volume. These changes and the 0.3 mm construction clearance are not clinical measurements. Source neighboring contacts remain unresolved. Existing Blender scenes were preserved while complete contours, peeling and finite tissue cuts were rendered and inspected.

Actual packed cells pass topology, self-crossing, normal orientation, volume closure and literal reciprocal-interface checks. The male Boolean overlap diagnostic retains a maximum signed residual of 1.301e-17 cubic metres within the recorded kernel-scale numerical guard; it is not represented as exact zero. See local v231–238 evidence. The capsule/cortex/medulla are individually selectable in TA98, can be peeled and use their own section colors. Whole-gland native functional references and renal adjacency survive replacement; source frame and laterality remain guarded. No new adrenal vessel junction is claimed.

Three of the eleven listed terms now have displayed first-version geometry: whole gland, cortex and medulla; capsule is additional. The eight other entries remain explicit obligations, with accessory glands conditional on variant evidence. Central veins and hilum, named surfaces/borders, histological contour refinement and wider packing are unfinished. This is local installation progress, not complete anatomical acceptance or finished priority 8.

## Central-vein construction and source interfaces

The current v240 inventory includes the complete authority entries for central vein A11.5.00.010 and external left/right suprarenal veins A12.3.09.011/.013; all three are leaves. Native female renal veins cross the adrenal exterior in 26 left/15 right face pairs but do not cross the constructed medulla. This is not a verified central-vein representation. Hidden male BodyParts3D suprarenal sources FJ3480/FJ3580 remain distinct reference frames, not installed native routes. The actual eight venous/gland placement views were inspected in existing Blender without saving or changing the original scene.

Before installation, construct an irregular finite intramedullary wall, a continuous blood space, a genuine hilar exit through the complementary gland tissues, and an external wall/lumen interface to the appropriate native parent. Left drainage is to the left renal vein; right drainage is to the vena cava. Existing renal branches must not be relabeled by proximity. Proximal sinusoidal tributary ports remain distinct from the main external venous junction. A provisional medial route or an annular tube alone does not establish the completed central-vein anatomy.

The v242 initializer uses clearance-constrained routes within the actual medullary volumes and two asymmetric longitudinal wall thickenings. Its 0.24 mm nominal inner radius and 65–155 micrometre wall are illustrative construction values, not measured clinical calibre. Both annular mouths remain open; diagnostic fluid closures are not anatomical wall disks. The candidate remains uninstalled until endocrine complements, hilar outlets, upstream tributary scope, adjoining walls/lumens, actual renders and anatomical course are reviewed.

References actually read: Barenthin 1975, PMID1229251 (primary light/electron microscopy; incomplete eccentric longitudinal muscle); the primary adrenalectomy series PMID18816449 (usual drainage and variation); and Histology Guide MH155a, human H&E specimen from the University of Minnesota, with the Medullary Veins view and slide metadata inspected in the browser. No pixels, UVs or tissue mesh are copied. The medRxiv 2023.01.04.23284197 preprint was read for corroboration only; its figure rendering was unavailable and it is not treated as peer-reviewed or as a clinical course measurement.

The four v242 finite-wall candidates replay from their actual banks and fit wholly within their installed medullae, with two open annular rims each, no detected wall self-crossings and one diagnostic blood-space component. Total compressed payload is 235,122 bytes with no new runtime images. The independent replay retains the signed Boolean residue and kernel-scale numerical guard rather than claiming exact zero overlap. All four v245 Blender sections were inspected with flat section normals, and the existing scene was preserved. They expose the material remaining task: the diagnostic medullary terminations still need actual upstream/hilar interfaces, and the external courses still need verified body-frame placement and parent lumens. The candidates are uninstalled; A11.5.00.010 remains unfinished in both catalogues. See local v240–246 evidence.
