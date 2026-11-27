# Spleen: TA98 inventory and first-version assembly plan

Priority 7 follows installation of the first pancreatic version in both TA98 models. Start in the male, then fit/transfer the complete assembly to the female. Preserve the skeleton, accepted stomach and current pancreatic carrier.

The whole assembly comprises a finite serosal coat, fibrous capsule, connected trabecular support and parenchyma. Red/white pulp, sinuses, penicilli and lymphoid nodules require distinct microscopic representations or an explicit detail model, not a generic outer shell. Surface impressions, poles and borders are surface/curve features; the hilum is a region with real vascular openings and cannot be modeled as a blocking solid. Accessory spleen is an optional variant, not a second default spleen.

Adjoining interfaces: diaphragm, stomach/short gastric routes, pancreatic tail, left kidney, colic flexure, splenic artery/vein and hilar branches. Keep the existing portal/splenic venous confluence fixed while defining separate hilar refit constraints. Use the current dependency forest. Review the complete source exterior and real sections before partitioning or constructing coats. Local reference meshes remain local and are not copied precisely or published.

| TA98 | Term | Parent | Representation | Pre-edit male source candidates | Pre-edit female source candidates |
|---|---|---|---|---|---|
| A13.2.01.001 | spleen | A13.2.00.001 | whole-organ assembly | `ZA:Spleen` (1330 triangles; catalogue source; runtime eligibility still to review) | `LOCAL:female:spleen:1038` (598 triangles; catalogue source; runtime eligibility still to review) |
| A13.2.01.002 | fibrous capsule | A13.2.01.001 | finite fibrous tissue coat | No distinct source representation | No distinct source representation |
| A13.2.01.003 | splenic trabeculae | A13.2.01.001 | internal connective trabecular tissue | No distinct source representation | No distinct source representation |
| A13.2.01.004 | splenic pulp | A13.2.01.001 | parenchymal tissue group | No distinct source representation | No distinct source representation |
| A13.2.01.005 | red pulp | A13.2.01.004 | parenchymal tissue | No distinct source representation | No distinct source representation |
| A13.2.01.006 | white pulp | A13.2.01.004 | lymphoid tissue distribution | No distinct source representation | No distinct source representation |
| A13.2.01.007 | diaphragmatic surface | A13.2.01.001 | outer surface region | No distinct source representation | No distinct source representation |
| A13.2.01.008 | visceral surface | A13.2.01.001 | visceral surface group | No distinct source representation | No distinct source representation |
| A13.2.01.009 | renal impression | A13.2.01.008 | surface impression | No distinct source representation | No distinct source representation |
| A13.2.01.010 | gastric impression | A13.2.01.008 | surface impression | No distinct source representation | No distinct source representation |
| A13.2.01.011 | colic impression | A13.2.01.008 | surface impression | No distinct source representation | No distinct source representation |
| A13.2.01.012 | pancreatic impression | A13.2.01.008 | surface impression | No distinct source representation | No distinct source representation |
| A13.2.01.013 | anterior extremity | A13.2.01.001 | extremity landmark/surface region | No distinct source representation | No distinct source representation |
| A13.2.01.014 | posterior extremity | A13.2.01.001 | extremity landmark/surface region | No distinct source representation | No distinct source representation |
| A13.2.01.015 | inferior border | A13.2.01.001 | boundary curve | No distinct source representation | No distinct source representation |
| A13.2.01.016 | superior border | A13.2.01.001 | boundary curve | No distinct source representation | No distinct source representation |
| A13.2.01.017 | splenic hilum | A13.2.01.001 | hilar region and openings; no solid plug | No distinct source representation | No distinct source representation |
| A13.2.01.018 | serosa | A13.2.01.001 | finite serosal tissue coat, with hilar exclusion | No distinct source representation | No distinct source representation |
| A13.2.01.019 | splenic sinus | A13.2.01.001 | microscopic vascular sinus wall/lumen | No distinct source representation | No distinct source representation |
| A13.2.01.020 | penicilli | A13.2.01.001 | microscopic arteriolar branch network | No distinct source representation | No distinct source representation |
| A13.2.01.021 | splenic lymphoid nodules | A13.2.01.001 | microscopic lymphoid nodules | No distinct source representation | No distinct source representation |
| A13.2.01.022 | accessory spleen | A13.2.01.001 | optional accessory-organ variant | No distinct source representation | No distinct source representation |

Authority: https://ifaa.unifr.ch/Public/EntryPage/ShowTA98.html. Full source identities, hashes, source topology and crossing checks: `reports/ta98-spleen-baseline-v198.json`. Inventory does not establish anatomical acceptance.

## Finite-tissue construction, before modification

Actual source exterior, hilar vessels, renal/pancreatic context and coronal sections were rendered in the existing Blender window (`spleen-source-render-v199.json`). Both organ surfaces are closed genus zero with zero tested crossings. The coarse visceral grooves and independently open vascular guides remain provisional. Preserve the source exterior, skeleton, pancreatic tail and repaired portal confluence; do not introduce a second spleen.

The first finite assembly has four tissue controls: serosa, fibrous capsule, trabeculae and remaining pulp. It uses one source carrier and canonical reciprocal interfaces. A nominal 0.24 mm combined coat, including 0.04 mm serosa, is a declared construction setting, not a measured patient thickness; the bounded offset field must pass the continuous-prism and complete stored-surface checks. The hilar exclusion follows the medial source surface near the actual hilar vascular branches and needs rendered review. It is a serosal termination, not an artificial solid representing the hilum. The capsule remains continuous with the internal fibrous support. Branching connective tissue partially subdivides pulp; it must not form sealed artificial organ lobules. The initial trabecular distribution is illustrative and must stay explicitly provisional. No new vessel or lumen claim follows from it.

The existing vascular guides do not yet have verified coupled inner/outer walls, and their terminal interfaces remain unfinished. Keep that obligation open; do not model computational disks as anatomy. Red/white pulp, lymphoid nodules, penicilli and sinuses remain a separate microscopic construction task. Reference morphology: [Columbia spleen histology](https://histologylab.ctl.columbia.edu/lab08/spleen/), including its actual nodule micrograph; [Leeds spleen histology](https://histology.leeds.ac.uk/home/lymphoid/lymph_spleen/). No reference image pixels or geometry are copied.

## Installed first finite version — both local TA98 models

`spleen-install-v203`, `spleen-peeling-v204` and `spleen-appearance-v205` install four finite tissue controls: serosa (.018), fibrous capsule (.002), pulp (.004), and trabeculae (.003). Together with the whole-organ assembly, this is **5 of 22 terms represented at this first construction level**, not 22 completed terms. The remaining **17 terms** above still need their own appropriate representations. The microscopic red/white pulp, sinus, penicilli and nodular morphology have not been fabricated or counted as finished.

The original male/female exteriors remain literal: 1,330/598 triangles, one closed genus-zero surface, with no duplicate organ. All eight stored tissue cells have zero tested self-intersections and valid closed boundaries; the pulp stays connected. Multiple trabecular extensions are allowed only when every component attaches to the same capsule. The constructed distribution is sparse and illustrative. The source grooves, polygonal serosal termination and unverified vascular guides remain anatomical review obligations. Hilar coupling is not accepted. Native vessels, skeleton, stomach, pancreatic geometry and visit coordinates remain fixed.

Each assembly uses one shared position/normal bank with reciprocal tissue interfaces and ordinary labels, components, visibility and peel controls. Smooth shading reduces the packed banks to **231,914 / 167,152 bytes**, **97,420 / 72,660 gzip bytes**, male/female. No image is added. Live section fills use the tissue's own appearance instead of the generic green lymphatic style. Independent actual-bank replay is in `reports/spleen-runtime-replay-v206.json`.

The source alternatives are suppressed and retained locally. A URL-state repair keeps catalogue-default hiding of a retired source from hiding its replacements; explicit user hides still follow the replacement organ. Tissue hide/peel, adjoining bowel peels, camera intent and shared explosion layout tests pass. This is local installation progress; it is not complete spleen acceptance or public deployment.


## Hilar interface source inventory before vascular reconstruction

All 21 applicable descendants of the splenic artery and vein are inventoried below from the authoritative term tree. These are source candidates, not proof of lumen continuity or finished branch anatomy. Preserve the repaired portal confluence, its shared source rims and the current whole-body coordinate frames. Classifier disks may close a numerical envelope but must never be exported as blood-space barriers. Source serial station ordering must be proved before deriving paired inner profiles. Finite vessel mouths close only their tissue annuli. Hilar contacts must be coupled to all four splenic tissues, with terminal continuations and unresolved microcirculation kept explicit.

| TA98 | Term | Parent | Male source candidates | Female source candidates |
|---|---|---|---|---|
| A12.2.12.040 | splenic artery | A12.2.12.012 | `ZA:Splenic artery` | `LOCAL:female:splenic_artery:1679` |
| A12.2.12.041 | pancreatic branches | A12.2.12.040 | `FJ3419` (hidden source)<br>`FJ3430` (hidden source)<br>`FJ3433` (hidden source)<br>`FJ3444` (hidden source) | No distinct source representation |
| A12.2.12.042 | dorsal pancreatic artery | A12.2.12.041 | `FJ3430` (hidden source) | No distinct source representation |
| A12.2.12.043 | inferior pancreatic artery | A12.2.12.041 | `FJ3444` (hidden source) | No distinct source representation |
| A12.2.12.044 | prepancreatic artery | A12.2.12.041 | No distinct source representation | No distinct source representation |
| A12.2.12.045 | greater pancreatic artery | A12.2.12.041 | `FJ3433` (hidden source) | No distinct source representation |
| A12.2.12.046 | artery to tail of pancreas | A12.2.12.041 | `FJ3419` (hidden source) | No distinct source representation |
| A12.2.12.047 | left gastro-omental artery | A12.2.12.040 | `FJ3501` (hidden source) | No distinct source representation |
| A12.2.12.048 | gastric branches | A12.2.12.047 | No distinct source representation | No distinct source representation |
| A12.2.12.049 | omental branches | A12.2.12.047 | No distinct source representation | No distinct source representation |
| A12.2.12.050 | short gastric arteries | A12.2.12.040 | No distinct source representation | No distinct source representation |
| A12.2.12.051 | splenic branches | A12.2.12.040 | No distinct source representation | No distinct source representation |
| A12.2.12.052 | posterior gastric artery | A12.2.12.040 | No distinct source representation | No distinct source representation |
| A12.3.12.028 | splenic vein | A12.3.12.001 | `ZA:Splenic vein` | `LOCAL:female:splenic_vein:1796` |
| A12.3.12.029 | pancreatic veins | A12.3.12.028 | No distinct source representation | No distinct source representation |
| A12.3.12.030 | short gastric veins | A12.3.12.028 | No distinct source representation | No distinct source representation |
| A12.3.12.031 | left gastro-omental vein | A12.3.12.028 | `ZA:Left gastro-omental vein` | No distinct source representation |
| A12.3.12.032 | inferior mesenteric vein | A12.3.12.028 | `ZA:Inferior mesenteric vein` | `LOCAL:female:inferior_mesenteric_vein:1754` |
| A12.3.12.033 | left colic vein | A12.3.12.032 | `FJ3400` (hidden source)<br>`FJ3429` (hidden source)<br>`FJ3495` (hidden source) | `HRA:VH_F_left_colic_vein` (hidden source) |
| A12.3.12.034 | sigmoid veins | A12.3.12.032 | `ZA:Sigmoid veins` | `LOCAL:female:sigmoid_veins:1753` |
| A12.3.12.035 | superior rectal vein | A12.3.12.032 | `FJ3657` (hidden source) | `LOCAL:female:superior_rectal_vein:1789` |

Measured source entries, mouth containment, crossings and hashes: `reports/splenic-hilum-v207.json`.

The male native artery is seven disconnected open tube strips. Their original serial sections can be recovered from connectivity even though their vertex table is not in station order. The first envelope union has one communicating blood volume, but six source branch-root disks only partly overlap the parent lumen. Those partial disks would act as false internal baffles if exported. Before installation, replace each internal source end closure by an explicitly constructed outer/inner junction sleeve ending strictly inside the parent lumen. This adds local tissue at the existing junction; it does not relocate source profiles, shrink their recorded section or move the organ. The source 0.88 lumen/outer radial ratio is a provisional construction parameter. Only true external mouths receive annular wall ends, with numerical fluid classifier disks outside the anatomical wall. Verify no temporary fluid disk survives, the eight mouths stay open, the whole stored wall is a single finite genus-seven boundary, and actual junction/cut renders look plausible. Preserve unresolved clinical calibre, adjacent walls and downstream microcirculation as review obligations.


## First native arterial wall installation

Male `ZA:Splenic artery` now uses the single finite eight-mouth network from `reports/splenic-arterial-wall-v211.json`, independently replayed in v213 and installed by v215. Six branch roots communicate; no temporary blood-space disks are in the tissue mesh. All source outer profile coordinates and every other anatomical part remain fixed. Eight annuli, one communicating blood space, strict topology and zero detected stored self-crossings pass. Seven actual packed Blender views and the installed live section were inspected. Source branch collars remain angular, calibre is a construction value, and splenic/venous tissue passages remain unfinished. This is a first improvement, not complete hilar or whole-spleen acceptance. Female retains its native artery and the four installed finite splenic tissues. No new image assets or duplicate vessel instance.
