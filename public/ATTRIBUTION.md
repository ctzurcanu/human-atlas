# Anatomy data attribution

BodyParts3D, © The Database Center for Life Science licensed under CC Attribution 4.0 International.

- License: https://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html (updated 2025-02-27)
- Dataset: https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html
- License terms: https://creativecommons.org/licenses/by/4.0/
- Source geometry: `isa_BP3D_4.0_obj_99.zip`, BodyParts3D 4.0.
- English names and relationships: IS-A and PART-OF concept, element, and inclusion tables from the same archive.
- Publication: Mitsuhashi et al. (2009), BodyParts3D: 3D structure database for anatomical concepts. https://doi.org/10.1093/nar/gkn613

Adaptations: axes and units converted from millimeters/Z-up to meters/Y-up; translated to rest at the stage; geometry simplified using meshoptimizer with 0.2% relative error limit per structure; normals quantized to signed 16-bit; packed into binary chunks; curated display system groupings and colors. The source contains 2,234 individual OBJ meshes; all remain represented. The combined hierarchy contains 3,432 named FMA concepts, which may reference multiple meshes. Original source identity is preserved in the manifest.

Source OBJ comments mention an older CC BY-SA 2.1 Japan license. The official current database license linked above supersedes that legacy text and explicitly permits redistribution and adaptation under CC BY 4.0.

BodyParts3D represents an adult male reference anatomy based on TARO MRI and anatomical illustration refinements. It is not a complete model of every possible human anatomical structure or variation. This interface is educational and is not a clinical tool.

## Female anatomy

Kristen Browne and Heidi Schlehlein, Human Reference Atlas / HuBMAP, *3D Reference Organ Set for Female v1.10* (2026), with eight pelvic pieces retained from v1.5 (2023). CC BY 4.0.

- Reference library: https://humanatlas.io/3d-reference-library
- Dataset: https://purl.humanatlas.io/ref-organ/united-female/v1.10
- Original GLB: https://cdn.humanatlas.io/digital-objects/ref-organ/united-female/v1.10/assets/3d-vh-f-united.glb
- Earlier pelvic geometry: https://cdn.humanatlas.io/digital-objects/ref-organ/united-female/v1.5/assets/3d-vh-f-united.glb
- License: https://creativecommons.org/licenses/by/4.0/

Female lower-limb muscles: Thor E. Andreassen, Donald R. Hume, Landon D. Hamilton, Karen E. Walker, Sean E. Higinbotham, and Kevin B. Shelburne, *Three Dimensional Lower Extremity Musculoskeletal Geometry of the Visible Human Female and Male*, Scientific Data 10, 34 (2023). Only female geometry is used. CC BY 4.0.

- Female dataset and license: https://digitalcommons.du.edu/visiblehuman/1/
- Publication: https://doi.org/10.1038/s41597-022-01905-2
- Underlying imagery: Visible Human Project, U.S. National Library of Medicine.

Geometry adaptation: slorksmo/Human-Atlas, revision `5bb5713aab18d7fe9380c3339eb09f173491ea06`: https://github.com/slorksmo/Human-Atlas/tree/5bb5713aab18d7fe9380c3339eb09f173491ea06

Adaptations inherited from that package: HRA node transforms baked into meter/Y-up stage coordinates; coincident vertices welded; normals averaged and quantized; geometry simplified with a 0.2% relative error bound per structure. The 76 female lower-limb muscle surfaces are fitted by group using rotation, uniform scale and translation to the HRA pelvis and leg bones. Reported bone-fit residuals range from 17 to 36 mm. This placement is approximate and is not a validated single-scan assembly.

Additional changes here: all 180 male-derived meshes in the adaptation are excluded; buffers are repacked so those surfaces are not distributed. Two HRA rectus femoris surfaces are replaced by their female-study counterparts to avoid duplicate anatomy. Brain regions are grouped under Nervous system; female-study muscles under Muscles. Each piece retains a source link. The result contains 1,038 selectable meshes and 1,253 searchable concepts, including 90 muscle pieces. Eight pregnancy-reference pieces remain hidden by default.

This is still a partial reference. It lacks the skull, ribs, shoulder, arm, hand and foot bones, and most upper-body muscles. It must not be described as equivalent in completeness to the BodyParts3D male atlas.

## Detailed male reference

Anatomy adapted by **Brian Pridgen** from **Z-Anatomy** and **Open 3D Model**, with BodyParts3D and other upstream contributors. Full upstream attribution and licensing qualifications are preserved verbatim in [MALE-ANATOMY-NOTICE.md](MALE-ANATOMY-NOTICE.md).

- Source anatomy: https://anatomy-atlas.brianp.chatgpt.site/
- Source notice: https://anatomy-atlas.brianp.chatgpt.site/NOTICE.md
- Z-Anatomy: https://github.com/Z-Anatomy/Models-of-human-anatomy
- Open 3D Model: https://anatomytool.org/open3dmodel
- Open 3D Model institutions: Leiden UMC, UMC Utrecht, Maastricht UMC, University of Leuven KULAK, Amsterdam UMC, Radboud UMC Nijmegen, and University of Ghent.

The adapted anatomy and catalog are CC BY-SA 4.0 except the two kidneys (Lissie Cowley, CC BY-NC 4.0) and the two cochleae and two vestibules (University of Dundee School of Medicine, CC BY-NC-SA 4.0). Retain each component’s license when redistributing these meshes or outputs containing them. The upstream notice preserves qualifications about the underlying white-matter sources; this import does not claim new independent clearance.

Changes for Human Atlas: source regional GLBs decoded; node transforms baked, including inverse-transpose normal transforms and winding correction where necessary; translated upward by 31 mm to sit on the stage. The optimized option simplifies geometry at a maximum 0.2% relative error per named structure. The full-resolution option preserves all 11,647,803 imported triangles and includes the derived outer body surface, with no mesh simplification. Both options retain this viewer’s system colors rather than source textures; normals quantized; binary buffers repacked and gzip compressed. Source-defined muscle attachment annotations, topographical regions, fascia, and 27 source-designated schematic structures are separate optional layers. The schematics retain their non-anatomical source placement and are hidden by default. Source names and groups become searchable concepts; unique spelling/spacing differences are resolved, duplicate representations are excluded, and unidentifiable geometry is omitted. Per-mesh source, license, and mirroring information remains in the manifest.

The import contains 4,391 individually selectable pieces, including all 769 Open 3D Model records. The source’s adapted upper limbs include mirrored structures; neither this import nor those adaptations establish independent anatomical validation. The original BodyParts3D atlas remains available because coverage and subdivision differ.

Geometry identity is pinned in `scripts/data/male-source-index.json`. The served abdomen-pelvis file includes the credited kidney restoration and differs from the stale hash in the source page index; both hashes and the verification note are retained. No proprietary viewer code, styling, or text definitions were imported. Human Atlas’s conversion script and UI code remain MIT; those terms do not override the separate anatomy licenses.
