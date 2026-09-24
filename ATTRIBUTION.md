# Human Atlas attribution and component licenses

The application code is distributed under [GNU GPLv3](LICENSE). Original upstream application code is copyright © 2026 ashemag; its [MIT notice](LICENSES/MIT-upstream.txt) remains in the repository. Anatomical geometry, metadata, fonts, and dependencies retain their respective licenses below.

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

Anatomy adapted by **Brian Pridgen** from **Z-Anatomy** and **Open 3D Model**, with BodyParts3D and other upstream contributors. Detailed upstream credits and licensing qualifications appear below.

- Source anatomy: https://anatomy-atlas.brianp.chatgpt.site/
- Z-Anatomy: https://github.com/Z-Anatomy/Models-of-human-anatomy
- Open 3D Model: https://anatomytool.org/open3dmodel
- Open 3D Model institutions: Leiden UMC, UMC Utrecht, Maastricht UMC, University of Leuven KULAK, Amsterdam UMC, Radboud UMC Nijmegen, and University of Ghent.

The adapted anatomy and catalog are CC BY-SA 4.0 except the two kidneys (Lissie Cowley, CC BY-NC 4.0) and the two cochleae and two vestibules (University of Dundee School of Medicine, CC BY-NC-SA 4.0). Retain each component’s license when redistributing these meshes or outputs containing them. The upstream notice preserves qualifications about the underlying white-matter sources; this import does not claim new independent clearance.

Changes for Human Atlas: source regional GLBs decoded; node transforms baked, including inverse-transpose normal transforms and winding correction where necessary; translated upward by 31 mm to sit on the stage. The optimized option simplifies geometry at a maximum 0.2% relative error per named structure. The full-resolution option preserves all 11,647,803 imported triangles and includes the derived outer body surface, with no mesh simplification. Both options retain this viewer’s system colors rather than source textures; normals quantized; binary buffers repacked and gzip compressed. Source-defined muscle attachment annotations, topographical regions, fascia, and 27 source-designated schematic structures are separate optional layers. The schematics retain their non-anatomical source placement and are hidden by default. Source names and groups become searchable concepts; unique spelling/spacing differences are resolved, duplicate representations are excluded, and unidentifiable geometry is omitted. Per-mesh source, license, and mirroring information remains in the manifest.

The import contains 4,391 individually selectable pieces, including all 769 Open 3D Model records. The source’s adapted upper limbs include mirrored structures; neither this import nor those adaptations establish independent anatomical validation. The original BodyParts3D atlas remains available because coverage and subdivision differ.

Geometry identity is pinned in `scripts/data/male-source-index.json`. The served abdomen-pelvis file includes the credited kidney restoration and differs from the stale hash in the source page index; both hashes and the verification note are retained. No proprietary viewer code, styling, or text definitions were imported. Human Atlas’s conversion script and UI code are distributed under [GPLv3](LICENSE); those terms do not override the separate anatomy licenses. The upstream MIT notice for original application code is preserved in [LICENSES/MIT-upstream.txt](LICENSES/MIT-upstream.txt).

---

## Detailed male source credits and licenses

The source anatomy was adapted by **Brian Pridgen, MD**. His original viewer
software is not included in Human Atlas. Human Atlas application code is
distributed under [GPLv3](LICENSE); the anatomy and third-party components
below retain their separate licenses.

### Anatomical geometry and anatomy-derived metadata

This free educational release includes separately identifiable components with
**different licenses**. It must not be described as uniformly CC BY-SA.

- The previously released adapted anatomy and catalogue remain under
  [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/), including
  Brian Pridgen's copyrightable contributions to those adaptations.
- Both kidneys by **Lissie Cowley** are under
  [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).
- Both cochleae and vestibules from **University of Dundee School of Medicine**
  are under [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
  Shared adaptations of these components retain that license.

The kidney and inner-ear components may be used only for noncommercial purposes
unless the relevant rights holders grant additional permission. Their inclusion
here does not grant permission to use them in commercial products or promotional
material. Other components retain their own permissions, including commercial
reuse where their licenses permit it. Viewer software terms do not restrict the
rights granted by these component licenses.

#### Open 3D Model

Open 3D Model — departments of anatomy of Leiden University Medical Center,
University Medical Center Utrecht, Maastricht University Medical Center and
University of Leuven KULAK, with Amsterdam UMC, Radboud University Medical Center
Nijmegen and University of Ghent. Funded by the Dutch Ministry of Education,
Culture and Science.

- [Upper limb, English labels](https://anatomytool.org/content/open3dmodel-upper-limb-english-labels)
- [Brachial plexus and branches, English labels](https://anatomytool.org/content/open3dmodel-brachial-plexus-and-branches-english-labels)
- Reviewed input edition: July 2025; CC BY-SA 4.0.

#### Z-Anatomy and upstream contributors

[Z-Anatomy — The libre 3D atlas of anatomy](https://github.com/Z-Anatomy/Models-of-human-anatomy),
CC BY-SA 4.0. The reviewed source notice is available at
[the pinned upstream license](https://github.com/Z-Anatomy/Models-of-human-anatomy/blob/bb293be456f6d2245191f37752d14ec76894740d/License.txt).
Its upstream attributions are retained here:

- **BodyParts3D — The Database Center for Life Science — CC-BY-SA 2.1 Japan**, as
  credited by Z-Anatomy. [Original model](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html);
  [license](https://creativecommons.org/licenses/by-sa/2.1/jp/).
  Kousaku OKUBO is credited for the original BodyParts3D model.
- Gauthier KERVYN: design, 3D and anatomy.
- “Brainder” and “White matter” from the University of Washington, identified
  by Z-Anatomy as reference/included and adapted material. Its notice does not
  specify separate license terms for these entries; this is not a claim of
  independent clearance of those underlying sources.
- Cranial Nerves and Foramina — University of Dundee, CAHID —
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- The upstream project also credits Marcin ZIELINSKI (Blender add-on), Lluis
  VINENT (Unity development), Ana Teresa BIGIO (Portuguese translation), and
  Carlos TORRES VILLAR (Spanish translation). Their application software is
  not bundled in this viewer.

BodyParts3D: Mitsuhashi N, Fujieda K, Tamura T, Kawamoto S, Takagi T, Okubo K.
BodyParts3D: 3D structure database for anatomical concepts. Nucleic Acids Res
2009;37:D782–5. [PMID 18835852](https://pubmed.ncbi.nlm.nih.gov/18835852/).
The upstream BodyParts3D credit above concerns model ancestry; it does not
assert that a historical measurement dataset is included.

#### Changes made for this viewer

Adaptation and integration by Brian Pridgen. Released adaptations include
combining the whole-body and upper-limb sources, mirroring the grafted upper
limb for the opposite side, preparing regional geometry chunks and a derived
body surface, and adapting names, catalogue metadata and display organization.
The viewer also crops regional skin and adjusts display visibility at runtime.
These are adaptations, not unmodified upstream releases. Source and mirroring
information are recorded per structure. Neither a mirrored structure nor a
display layer establishes independent anatomical validation.

#### Educational components restored in this release

The six kidney and inner-ear meshes were recovered from the recorded Z-Anatomy
source edition, preserving their source placement. Export changes comprise
applying source object transforms/modifiers, coordinate conversion to glTF,
and assigning separate named objects and license metadata. They are not fused
with other meshes. Existing released mesh geometry was retained unchanged.
Model conversion and integration by Brian Pridgen; the component licenses above
continue to apply. This restoration does not establish anatomical validation.

#### Material not included

Historical BodyParts3D-derived calibre measurements
and literature measurement tables are not bundled. Z-Anatomy's Wikipedia
text definitions are not supplied as a definitions dataset in this release.

### Rendering libraries and fonts

- [three.js](https://github.com/mrdoob/three.js), r160: MIT; copyright 2010–2023 Three.js Authors.
- [meshoptimizer](https://github.com/zeux/meshoptimizer): MIT; copyright 2016–2026 Arseny Kapoulkine.
- [IBM Plex](https://github.com/IBM/plex): SIL Open Font License 1.1; copyright 2017 IBM Corp., reserved font name “Plex”.

Dependency permission, copyright, warranty, and font notices are provided
with their respective packages.
The Open 3D Model GPL-3.0 viewer and Z-Anatomy application code are not used.

### Attribution when sharing anatomy

Retain the relevant upstream credits, license links and modification notices,
and identify your further changes. A reasonable credit for the combined
anatomy is “Anatomy adapted by Brian Pridgen from Open 3D Model and Z-Anatomy,
with BodyParts3D and other upstream contributors; component licenses include
CC BY-SA 4.0, CC BY-NC 4.0 (kidneys, Lissie Cowley), and CC BY-NC-SA 4.0
(inner ear, University of Dundee School of Medicine); see credits.” Include this document or a link to its hosted
copy. Videos, screenshots, downloads and other shared outputs containing the kidney
or inner-ear models must retain their attribution and applicable license
notices and comply with their NonCommercial terms.
This example does not replace the applicable license terms.

### Intended use and warranties

Educational and reference use; the viewer is not intended for diagnosis,
clinical decision-making or intraoperative guidance. Anatomy may be incomplete
or inaccurate. No upstream institution or contributor is represented as
endorsing this viewer. Material is provided as is, subject to the warranty
and liability provisions of its applicable license. This intended-use notice
adds no restriction to rights granted by the anatomy licenses.
