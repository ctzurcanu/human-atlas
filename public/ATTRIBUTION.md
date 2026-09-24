# Human Atlas attribution and component licenses

The application code is distributed under [GNU GPLv3](LICENSE.txt). Original upstream application code is copyright © 2026 ashemag; its [MIT notice](MIT-upstream.txt) remains available. Anatomical datasets and other components retain their own licenses.

## Detailed male: Z-Anatomy

Source: [Z-Anatomy original Blender archive](https://github.com/Z-Anatomy/Models-of-human-anatomy), CC BY-SA 4.0, with contributors identified in the [publisher license](https://github.com/Z-Anatomy/Models-of-human-anatomy/blob/master/License.txt). Source ZIP SHA-256: `e029688545627bd0214b269e1063143abb580aad72b2c2445d6d8a9a0d9da736`; `Startup.blend` SHA-256: `9f08a17ea0115fed80b2a73ecdf0a1bc2ab2f6956f37c593ce23d513ea35afcd`.

Two kidneys by Lissie Cowley retain CC BY-NC 4.0. Two cochleae and two vestibules by the University of Dundee School of Medicine retain CC BY-NC-SA 4.0. The source also credits University of Dundee CAHID cranial nerve material under CC BY 4.0. These component licenses are recorded on the corresponding mesh entries.

Adaptation by Human Atlas: original Blender objects, collection membership, materials, and world transforms are read directly; coordinates are converted to Y-up; material-defined surfaces are packed without triangle reduction; source material colors are retained alongside display colors. Reference lettering and diagrams are excluded. The original archive has no anatomical image texture files. This source contributes 3,817 selectable surfaces and 3,213 concepts to the detailed male build.

## Detailed male supplement: Open 3D Model upper limb

Source: [Open 3D Model publisher source files](https://anatomytool.org/open3dmodel-create), from Leiden UMC and its [credited project partners](https://anatomytool.org/open3dmodel-credit), CC BY-SA 4.0. The original [upper-limb GLB archive](https://caskanatomy.info/open3dmodelfiles/upper-limb/upper-limb-glb.zip) has SHA-256 `5af0190a6d7bf47393447ac30021e4f3ba619721c7f3a620c39a895947078432`; its GLB has SHA-256 `e440c84c794239d1850e62b4ede0195d81bcff4f0078c7945528c388ab72fdb4`.

Adaptation by Human Atlas: Blender decodes the original GLB directly. Meshes with source names already represented by Z-Anatomy are omitted. Additional right-side arm, forearm, hand, pectoral, and brachial neurovascular meshes are converted into the same coordinates. Each additional right-side mesh is mirrored to the left, with that derivation recorded on the mesh. The result adds 588 selectable surfaces and 588 concepts, bringing the combined detailed male atlas to 4,405 surfaces and 3,801 concepts. The source GLB embeds 84 images; the current viewer uses source material colors in its common display palette and does not render those images.

## Female: Human Reference Atlas and Visible Human Female

Source: [Human Reference Atlas united-female v1.10](https://purl.humanatlas.io/ref-organ/united-female/v1.10), by Kristen Browne and Heidi Schlehlein / HuBMAP, CC BY 4.0. Original GLB SHA-256: `95f0c3d2f918582608692ca1139e8bdb18c147a16470e9ee9af8b276bd77c422`.

Additional source: [Human Reference Atlas united-female v1.5](https://purl.humanatlas.io/ref-organ/united-female/v1.5), same authors and license. Original GLB SHA-256: `472567a56896b9b7890508da6501fbf858e56aaa30745365f7a71ade782b529c`. Eight pubis and ischium bone surfaces in v1.5 are absent from v1.10.

Additional source: [University of Denver Visible Human Female 3D STL models](https://digitalcommons.du.edu/visiblehuman/1/), CC BY 4.0. Original `Final 3D STL Models-stl.zip` SHA-256: `9886eda040f6087bbb65182530f8bd262be456f366b155025c12b3dc56da774e`. The source archive remains in the ignored local source folder and is not distributed with the site.

Adaptation by Human Atlas: HRA source nodes, labels, ontology identifiers, transforms, and material colors are read directly. HRA mesh geometry is simplified with a 0.2% relative error bound and packed for the browser; neither HRA GLB has image textures. The original Denver archive contains 76 lower-limb muscle STLs. We independently align them to HRA femur, tibia, and fibula geometry and include 74; the two rectus femoris STLs duplicate existing HRA meshes. Per-bone fit measurements and transforms are in the model manifest. Female anatomy has 1,030 selectable meshes and 1,234 concepts. The separate Embryo model contains the eight original HRA placenta-group meshes: basal plate, chorionic plate, placental vessels, amnion, umbilical cord, two umbilical arteries, and umbilical vein. The HRA source does not include an embryo or fetal body. Female upper-body muscle and skeleton coverage remains partial, and the separate source sets should not be treated as a single-person scan.

## Standard male: BodyParts3D

BodyParts3D, © The Database Center for Life Science, [CC BY 4.0](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/lic.html). Source: [BodyParts3D 4.0 download](https://dbarchive.biosciencedbc.jp/en/bodyparts3d/download.html), `isa_BP3D_4.0_obj_99.zip` and its official IS-A and PART-OF tables. See Mitsuhashi et al. (2009), [BodyParts3D: 3D structure database for anatomical concepts](https://doi.org/10.1093/nar/gkn613).

Adaptation by Human Atlas: source millimeter/Z-up coordinates are converted to meter/Y-up; meshes are simplified with a 0.2% relative error bound; normals are quantized; mesh data are packed in binary chunks; display systems and colors are curated. All 2,234 source OBJ meshes are represented, with 3,432 named FMA concepts. Source identity and original file hashes are recorded in the manifest. Some OBJ comments state an older license; the publisher's current license page above applies to the current dataset.

## Human eukaryotic cell

Source: [Human Eukaryotic Cell and its main components](https://sketchfab.com/3d-models/human-eukaryotic-cell-and-its-main-components-e769645f40fc47c38b89bcfc6b8b4d44) by [E-learning UMCG](https://sketchfab.com/eLearningUMCG), [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/). Source GLB SHA-256: `826a9b440515176cfa1a5e8d567d73a215dd3ce48c5602a35fdc18dcb18dc2ea`.

Adaptation by Human Atlas: the 20 original named meshes, their vertex colors, and source transparency are converted to selectable cell components with function and compartment groupings. The source GLB and adapted cell assets retain CC BY-NC-SA 4.0, separate from the application code license.

## Libraries and fonts

- [three.js](https://github.com/mrdoob/three.js): MIT, Three.js Authors.
- [meshoptimizer](https://github.com/zeux/meshoptimizer): MIT, Arseny Kapoulkine.
- [IBM Plex](https://github.com/IBM/plex): SIL Open Font License 1.1, IBM Corp.

The viewer is for education and reference. It is not intended for clinical decisions. Preserve the applicable source credits, component licenses, and modification notices when sharing anatomy or derived images.
