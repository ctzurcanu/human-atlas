# TA98 tongue inventory

The table records the original source inventory. Local construction stages below are not installation or anatomical acceptance.

The authority contains **36 relevant terms**, including the separate muscle group, lingual tonsil and palatoglossus cross-references. Their canonical parents are preserved. The thyroglossal duct is a variant and is not required in the ordinary adult assembly.

Both models have source evidence for only **5 terms**; **31 terms per model have no matching source geometry**. Source evidence includes container membership and hidden, misregistered sources, so it is not a completion count.

| TA98 | Term | Parent | Representation | Male source meshes | Female source meshes |
|---|---|---|---|---|---|
| A04.1.05.001 | muscles of tongue | A04.1.00.001 | whole-organ | `LOCAL:male:l_styloglossus:1030` (hidden); `LOCAL:male:r_styloglossus:981` (hidden); `ZA:Genioglossus muscle.l`; `ZA:Genioglossus muscle.r`; `ZA:Hyoglossus muscle.l`; `ZA:Hyoglossus muscle.r` | **missing** |
| A05.1.04.001 | tongue | A05.1.00.001 | whole-organ | `ZA:Tongue` | `LOCAL:female:tongue:1089` |
| A05.1.04.002 | body of tongue | A05.1.04.001 | tissue-region | **missing** | **missing** |
| A05.1.04.003 | root of tongue | A05.1.04.001 | tissue-region | **missing** | **missing** |
| A05.1.04.004 | dorsum of tongue | A05.1.04.001 | surface-region | **missing** | `HRA:VH_F_dorsal_tongue` (hidden) |
| A05.1.04.005 | anterior part | A05.1.04.004 | surface-region | **missing** | **missing** |
| A05.1.04.006 | posterior part | A05.1.04.004 | surface-region | **missing** | **missing** |
| A05.1.04.007 | inferior surface of tongue | A05.1.04.001 | surface-region | **missing** | **missing** |
| A05.1.04.008 | fimbriated fold | A05.1.04.007 | tissue-region | **missing** | **missing** |
| A05.1.04.009 | margin of tongue | A05.1.04.001 | boundary-curve | **missing** | **missing** |
| A05.1.04.010 | apex of tongue | A05.1.04.001 | point-landmark | **missing** | **missing** |
| A05.1.04.011 | mucous membrane of tongue | A05.1.04.001 | wall-layer | **missing** | **missing** |
| A05.1.04.012 | frenulum of tongue | A05.1.04.001 | tissue-region | **missing** | **missing** |
| A05.1.04.013 | papillae of tongue | A05.1.04.001 | microscopic-structure | **missing** | **missing** |
| A05.1.04.014 | filiform papillae | A05.1.04.013 | microscopic-structure | **missing** | **missing** |
| A05.1.04.015 | fungiform papillae | A05.1.04.013 | microscopic-structure | **missing** | **missing** |
| A05.1.04.016 | vallate papillae | A05.1.04.013 | microscopic-structure | **missing** | **missing** |
| A05.1.04.017 | foliate papillae | A05.1.04.013 | microscopic-structure | **missing** | **missing** |
| A05.1.04.018 | midline groove of tongue | A05.1.04.001 | boundary-curve | **missing** | **missing** |
| A05.1.04.019 | terminal sulcus of tongue | A05.1.04.001 | boundary-curve | **missing** | **missing** |
| A05.1.04.020 | foramen caecum of tongue | A05.1.04.001 | opening | **missing** | **missing** |
| A05.1.04.021 | thyroglossal duct | A05.1.04.020 | tissue-region | **missing** | **missing** |
| A05.1.04.022 | lingual tonsil | A13.2.02.001 | tissue-region | **missing** | **missing** |
| A05.1.04.023 | lymphoid nodules | A05.1.04.022 | microscopic-structure | **missing** | **missing** |
| A05.1.04.024 | lingual septum | A05.1.04.001 | tissue-region | **missing** | **missing** |
| A05.1.04.025 | lingual aponeurosis | A05.1.04.001 | tissue-region | **missing** | **missing** |
| A05.1.04.101 | genioglossus | A04.1.05.001 | tissue-region | `ZA:Genioglossus muscle.l`; `ZA:Genioglossus muscle.r` | `LOCAL:female:genioglossus:239` |
| A05.1.04.102 | hyoglossus | A04.1.05.001 | tissue-region | `ZA:Hyoglossus muscle.l`; `ZA:Hyoglossus muscle.r` | `LOCAL:female:l_hyoglossus:242`; `LOCAL:female:r_hyoglossus:255` |
| A05.1.04.103 | chondroglossus | A05.1.04.102 | tissue-region | **missing** | **missing** |
| A05.1.04.104 | ceratoglossus | A05.1.04.102 | tissue-region | **missing** | **missing** |
| A05.1.04.105 | styloglossus | A04.1.05.001 | tissue-region | `LOCAL:male:l_styloglossus:1030` (hidden); `LOCAL:male:r_styloglossus:981` (hidden) | `LOCAL:female:l_styloglossus:243`; `LOCAL:female:r_styloglossus:256` |
| A05.1.04.106 | superior longitudinal muscle | A04.1.05.001 | tissue-region | **missing** | **missing** |
| A05.1.04.107 | inferior longitudinal muscle | A04.1.05.001 | tissue-region | **missing** | **missing** |
| A05.1.04.108 | transverse muscle | A04.1.05.001 | tissue-region | **missing** | **missing** |
| A05.1.04.109 | vertical muscle | A04.1.05.001 | tissue-region | **missing** | **missing** |
| A05.1.04.110 | palatoglossus | A04.1.06.001 | tissue-region | **missing** | **missing** |

## Current source limitations

- Male: the active whole tongue has 343 vertices / 682 triangles. Paired genioglossus and hyoglossus are present. The two native styloglossus sources are hidden and explicitly await alignment/overlap review; do not simply unhide them.
- Female: the native whole tongue has 205 vertices / 304 triangles. Genioglossus and paired hyoglossus/styloglossus are present. The hidden HRA dorsal source is large and misregistered; it does not establish an installed dorsum.
- The existing oral-muscle plan explicitly leaves intrinsic muscles, palatoglossus, chondroglossus/ceratoglossus boundaries, attachments and fibre architecture unverified or missing.
- Surface regions use owned surface subsets, sulci/margins use curves, the apex uses a point, and the foramen caecum is a blind surface pit. None should become an unrelated solid.
- Mucosa needs finite tissue ownership. Septum/aponeurosis are connective tissue. Papilla types need correct scale and distribution, preferably grouped relief rather than a mesh for every repeated element.
- The male assembly is completed first with the skeleton fixed, then transferred and fitted to female mandibular, hyoid, styloid and palatal attachments. Existing source muscle envelopes alone do not establish those attachments.
- Keep the body assembly economical. A separate detailed organ model can own additional meshes and hierarchies, with its count excluded from the body model.

Full source records, catalogue/binary hashes and individual missing/unverified flags are in [reports/ta98-tongue-inventory.json](reports/ta98-tongue-inventory.json). Reproduce with the bounded preparation script:

```sh
python3 -B -m geometry_tools.runtime --seconds 20 -- python3 -B server/prepare-tongue-inventory.py
```

## Current local construction

The intrinsic proposal has four compact macroscopic muscle-region carriers: one superior longitudinal, two inferior longitudinal and one shared transverse/vertical territory. Transverse and vertical retain separate TA98 controls and directional fibre materials; the shared mesh is counted once. These are approximate muscle envelopes, not exclusive histological tissue cells. The latest side view still requires crease-normal and regional-boundary corrections. The attempted bank preparation stopped at the bone gate and did not install anything.

The complete bone audit found **280 hyoid crossings in the original whole-tongue source**, explaining the basal crossings in the first derived intrinsic and mucosal candidates. A supported, globally invertible protected cage now lifts only the basal region by at most 6 mm: actual stored geometry has zero self/bone crossings, all 145 protected source vertices and all 682 face indices remain literal, and minimum measured hyoid clearance is 1.58 mm. The corrected envelope is 41.547 mL versus the original 45.383 mL. Fixed bone surfaces were not capped, repaired or moved. Actual right, dorsal and sagittal Blender views were inspected; the result is accepted as a macroscopic construction carrier, not as the completed tongue. Reproducible source-bound evidence is in the local `basal-whole-protected-v1` stage.

The paired genioglossus candidate removes the original two self-crossings per side, closes the mandibular attachment gap and removes hyoid penetration. Actual stored Float32 geometry has zero strict self/jaw/hyoid crossings. Bony contact is explicitly adjacent-unjoined with a measured 0.497–0.536 µm gap and a 1 µm numerical tolerance; this does not establish microscopic enthesis. The 1,261 upper source faces and their corner normals are literal. Intralingual placement must be checked after the whole-envelope correction.

Four intrinsic carriers have been regenerated from the corrected whole: **9,106 triangles**, with all six intrinsic-pair checks and all 32 fixed-bone checks complete and crossing-free. Separate corner-normal charts retain their geometry and address crease shading; the regenerated full assembly still requires visual/interface review. Six typed surface/point/curve controls have also been regenerated on the corrected source, with runtime bank binding deliberately pending rather than falsely pointing at the old geometry bank.

Finite mucosa and dorsal lingual aponeurosis cells share one reciprocal interface bank, totaling **34,512 bytes without image textures**. After regeneration from the corrected whole, individual cells, their mutual occupied-volume checks and all 16 bone pairs pass. However, the full closed mucosa coats computational basal boundaries crossed by the extrinsic muscle envelopes, and the aponeurosis still crosses the superior-longitudinal carrier. These are unresolved assembly defects, not permitted tissue interdigitation: define the exposed mucosal extent and entry boundaries, fit intralingual envelopes and establish the aponeurotic interface before installation. Consequently these cells are **uninstalled, incomplete assembly candidates**, not completed TA98 layers. Thickness fields are modeling estimates informed by [May 2022 human histology](https://onlinelibrary.wiley.com/doi/10.1111/joa.13637), not measurements of this source tongue. Papillae, septum, body/root divisions, remaining extrinsic subdivisions and their interfaces still need construction/review.

The subsequent exposed-mucosa candidate replaces that incorrect full basal coating with one connected **348-triangle carrier patch**, finite thickness and one reflection rim. An exposed-side intralingual correction removes all measured mucosal crossings for the four GG/HG envelopes while retaining actual below-reflection source coordinates and corner normals. Native HG hyoid-origin defects remain separate unresolved work. Six typed controls now reference the real exposed patch: the inferior surface means the free ventral body/blade, excluding the computational basal stem. Every cut-facet ancestor is recorded; no extra solids are added for these controls.

The aponeurosis–superior-longitudinal complement now has **122 literal reciprocal facets** and zero actual stored pair crossings. One bounded common-boundary triangulation repair fixes Float32 subdivision defects without moving source controls or changing thickness. This passed geometric interface still requires common storage with the exposed mucosa. Current actual Blender surfaces/sections show coarse angular intrinsic envelope features, so numerical success is not complete anatomical or appearance acceptance.

An original paired **palatoglossus** candidate has **600 triangles and 14,832 estimated packed bytes per side**, without images or UVs. Its flattened belly, palatal fan and lingual taper follow the [adult cadaver investigation](https://pubmed.ncbi.nlm.nih.gov/281278/); dimensions are explicit modeling estimates. Native-palate ray witnesses locate the full fan inside the actual finite palate, and each distal cap is inside the corrected tongue. Independent whole-cap, single connected attachment-footprint, self, normal and fixed-context checks pass. Actual front/right Blender surfaces and sagittal section were inspected. Faucial-pillar mucosal continuity remains pending; the candidate is **not installed or accepted as a completed tongue**. Local stages and inspected-image hashes retain these qualifications.

The emitted `coupled-male-v3` assembly was rejected after actual Blender surface and section review: its repaired HG footprint produced broad planar bands, and intrinsic Boolean regions remained angular. Its passing numerical installation job does not authorize applying that rejected geometry. The replacement `native-hg-short-origin-v1-bilateral` preserves the original curved bellies, uses short cortical transitions, and removes the false long posterior footprint. Both actual stored muscles have zero tested self crossings; all sixteen bone pairs, four covering pairs and the bilateral pair pass. Four actual Blender views were inspected and accepted for the next coupled candidate, with attachment territory and histological limits retained. Neither ordinary atlas has changed.

The source-owned papillary detail field now has real filiform and fungiform relief and three additional typed controls on one mucosal host. It is one small anterior field, not full dorsal coverage; its detailed geometry remains separate from the economical body candidate. The PG pair packs into **30,512 bytes / 1,200 triangles**, but its actual lingual entry crosses the finite mucosa in **52 face pairs per side**. The addon explicitly reports `assembly_ready: false`; a cap inside a whole-organ envelope does not resolve that covering interface. The separate palatoglossal arch must retain its own A05.2.01.005 identity rather than being relabeled as tongue mucosa. Current saved geometry, renders, witnesses and remaining limits are in the local tongue workbench.

The partial male revision `coupled-male-v4-curved` is now installed locally with **17 of 35 ordinary-adult TA98 controls on 12 unique physical meshes**. The curved inferior-longitudinal and shared transverse/vertical reference envelopes replace the angular candidates; repaired native-based HG bellies replace the false posterior fans. Actual emitted surfaces/sections, all 115 installation constraints, served JSON/bank hashes and systems/TA98 control checks pass. The two tongue banks total **1,208,449 bytes**, with no added image textures. Skeleton and unrelated organs remain literal; the female source is unchanged. This is not a completed tongue: 18 adult terms, full papillary distribution, precise muscle territories/entheses, septum, body/root landmarks, palatoglossal covering and the floor-of-mouth/frenulum remain unfinished. Failed arch and FOM candidates are excluded. The local rollback is `male-tongue-v4-curved-partial/manifest.json`; subsequent work must rebind against source SHA `1cdb70649f10bdb523beb38f1c0d09167e2344a02c3d8ce74d7f40739a777d39` rather than overwrite it from the old source inventory.

The female core revision is now installed locally with **14 of 35 ordinary-adult controls on six physical meshes**, at 0.846582 of the male linear size with a protected lower-region fit for the inferior longitudinal pair. Mucosa, aponeurosis and four intrinsic carriers pass actual encoded self/normal, fixed-context and reciprocal-interface checks; 92 correctly classified constraints, native Blender surfaces/section, served geometry banks and Systems/TA98 controls were verified. Five original female extrinsic muscles remain retained but suppressed because their transferred/native fits crossed the covering; their three terms are excluded from the new assembly. **21 adult terms remain unfinished**, including those attachments. Skeleton, male anatomy and existing visit coordinates remain literal; the pre-existing stale visit is not certified or refreshed. The rollback is `female-tongue-core-partial-v1/manifest.json`, and current female source SHA is `8990a73098b924aa51bd4f23a0f5e52ccd2518fbb8365740f954aebe85f21a4e`. This is a partial macroscopic preview with coarse exterior/reference envelopes, not a completed tongue or accepted histology.
