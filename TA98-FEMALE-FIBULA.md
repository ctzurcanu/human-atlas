# Female fibula: full TA98 construction inventory

Both original native female fibulae are staged without coordinate, topology, normal, UV or material changes. The source master SHA-256 is `7e17153e8d56e59407a72a972e55ac770f59f48b67a4f328f5e0e0b739e98106`. Each side contains 1,168 source triangles. This hash identifies the preserved pre-installation checkpoint. The independently registered HRA bones are excluded from the native assembly.

Authority: [IFAA TA98, skeletal system](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/02%20TA98%20EN.htm). All seventeen terms below were inventoried before any proposed partition. Four regional exteriors per side are now installed in the local review model; anatomical boundary acceptance and twelve fine-feature representations remain incomplete.

| TA98 | Term | Parent | Representation | Native independent mesh |
| --- | --- | --- | --- | --- |
| A02.5.07.001 | fibula | A02.5.03.001 | whole bone exterior | left/right whole source only |
| A02.5.07.002 | head | A02.5.07.001 | material bone region | installed bilateral exterior; provisional boundary |
| A02.5.07.003 | articular facet | A02.5.07.002 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.004 | apex | A02.5.07.002 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.005 | neck | A02.5.07.001 | material bone region | installed bilateral exterior; provisional boundary |
| A02.5.07.006 | shaft | A02.5.07.001 | material bone region | installed bilateral exterior; provisional boundary |
| A02.5.07.007 | lateral surface | A02.5.07.006 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.008 | medial surface | A02.5.07.006 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.009 | posterior surface | A02.5.07.006 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.010 | medial crest | A02.5.07.009 | boundary curve | missing; owning whole bone available |
| A02.5.07.011 | anterior border | A02.5.07.006 | boundary curve | missing; owning whole bone available |
| A02.5.07.012 | interosseous border | A02.5.07.006 | boundary curve | missing; owning whole bone available |
| A02.5.07.013 | posterior border | A02.5.07.006 | boundary curve | missing; owning whole bone available |
| A02.5.07.014 | lateral malleolus | A02.5.07.001 | material bone region | installed bilateral exterior; provisional boundary |
| A02.5.07.015 | articular facet | A02.5.07.014 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.016 | malleolar fossa | A02.5.07.014 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |
| A02.5.07.017 | malleolar groove | A02.5.07.014 | owning-bone surface feature; not a separate solid | missing; owning whole bone available |

Native sources: `LOCAL:female:l_fibula:670` and `LOCAL:female:r_fibula:666`. All planned regional boundaries belong to these source surfaces. No male cut level or donor dimension is transferred.

The planned tissue assembly consists of head, neck, shaft and lateral malleolus. Articular facets must follow reviewed native surfaces; apex, fossa and groove must be distinguished from filled solids. Surface borders and crests require explicit feature representations. A whole exterior does not establish cortical/trabecular tissue, marrow or periosteal layers.

Adjoining interfaces to review:


- head articular facet / superior tibiofibular joint
- head / neck / shaft shared source tissue
- shaft / lateral malleolus shared source tissue
- distal talar facet / talus and inferior tibiofibular syndesmosis
- interosseous border / membrane of leg
- head / biceps femoris tendon and fibular collateral ligament
- neck / common fibular nerve
- posterior malleolar groove / fibular tendons

Whole-source evidence: exact-coordinate topology shows one closed genus-zero component per side, zero open/nonmanifold/winding/duplicate/degenerate/vertex-link defects. All triangle AABB overlaps were checked with float64 narrow-phase tests at 1e-7 m tolerance; zero proper/coplanar crossings were detected. Thirty-six actual sections are closed. Sections through the curved long bone can form multiple closed contours; these are not filled or merged. Source positions, indices, packed normals, UVs and materials match independently, including 52 native neighboring guide records.

Both actual textured source bones were inspected, with two oblique directions on the left, and the complete section montage was inspected. These passes do not certify compartment packing, attachment footprints, fine feature identity or anatomical boundary accuracy. Whole source repair is not indicated by the measured topology; native landmark and adjoining-interface review must precede regional cuts. All seventeen terms remain anatomically unaccepted.

The live review also exposed a missing arterial synonym: fibular artery equals peroneal artery. The relationship resolver now recognizes both while preserving side, tissue-system and native-frame guards. The source fibular shaft receives nutrient/periosteal branches from that artery; proximal contribution is distinguished from anterior tibial branches. [NCBI fibula reference](https://www.ncbi.nlm.nih.gov/books/NBK470591/?report=printable). No small nutrient branch or nutrient foramen was fabricated.

Reproduce source evidence with `scripts/prepare-ta98-female-fibula.py` then `scripts/validate-ta98-female-fibula-source.py`. Reports: `reports/ta98-female/fibula-{parts-plan,source-review,source-preservation}.json`; sections: `outputs/ta98-female-fibula-sections.png`. Relationship regression: `node --test tests/fibular-supply.test.mjs`. Local review sources remain local, with redistribution rights unverified.

## Native landmark measurements and related nerves

`review-ta98-female-fibula-landmarks.py` measures 180 closed transverse sections per bone. The actual area-profile rendering was inspected. Both native exteriors extend from Y 0.056510 to 0.487120 m; midshaft median section area is approximately 93.83 mm². The distal enlargement peaks near Y 0.0783 m and the proximal enlargement near Y 0.4725 m. These are measured source landmarks, not automatically inferred anatomical cuts. Subsequent native adjoining-interface review constrains the provisional regional candidate below; anatomical boundary acceptance remains pending. [Gray plate 258](https://commons.wikimedia.org/wiki/File:Gray258.png) was visually inspected as an anatomical reference; its donor dimensions are not transferred.

The long muscular source records labeled calcaneofibular ligament span Y 0.044653–0.310753 m and approach the fibula at Y 0.159062 m. This extent is inconsistent with an ankle ligament; the landmark report explicitly excludes them from boundary decisions. Other nearby guides remain provisional too: the measurements use nearest original vertices, not continuous surface contact or attachment footprints. No mesh was moved or relabeled to force agreement.

Both native common peroneal nerves exist (`LOCAL:female:r_common_peroneal_nerve:817` and `LOCAL:female:l_common_peroneal_nerve:821`). The relationship resolver now recognizes common, deep and superficial peroneal as the source synonyms of the corresponding fibular nerves. This restores the reciprocal sciatic → common → deep/superficial routes on each side, while excluding contrary-side, suppressed, wrong-system and independently registered donor records. Branching nomenclature reference: [NCBI Bookshelf NBK532968](https://www.ncbi.nlm.nih.gov/books/NBK532968/). This semantic branching relationship does not certify the actual bifurcation surfaces or nerve course against the neck.

The actual left common nerve inspector was reviewed and shows both terminal branches. Its existing coarse native surface remains anatomically unaccepted. No geometry, texture, material group or network asset was added. Four targeted fibular regressions, TypeScript, production build and complete male/female directed-relation reciprocity pass (4,026 female links; 16,084 male links). Source fibular surfaces remain unchanged. Whole-body TA98 completeness, regional boundaries and packing are still incomplete.

## Local regional exterior candidate

The native source head/neck and ankle were inspected anteriorly and posteriorly with tibia, talus, head attachments and the original peroneal nerve surfaces. `stage-ta98-female-fibula-interfaces.py` stages six exact original nerve records without changing the whole-source checkpoint. The local page `web/female-fibula-review.html` displays the source/candidate, original textures, optional review colors and native context; it is a review tool, not required by the main model.

`construct-ta98-female-fibula-regions.py` creates eight regional exterior meshes: head, neck, shaft and lateral malleolus on each side. Requested Y interfaces are 0.112, 0.448 and 0.462 m, with at most 30 µm numerical adjustment to avoid source vertices. The proximal taper and head enlargement constrain the neck; distal enlargement and native talar/syndesmotic extent constrain the malleolus. These planar boundaries are explicitly provisional. No male proximal cuts or donor dimensions were used.

The candidate preserves every original exterior face exactly once, retaining native coordinates, UVs, normal interpolation and materials, without overlap or added caps. Original source triangles total 2,336; conforming seam subdivisions produce 2,748 candidate triangles (+412). Existing textures are reused. The original candidate remains uncompressed staging data. Its compact derivative is now installed locally, as described below; original exterior textures and materials are reused.

`validate-ta98-fibula-regions.py --female` independently checks full source-face ownership, fragment overlap, barycentric positions and attribute interpolation. Maximum measured surface error is below 3.2e-17 m; UV error is below 1.2e-16. The complete assemblies remain closed genus zero at double and float32 precision, with zero measured proper/coplanar self-crossings at the existing 1e-7 m tolerance. The same shared preservation verifier still passes the original male candidate.

Actual candidate head/neck and distal renders on both sides were inspected, with original textures and region review colors. Thirty-six whole cuts and twelve cuts straddling the seams were inspected; contours remain closed, including the inherited separate contours of the curved whole bone. Native attachment packing, fine feature identity and planar boundary accuracy remain unaccepted. The original staging candidate stays local; its compact derivative is installed in the local review catalogue. All twelve remaining facet/apex/surface/border/crest/fossa/groove terms remain explicitly missing, and these four regional exteriors do not establish separate internal tissue volumes or complete anatomical coverage.

Candidate: `.local-models/ta98-review/female-fibula-regions-candidate.json`. Evidence: `reports/ta98-female/fibula-regions-{construction,preservation,seam-sections}.json`. Cross-section renderer: `scripts/render-ta98-female-fibula-region-sections.py`. Continue with reviewed surface-feature ownership while retaining missing features and provisional boundary status.

## Compact buffers and whole-assembly relationship gate

`pack-ta98-female-fibula-regions.py` independently decodes its compact buffers and checks the candidate positions, remapped indices, quantized normals and UVs. All eight region records total 78,272 bytes, or 39,353 bytes gzipped. Compaction reduces unused staging vertices from 6,120 to 1,742 across the eight arrays; this is not decimation of the anatomy. Both packed whole assemblies remain closed genus zero. The manifest provides native-frame record IDs, exact term identities with explicit provisional coverage notices, whole-source replacement sets and source assembly identities. This packing step did not change the main catalogue. The subsequent local installation is documented below; no buffer is deployed publicly.

Integration review found that whole-bone functional queries previously consulted only the visceral supply table. A fully assembled fibula therefore lost arterial supply in the inspector after regional subdivision. The resolver now supports the stored complete native fibula assembly: all selected records must belong to its exact native assembly, and the selection must cover every active assembly record. Partial selections, forged partial whole selections, selection-set queries, opposite-side vessels and independently registered donors cannot inherit these links.

The actual male left whole-fibula inspector now shows the fibular artery and anterior tibial artery, with notes distinguishing shaft nutrient/periosteal supply from proximal head contribution. Both actual male assemblies pass the whole-query regression; female regional assembly fixtures pass the same safeguards. These whole-assembly semantic queries do not grant every individual regional surface the same supply or certify mesh contacts. Six targeted fibular regressions, TypeScript, requested-selection identity across 10,221 nonempty concepts, existing gastric/pelvic functional regressions and production build pass.

Packed manifest: `.local-models/ta98-review/female-fibula-regions-packed.json`; evidence: `reports/ta98-female/fibula-regions-packed.json`. Local installation now preserves concept replacements, old source-piece links, explicit region-to-whole parents and region-specific attachments. Actual main-view and component review is documented below. Fine features and anatomical boundary acceptance remain incomplete.

## Installed local runtime and main-view review

`install-ta98-female-fibula-regions.py` installs the eight compact exterior records and suppresses the original two whole meshes, retaining their provenance. Every existing concept membership is remapped; legacy numeric source-piece URLs resolve to the four-piece whole assembly. Original materials and texture files stay unchanged. The backup `.local-models/ta98-review/female-before-fibula-regions.json` retains the source checkpoint; source preparation/packing scripts deliberately reject the changed master and must not be rerun against it as if it were the original source.

`validate-ta98-female-fibula-export.py` independently verifies actual packed attributes, complete assembled topology, all unrelated part/concept/material records and gzip replay. Installed master SHA-256: `e2fb4d326abdbd0b9c2cc34d6e0e867daf33f869f46fa1243cd1306886ac67d1`. The runtime contains 2,780 active records; the eight regional records replace two originals, adding 412 seam triangles and zero textures. Their pre-repack buffer is 78,272 bytes / 39,353 gzipped. Runtime byte preservation also passes. These resource figures do not claim a measured frame rate.

The relationship graph assigns superior tibiofibular articulation, fibular collateral ligament and biceps femoris attachments to the head, talar articulation and ankle ligaments to the malleolus, and fibular artery supply and named muscle origins to the shaft. The common fibular nerve is adjacent to the neck, not claimed as bone innervation. Whole-assembly arterial queries retain their completeness gate. Links are reciprocal, ipsilateral and restricted to the native source frame; ambiguously mislabeled calcaneofibular source meshes are excluded. Actual attachment footprints and packing remain unverified.

Actual installed left and right whole renderings show continuous textured exteriors and four standard leader labels per side. The legacy whole link, Show components and head selection were reviewed; both component expansion and head selection preserve the exact serialized camera pose. The head inspector exposes its whole parent and seven related structures with proximal supply notes. Eight targeted fibular tests now include the actual runtime catalogue without a reconstructed fixture. General component checks, requested identity across 10,231 nonempty male/female concepts, and the complete 7,474-term hierarchy checks pass. The hierarchy generator now gives a unique reviewed shared-source identity precedence over an old name concordance, correcting the existing female proximal subclavian selector mismatch without changing geometry.

Viewer evidence: `outputs/ta98-female-fibula-installed-components.png` and `outputs/ta98-female-fibula-installed-right.png`. Installation/export evidence: `reports/ta98-female/fibula-{local-installation,export-validation}.json`. Twelve fine-feature terms in the inventory remain explicitly missing. Planar region boundaries, internal tissue volumes and whole adjoining anatomical packing remain unaccepted. This is local provisional construction progress, not completed fibular anatomy or completed whole-body TA98 coverage.
