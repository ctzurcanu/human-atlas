# Female TA98 development model

The local viewer now has **Female · TA98 review**:

<http://localhost:3016/?model=local-female-ta98>

This is a real independently selectable female mesh catalogue with a lossless on-demand runtime. It is an initial source assembly, **not completed TA98 anatomy**. The unchanged primary female skeleton establishes its coordinate frame. Female source anatomy is preserved; no male genitalia or male urethra is manufactured to meet a term count. The current primary body silhouette is moderately fit and was inspected in the dark viewer. No phenotype deformation has been made.

## Geometry and coverage

| Measure | Current result |
|---|---:|
| Full authority inventory | 7,474 terms |
| Applicable female/shared material terms | 5,903 |
| Male-specific terms excluded | 188 |
| Source-exact material correspondences | 749 |
| Material terms with at least one image-textured source part | 516 |
| Material terms whose every currently exact-mapped part has UVs and an image texture | 467 |
| Material terms without an exact source representation | 5,154 |
| Newly constructed distinct terms | 0 |
| Anatomically accepted complete terms | 0 |
| Primary female meshes | 1,950 |
| Hidden HRA donor meshes | 822 |
| Total meshes / named concepts | 2,772 / 3,591 |
| Primary / total triangles | 1,416,440 / 5,175,582 |

The 749 correspondences describe source availability. They do not certify every region/layer, topology, interfaces, packing, texture coverage or anatomical correctness. The strict texture count requires every exact-mapped part of the term to have both UVs and an image map. It does not prove that the texture is anatomically correct or that hidden HRA donors are registered. Copied source meshes and selection aliases do not count as newly constructed anatomy. Missing terms are retained in the full TA98 hierarchy.

## Source frame and placement order

The required phenotype is a **fit, idealized adult female**, with balanced proportions and moderate muscular definition, using the actual reference body female exterior and coordinated skeleton as the reference. This is a firm acceptance requirement, reaffirmed by the user. Compare front, side and back silhouettes against that source before accepting a shape change. Source proportions currently remain unchanged; this requirement does not certify the eventual completed model. Any proportion change must transform affected attachments and internal tissues together and pass the same packing checks.

The primary `.local-models/female.json` has 258 skeletal meshes and only two explicitly named fascial meshes: left/right plantar fascia. Its trunk fascia is absent. The first placement report records every native skeletal/fascial guide and the unchanged bounds; bounds are broad-phase guides, never certified tissue volumes.

Use native female bones first, then sourced/validated fascia and adjacent structures. For the female pelvis, preserve the actual pelvic bones, then establish bladder, rectum, vagina, ovaries, uterine tubes and vascular attachments before fitting the uterine assembly. A donor's bounding-box similarity does not establish registration. Never rescale a uterus independently to make it fit, treat another organ's lumen as free space, or create a false closed fascial volume from absent sheets.

HRA source donors for missing terms are hidden by default and retain explicit `female-source-registration-and-packing-pending` status. Their source frame is not presumed to agree with the primary female source. A source aggregate may include an alternative representation of a primary organ; these alternatives remain hidden. No merged donor assembly is accepted.

## Concrete uterine construction and review

The full 32-term authoritative uterine inventory is in `reports/ta98-female/uterine-assembly-inventory.json`, including representation classes, parentage, exact source IDs, missing tissue layers, glands, openings and interface evidence. The native pelvic cervix had inherited a false dental-neck correspondence (`A05.1.03.030`); it now resolves to uterine cervix (`A09.1.03.010`) in the review catalogue and source-aware hierarchy. This correction is not newly constructed anatomy.

The real native uterus/cervix, 15 pelvic skeletal guides and neighboring bladder, rectum, vagina, ovaries, tubes and vessels were staged (37 objects), image-textured and rendered together. Actual sagittal, axial and coronal triangle cuts were inspected. The native uterus is a closed genus-two wall with three openings; the cervix is genus one. Those genera do not establish histology or correct packing. Native left uterine-tube geometry has 9 nonmanifold edges, 3 duplicate faces and 6 degenerate faces. Source crossings include rectum/cervix 29 face pairs and cervix/vagina 2.

A direct Boolean union failed with 15 nonmanifold edges and remains rejected. A concrete 4 mm measured annular junction reconstruction preserves native geometry outside the collar, clips exact source edges, interpolates original UVs and joins the actual inner/outer wall contours without caps. The resulting whole wall has 1,284 triangles, one genus-two component and zero open/nonmanifold/winding/duplicate/degenerate faces or self-crossings. Front/back/side renders and 12 actual pelvic sections were inspected. It remains **uninstalled and anatomically unaccepted**: rectum/vagina crossings remain and the collar creates five additional vaginal artery crossing pairs requiring repair/interface classification.

Two independently selectable tissue regions, fundus (`A09.1.03.002`) and body (`A09.1.03.003`), were concretely cut from that repaired native source wall. This is a staged dissection, not source copying. The regional plane is explicitly provisional, derived from actual bilateral tube-interface coordinates; its 53 micrometre numerical displacement avoids float32 unstable sliver cells and does not assert an established anatomical boundary. All 1,284 whole faces are retained once or divided into nonoverlapping fragments, with one common exact vertex/normal array and original/interpolated source UVs. The 1,572-face regional assembly passes double and float32 topology, zero self-crossings, source surface/UV interpolation and independent positive-area overlap checks. Front/back/side and pelvic renders were inspected. **Two provisional dissected terms; zero newly installed or accepted terms.** No wall layer, cavity, opening or gland is falsely counted as a newly manufactured solid.

The individual reference body female exterior envelope and pelvic/femoral landmarks are recorded in `reports/ta98-female/reference-body-body-reference.json`, including native manifest and catalogue hashes. Its proportions remain unchanged, following the requested fit/idealized female source reference. The importer's 1.75 m coordinate height is a viewer convention, not subject anthropometry. Missing pelvic/trunk fascia remains explicit.

The left uterine tube now has a concrete source repair staged: six exactly collapsed zero-area faces and both faces of three proven opposite coincident pairs are removed. All retained source coordinates, UVs, native Int16 normals and image material are unchanged; signed wall volume is unchanged. The result has 952 triangles, closed genus one, zero topology defects at double/float32, and zero self-crossings. Its 14 exact uterine-interface coordinates remain unchanged. Actual whole/detail textured renders and six paired source/repaired triangle sections were inspected; the defective section's degree-four branches are eliminated. This repairs source integrity, not all fimbrial material representations. Packing still exposes inherited ovary/tube crossing (42 pairs) and native arterial/venous contacts (4/8), so it stays uninstalled.

A 36-object coupled whole pelvic trial includes the repaired uterine/cervical wall and left tube with unchanged native skeletal references. Its textured pelvic render and 12 actual sections are saved. The cervix/vagina source shares 25 exactly opposite coincident wall-face pairs and 38 coordinates. Cancelling both copies preserves UVs and removes duplicated tissue but leaves 13 nonmanifold contact edges in a separate ring. A 20 micrometre local contact-ring clearance trial preserves the actual joining seam and passes manifold/genus/float32 checks, **but introduces 25 self-crossing pairs and is rejected**. No source gate is relaxed and neither cervicovaginal trial is installed.

Next concrete construction: reconstruct the measured coupled rectum/cervix/vagina interfaces without closing real lumina; classify native vascular tissue relationships and eliminate the newly introduced collar artery contact. Only then can individual regional boundaries and missing material layers/ligaments be accepted. HRA donors must be registered against measured female landmarks and texture evidence before use.

Local artifacts:

- `outputs/ta98-female-uterine-interface-candidate.blend`
- `outputs/ta98-female-uterine-regions-candidate.blend`
- `outputs/ta98-female-uterine-regions-{front,back,side,pelvis-oblique}.png`
- `outputs/ta98-female-uterus-interface-sections.png`
- `outputs/ta98-female-left-tube-repair.blend` and `ta98-female-left-tube-{whole,defect}-{source,repaired}.png`
- `outputs/ta98-female-left-tube-sections.png`
- `outputs/ta98-female-pelvic-whole-candidate.blend` and `ta98-female-uterus-pelvic-whole-sections.png`
- `reports/ta98-female/left-tube-repair-{plan,validation,packing,sections}.json`
- `reports/ta98-female/cervicovaginal-interface-{plan,validation}.json`
- `reports/ta98-female/fornix-contact-{plan,validation,packing}.json`
- `reports/ta98-female/uterine-interface-{validation,packing,sections}.json`
- `reports/ta98-female/uterine-regions-{packing,independent-validation}.json`

## Runtime evidence

The runtime keeps every unsuppressed source mesh and its positions, normals, indices, UVs/colors and materials bit-for-bit. Hidden donor chunks are separate. Current total geometry is 216,407,376 raw bytes / 92,625,153 gzip bytes in 243 chunks. The full compressed runtime catalogue is 379,319 bytes, plus geometry and textures. Whole-body downloads remain substantial; these totals do not establish efficiency acceptance.

Inspected dark whole-body browser state: 79 required/loaded chunks, zero pending, 79 render calls, 1,397,132 submitted triangles. Hidden donor geometry was not loaded. The warm isolated primary uterus used one required chunk, one render call and 4,672 submitted triangles; prior loaded buffers remained cached. FPS, cold-load timing and GPU memory remain unmeasured. Browser review evidence is in `reports/ta98-female/browser-review.json`.

```sh
npm run build:ta98-female-review
npm run test:ta98-female-review
python3 scripts/validate-ta98-runtime.py --model female
node scripts/validate-ta98-hierarchy.mjs
node --experimental-strip-types scripts/validate-on-demand-chunks.mjs
npm run check
```

The runtime is the default Female TA98 local view; `&delivery=0` opens the unchanged review catalogue. The local development server rejects stale female runtime catalogues and serves the current catalogue compressed. The original Male TA98 mapping is unchanged after hierarchy regeneration, and its runtime verifier still passes.

## Licenses and acceptance

The local female source has unverified redistribution rights and stays local. The HRA catalogue preserves its recorded CC BY 4.0 provenance and source hashes; source identity and license records are preserved, not independently adjudicated. No local review model was deployed, committed or published. Anatomical acceptance and the full male/female coverage objective remain incomplete.

Evidence files: `terminology-and-source-inventory.json`, `coverage.json` (per-term strict texture and registration status), `placement-reference-inventory.json`, `source-preservation-check.json`, `runtime-repack.json`, and `browser-review.json`, all under `reports/ta98-female/`.

## Rebuilt female pelvic context (2026-10-02)

The rebuilt source is now staged as 121 unchanged local meshes with native pelvic bone and fascial guides, nearby bladder/rectal/vaginal tissues, tubes, ovaries, vessels and nerves. See `reports/ta98-ca2024/female-pelvic-context-review.md`. The complete 32-term uterus comparison and interfaces precede staging. Exact scene/resource preservation passed; both textured views and all nine source sections were inspected. Contact diagnostics remain provisional and this source frame is not registered to the required reference body frame. Neither this staging nor the earlier whole-wall candidate adds accepted anatomy.

## Rebuilt vaginal whole-wall trial

The exact sagittal cap cancellation trial preserves all vertex/image attributes and removes 3,036 internal partition triangles. It fails whole-wall acceptance because complete mucosal/muscular assemblies contain disconnected patches and remain open/nonmanifold. The 63-term full uterine/tubal/vaginal inventory, inspected renders/sections, unchanged cervical/tubal crossing evidence and exact component source-face map are linked from `reports/ta98-ca2024/female-pelvic-midline-review.md`. Keep the candidate local and uninstalled; no term is newly accepted.

## Exact sagittal cap-remnant continuation

Six detached patches are proven opposite, same-wall/layer zero-thickness cap remnants despite different left/right triangulations. Exact rational face-by-face coverage and independent index/resource preservation support removing their 340 triangles in a separate local candidate. Complete mucosal/muscular component counts improve to one/two respectively. The nonplanar 380-triangle patch remains, along with open/nonmanifold walls and interface defects. All four surface views and nine section panels were inspected. See `reports/ta98-ca2024/female-pelvic-cap-remnant-review.md`; no new anatomy is accepted or installed.


## Native cervix hierarchy correction (2026-10-03)

The runtime mesh already has system `reproductive` and the review catalogue correctly identifies uterine cervix. The shared viewer terminology still matched the ambiguous source name "Cervix" to dental neck (`A05.1.03.030`), moving it into the digestive browser. A reviewed reproductive-only mapping now resolves `LOCAL:female:cervix` to uterine cervix (`A09.1.03.010`) and uterus ancestry. The reproductive subbranch uses authoritative female/male TA98 prefixes before fallback names, placing this leaf under Female genital system. Geometry and local catalogue bytes are unchanged.

The user's exact view fingerprint `252.1vw60qe.1w76eks` hides cervix index 1451 (`14b`) in range `z5-14b` while leaving uterus/vagina visible. Changing that range to `z5-14a` restores only the cervix. Existing saved URLs retain their explicit hidden sets; the hierarchy correction does not silently change them. `uterovaginal-view-visibility.json` binds this diagnosis; `test:hierarchies` now checks native female/runtime cervix identity, reproductive membership, female path and absence of dental ancestry. No claim of anatomical junction repair follows from restoring visibility.

## registered collection FEMALE circulation (female-circ-v2, replaces the male retarget)

The review route `?model=local-female-ta98&circulation=registered-collection-review` now loads `female-circ/female-circ-v2`.

**Source.** The registered collection `human_20_female` cardiovascular system: 1,762 parts with its own registered-collection female heart. It comes from `fetch/viewer/data/female/catalog.json`, captured from `fetch/content/registered collection/female-complete`. It includes the complete female pelvic vasculature: uterine, ovarian, vaginal, clitoral and labial vessels. All 444 native reference body vessel and heart parts are hidden, so there are no doubles and no bridges. The male-retarget build (female-circ-v1, bridges, polish) is superseded.

**Pipeline** (`.local-models/geometry-workbench/bdf-circ/scripts/bdf_chain.sh`):
1. `bdf_extract`: decode BDA meshes with the hierarchy.
2. `bdf_anchors`: global similarity (cm → m, scale 0.00888), then 163 similarity anchors (median residual 2.2 mm). The anchors pair registered-collection female and reference body female bones by name, plus organs: heart, lungs, trachea, liver, kidneys, bladder, uterus, vagina, thyroid, oesophagus, pancreas, gallbladder, ascending and descending colon, and the skull composite.
3. `bdf_transfer`: blended anchor similarities. The registered-collection heart group moves rigidly with the heart anchor. The displacement is averaged over a 5 mm source-space ball, which removes the long-tube twisting that folded the radicular vessels.
4. `bdf_join`: global continuity solve over every source contact under 2 mm, plus branch-origin gaps up to 10 mm. Heart-attached joints are weighted 20×. Bone and skin are soft constraints in the same solve.
5. `bdf_polish`: ARAP fold polish that pins every source-contact vertex.
6. `bdf_install`: install the result.
7. `bdf_check` / `bdf_gap`: checks.

**Metrics:**
- Source joints: 36,330. 4 have a gap over 2 mm (max 2.7 mm); none are over 5 mm.
- Every great vessel and pulmonary vein meets the registered-collection female heart within 0.1 mm.
- Folds: 0.36%.
- Bone: 22.4k vertices inside, median 0.3 mm deep.
- Skin: 528 vertices outside.
- Brain: 5% of vertices are more than 10 mm from a vessel.
- 703 of 13.5k registered-collection heart vertices lie inside the reference body lungs.

**Organ-fit pass.** New anchors: left and right cerebrum (registered-collection gyri, lobes and white matter), cerebellum, pons, medulla, spinal cord, jejunum, ileum, caecum, appendix, spleen and stomach. Organ fits may rotate up to 60°. The stomach and the transverse and sigmoid colon fit at 32–46° rotation and twisted the mesenteric vessels, so they are left out of the transfer. Transfer displacement is averaged over 8 mm.

The continuity solve now also has organ constraints, applied at deformation-node level so thin tubes are not flattened:
- Hepatic, renal and segmental pulmonary vessels stay inside the liver, kidneys and lungs.
- Cerebral, splenic, vasa recta and uterine/vaginal branches stay within 2–3 mm of the brain, spleen, small bowel and uterus.

Vessel vertices inside or within 5 mm of their organ, before → after:

| Organ | Before | After |
|---|---|---|
| Liver | 89% | 98% |
| Kidneys | 93% | 100% |
| Lungs | 97% | 100% |
| Uterus | 90% | 95% |
| Cerebral | 86% | 96% |
| Spleen | 85% | 98% |
| Small bowel | 33% | 80% |

Other metrics:
- Brain: 1% of cerebrum/cerebellum vertices are more than 10 mm from a vessel (was 5%).
- Continuity: 4 joints over 2 mm, none over 5 mm.
- Folds: 0.77% (was 0.36%), the cost of placing mesenteric and omental vessels on the reference body bowel.
- Bone: 23k vertices inside, median 0.3 mm deep.
- Skin: 668 vertices outside.

**Smooth-field pass.** Folds came from the transfer itself. Neighbouring anchors (for example the uterus and the sacrum) send nearby registered-collection points along crossing paths, which folds space. Averaging over nearest neighbours could not fix this, because in dense regions it only spans a few millimetres. The transfer displacement is now averaged on a fixed 20 mm grid (normalized Gaussian, trilinear sampling). Local fit is then recovered by the solve: joints, bone, skin, and organ constraints at node level. The bowel constraint weight is 15.

Metrics (previous → now):

| Metric | Previous | Now |
|---|---|---|
| Folds | 0.77% | 0.27% |
| Source joints over 2 mm | 4 | 0 |
| Vertices outside skin | 668 | 422 |
| Small-bowel vessels within 5 mm | 80% | 86% |

- Liver, kidney, lung, uterus, cerebral and spleen vessels are at 96–100% within 5 mm.
- Brain: 1% of vertices are more than 10 mm from a vessel.
- Bone: 25k vertices inside, median 0.3 mm deep.
- Heart–lung overlap: 703 registered-collection heart vertices sit inside the reference body lungs (median 3 mm deep). The native reference body heart overlaps them more (median 15 mm), so the oversized reference body lungs cause the overlap, not the registered-collection heart fit.

**Whole-organ audit pass.**
- `bdf_organfit` now audits 17 organ groups.
- All anchors are back in the transfer, now that the 20 mm smooth field absorbs their disagreement (38 transfer folds). An adrenal anchor was added.
- Anchor subsampling is seeded per anchor, so adding an anchor no longer changes the others' fits. Pons and medulla fits rotated 28–39° and are rejected.
- New node-level organ constraints: colon, rectum, thyroid, adrenals, bladder, pancreas, gallbladder and oesophagus. Ovarian vessels get a terminal-only constraint: nodes within 4 cm of the ovary.
- A stomach constraint was tried and dropped: the small reference body stomach folded the gastroomental arcade.
- Pattern pitfalls found and fixed: "callosomarginal" matched "marginal artery", "arcuate" caught the foot, and "interlobar" caught the pulmonary trunk branch.

Share within 5 mm, before → after:

| Organ | Before | After |
|---|---|---|
| Stomach | 17% (median 40 mm off) | 25% (9 mm off) |
| Colon | 25% | 37% |
| Rectum | 56% | 75% |
| Thyroid | 50% | 64% |
| Adrenals | 45% | 77% |
| Bladder | 32% | 73% |
| Ovaries | 0% | 53% |
| Pancreas | 45% | 72% |
| Gallbladder | 48% | 89% |
| Oesophagus | 86% | 95% |

Liver, kidneys, lungs, uterus, brain, spleen and small bowel stay at 88–100%.

Other metrics:
- Continuity: 1 joint over 2 mm (2.3 mm), none over 5 mm.
- Folds: 0.49%.
- Skin: 522 vertices outside.
- Brain: 2% of vertices are more than 10 mm from a vessel.

**Identity and search pass.**
- **TA98 matching.** `bdf_ta98.py` matches unlabelled registered-collection female vessel names against the TA98 authority database (A12 only). It accepts exact name or synonym matches, side-aware names, and "X of Y" names where TA98 term X has Y as its exact parent or ancestor. It also maps head aliases (medial/lateral dorsal cutaneous branch → medial/lateral cutaneous branch; numbered anterior intercostal artery → anterior intercostal branches) and hepatic artery segment branches (→ anterior, posterior or lateral segmental artery).
- **Strict rules.** Containment matching was rejected after it mapped lumbar and anterior intercostal branches to posterior intercostal terms. No parent-term fallback is used.
- **Result.** 75 more vessels have TA98 terms, including every female pelvic vessel: uterine, ovarian, vaginal, clitoral, labial and the bulb of the vestibule. In total, 1,392 of 1,762 parts (79%) have a TA98 term, recorded as `ta98Match`. The remaining parts are registered-collection-specific venous branches and perforators that TA98 does not name.
- **Search.** The installer now adds one concept per registered-collection part with its TA98 term. Before this, viewer search found none of the registered-collection vessels, because the native vessel concepts pointed only at hidden parts. Search, selection, the TA98/THA/FMA/Latin identifiers and the TA98 hierarchy were checked in the viewer on the left uterine artery.
- **Names.** registered-collection name typos are corrected at install: Overian, atruim, hepartic, longitudal, postarior, vesic.

**Axial trunk pass.**
- **The defect.** A side-by-side render of the registered-collection source against the fit showed the abdominal aorta and IVC pinched and crooked toward the bifurcation, with the iliac veins bent into hooks. Soft-organ anchors (bowel, uterus, bladder, kidneys) and the 20 mm smoothing grid were dragging retroperitoneal trunks.
- **The fix.** The axial trunks (thoracic and abdominal aorta, IVC, common/external/internal iliacs and their divisions, azygos/hemiazygos, median sacral) are now transferred with bone anchors only. They are excluded from the general smoothing grid and get their own 40 mm trunk-only smoothing (`SMT`). The front view now matches the registered-collection source's smooth Y bifurcation.
- **Not a defect.** In the side view the aorta follows the reference body lumbar spine, whose L5 sits 25 mm further back than BD's (deeper lordosis). The iliacs therefore run forward to the groin, where the femoral heads match within 6 mm. This reflects the reference body skeleton.
- **Checks.**
  - Head, neck and limbs keep the registered-collection vessel shapes in the comparison renders.
  - Continuity: 3 joints over 2 mm, none over 5 mm.
  - Every great vessel meets the registered-collection heart.
  - Folds: 0.52%. Bone: 25k vertices inside, median 0.3 mm deep. Skin: 522 vertices outside.
  - Organ fit is unchanged.

**Full TA98 coverage.** All 1,762 registered-collection female circulation parts are now linked to TA98.

- **1,398 exact terms.** New in this pass:
  - aliases (perforator → perforating; gastro-omental, caecal and oesophageal spelling);
  - the unsided hepatic segment II branch;
  - curated synonyms: anterior labial arteries → anterior labial branches, radicular arteries → radicular branches.
- **364 parent placements** (`ta98Kind: 'parent'` on the concept; the part has no `ta98Term`), so the viewer files each vessel under its nearest TA98 vessel without claiming exact identity. Sources:
  - the "X of Y" name chain (posterior intercostal veins, lumbar veins, circumflex femoral veins, subcostal…);
  - a curated anatomical hint table:
    - cecal veins → ileocolic vein;
    - portal segment branches → right/left branch of the portal vein;
    - inferior polar artery → renal artery;
    - alveolar and infraorbital veins → pterygoid plexus;
    - hypophyseal veins → portal veins of hypophysis;
    - vasa recta → jejunal/ileal arteries;
    - foot and leg veins → anterior/posterior tibial veins or great saphenous vein;
  - registered-collection group names (pulmonary veins, veins of spinal cord).
- **Exclusions.** Male-only TA98 terms are excluded. No vessel falls back to a generic "arteries", "veins" or "cardiovascular system". Methods are recorded per concept (`ta98Method`).
- **Viewer check.** The left dorsalis pedis vein shows under IVC › … › popliteal vein › anterior tibial veins.

**Breast pass.**
- **Anchors.** registered-collection breast tissue (fatty tissue and axillary tail) is paired with the reference body "Breast fat (left/right)" on each side: scale 1.08, rotation 13°.
- **Surface constraint.** A new `surf` node-level constraint keeps the lateral and medial mammary branches (arteries and veins) on the breast surface, within 4 mm. This reproduces the registered-collection source relation: 0% inside registered-collection fatty tissue, 93–99% within 5 mm. A looser "near" rule was tried and rejected: it let the vessels sink inside and moved their convergence 21 mm off the breast apex.
- **Not usable.** The hidden HRA nipple sits 24 mm off the reference body skin and 9 cm too low, so it isn't used. BD's own convergence point is its breast apex (3 mm off); the registered-collection nipple is 22 mm lateral of it.
- **Result.**
  - 95–100% of mammary vessels are within 5 mm of the reference body breast; none are outside the skin.
  - The anterior convergence is 5–7 mm from the reference body breast apex, 5 mm under the skin.
  - Continuity is intact: none over 5 mm except one artery–vein neighbour pair (external iliac, 5.3 mm).
  - Folds: 0.52%.
- **Audit.** `bdf_breast.py` checks this.

**Nipple pass.**
- **Where the reference body nipple is.** The nipple is on the reference body skin at (±0.10, 1.288, 0.12). Two independent signals agree: a dense 3.2 mm skin bump, and the darkest areola pigment in the skin texture (`female-texture-11`, sampled through the skin UVs).
- **What was wrong.** That is 2 cm lateral and 2 cm above the breast-fat apex. The HRA breast set (nipple, areola, areolar tubercles, lactiferous sinuses and ducts, mammary lobes, suspensory ligaments, interlobar fat) sat about 10 cm too low and 4 cm too deep, which is why it was hidden.
- **The fit** (`bdf_breastparts.py`), per side:
  - a similarity ICP of the HRA interlobar fat onto the reference body breast fat, with the HRA nipple tip weighted onto the skin areola centre (scale 0.79–0.83, rotation 9–11°);
  - the areola is projected 0.3 mm onto the skin, and the nipple, tubercles and sinuses follow it (nipple tip within 3 mm of the areola centre);
  - lobes, ducts, sinuses, ligaments and fat are pushed at least 1.5 mm under the skin, with smoothed displacement (no vertex outside).
- **Visibility.** Nipple, areola, tubercles, lobes, ducts and sinuses are now visible. The HRA interlobar fat duplicates the reference body breast fat and stays hidden, as do the ligaments.
- **Bug avoided.** The suspensory ligaments of the ovary and of the lens matched a loose "suspensory" pattern; the pattern now names only breast parts.
- **Checks.** In the viewer the nipples sit on the breast mounds with nothing breaking the skin. The mammary vessels keep their registered-collection relation to the nipple (nearest 8 mm vs registered-collection 12 mm; median 52 vs 48 mm).

## registered collection FEMALE nervous system (same catalogue: female-circ-v2)

**Route.** `?model=local-female-ta98&nerves=registered-collection-review` (alias of the circulation route). The catalogue holds the circulation and the nerves.

**Source and replacement.** 1,487 parts of the registered collection `human_20_female` nervous system: 155 central (brain, brainstem, spinal cord) and 1,332 peripheral (cranial and spinal nerves, plexuses, ganglia, the eyes' inner structures). They replace all 237 native reference body female nerve parts, which are now hidden.

**Pipeline** (all in `bdf_chain.sh`):
1. `bdn_extract`.
2. `bdf_transfer`, with the same 181 bone/organ anchors and the 20 mm smooth field (`SYSF=nervous_system`).
3. `bdf_join` continuity solve. Organ groups are off. Ocular and intra-osseous structures (labyrinth, cochlea, alveolar and mental nerves) are exempt from the skin/bone push.
4. A new bone-proximity constraint (`BPROX`): nodes within 10 mm of bone in the registered-collection source are pulled back to their source distance when they drift more than 3 mm.
5. `bdf_polish`.

**TA98.** 1,177 terms come from the male registered-collection reviewed mapping, matched by objectId (names differ between the sexes). 230 come from the authority matcher (A14/A15/A13), with nerve rules: spinal roots, rootlets and rami; region-specific posterior-ramus branches; clunial, subscapular, mammary and phrenic nerves; and CN III/IV nuclei. 78 are parent placements. Every nerve is labelled.

**Metrics:**
- Continuity: 26 joints over 2 mm (all source near-contacts), none over 5 mm.
- Folds: 0.13%.
- Brain cortex vs reference body cerebrum: 2.3–3.3 mm median (both directions). Cerebellum: 3.5 mm. Spinal cord: 1.2 mm.
- Nerve-to-bone distances keep the registered-collection relation: sciatic 19 → 25 mm, femoral 25 → 22, median 10 → 8, ulnar 8 → 9. The radial nerve drifted 10 → 21 mm; the bone-proximity constraint brought it back to 16.
- Viewer: search, selection and TA98 checked on the sciatic nerve (A14.2.07.046, under the sacral plexus).

**Nipple-centred breast vessels.**
- **The problem.** The registered-collection breast anchor placed BD's nipple (anterior tip of its lactiferous ducts) 2.6 cm medial and above the true reference body nipple, so the mammary vessels radiated from the wrong centre.
- **The fix.** The breast anchors now carry the nipple as a weighted landmark (60% of the matched points) onto the reference body skin nipple (areola pigment + skin bump). Mapped nipple error: 0 mm. Shape residual: 7.7/5.0 mm.
- **Vessel distance from the nipple** (fit vs BD):

  | | Fit | registered-collection |
  |---|---|---|
  | Minimum | 11 mm | 12 mm |
  | q5 | 15–16 mm | 16 mm |
  | q25 | 31–32 mm | 29–30 mm |
  | Median | 49–51 mm | 48–49 mm |

- **Other breast checks.** 95–100% of mammary vessels are within 5 mm of the breast, none outside the skin. In the viewer they radiate to the nipple-areola hub.
- **Full rebuild.** Circulation and nerves rebuilt with the new anchors. Continuity: circulation 4 joints over 2 mm, nerves none over 5 mm. Folds: 0.50% circulation, 0.13% nerves. Organ fit unchanged. Audit: `bdf_nipple2.py`.

**Vessel refinement.**
- **Solve settings.** The circulation solve now runs 12 rounds with bone/skin weight 12, plus the bone-proximity constraint (46.7k nodes that sat within 10 mm of bone in BD; 10.1k → 1.8k drifted).
- **Installed result.**
  - Bone: 28.3k → 26.1k vertices inside (median 0.4 mm deep).
  - Skin: 514 → 477 vertices outside.
  - Continuity unchanged: 4 joints over 2 mm; one artery–vein neighbour pair at 5.5 mm.
  - Folds: 0.55%.
  - Nipple relation unchanged.

## registered collection FEMALE lymphatic system (same catalogue: female-circ-v2)

**Route.** `?model=local-female-ta98&lymphatics=registered-collection-review`.

**Source and replacement.** 306 parts of the registered collection `human_20_female` lymphatic system: vessels, node groups, trunks, the thoracic and right lymphatic ducts, thymus and spleen. They replace the 57 native lymphatic parts, which are now hidden.

**Pipeline.** The same anchors and smooth field, then the continuity solve with bone proximity. One new cross-system constraint (`XREF`): every lymph vertex within 5 mm of a blood vessel in the registered-collection source is tied to that vessel's fitted position, keeping its registered-collection offset. This applies to 189k vertices (77%).
- Lymph–vessel source contacts over 5 mm apart: 69 → 4 (over 10 mm: 31 → 0).
- The remaining offsets are the paramammary nodes (7–9 mm from the lateral thoracic veins, which moved with the breast fit) and one external iliac node (6.5 mm).

**Metrics.**
- Within-lymph joints: all at most 2.5 mm, except 3 between 5 and 10 mm after the vessel ties.
- Folds: 0.76%.
- Bone: 1.7k vertices inside. Skin: 111 outside.

**TA98.**
- 195 terms carried over from the male registered-collection lymphatic mapping (by objectId).
- 5 new exact terms: para-uterine and paravaginal nodes, prelaryngeal nodes.
- 106 parent placements. Lymphatic vessels are filed under their draining node group (gastric, external iliac, superficial inguinal, colic…) or under the regional groups (upper/lower limb, craniocervical, thoracic).
- Every part is labelled.

**Lymphatics kept inside the skin at the hips.**
- **Cause.** The superficial inguinal lymphatic vessels have long edges (up to 31 mm) that cut across the hip concavity. Their vertices tested inside, but the edges were up to 22 mm outside the trunk skin, visible between hip and hanging arm.
- **Why vertex checks missed it.**
  - The watertight skin proxy bridges the arm–hip gap.
  - Against the real skin, the hand's skin made those points read as "inside".
  - Pushing to the nearest skin point moved them into the hand.
- **Fix.**
  1. The lymphatic solve tests edge points against the real skin (`EDGESKIN`).
  2. A post-pass (`bd_skinfix.py`) works on the trunk skin only (arm/hand skin faces removed via nearest-bone labelling). It subdivides affected trunk lymphatic parts twice (1-to-4) and pulls outside points 1 mm under the skin with smoothed displacement. Elbow and limb structures use the full skin.
  3. The installer applies the result with `SKINFIX`.
- **Result.** Hip lymphatics are 0 mm outside (were 22 mm). The worst remaining is 1.9 mm, on the upper back near the shoulder opening of the trunk skin. The viewer shows no protrusion at the hips.

**Lymphatic fit and labels, refined.**

*Fit.*
- **Vessel ties per node.** Per-vertex vessel ties stretched nodes that lie between two vessels (up to 9× edge stretch, 0.78% folds). Ties are now per deformation node, using the average accompanying-vessel motion (`XNODE`, weight 8), plus light per-vertex ties (weight 3) only on true contact vertices (within 1 mm in BD).
- **Priority over bone proximity.** Vessel ties now override bone-proximity targets. Bone proximity had let the sacral nodes slide 37 mm along the sacrum away from the median sacral vessels.
- **Result.**
  - Folds: 0.57%.
  - Internal lymph joints over 5 mm: 2.
  - Lymph–vessel contacts over 5 mm: 6, none over 10 mm.
  - Hips kept inside the skin by the skin fix (22.7 → 2.0 mm).

*Labels* (`bdl_ta98_refine.py`).
- **Exact node terms refined to TA98 subgroups:**
  - lateral/medial/intermediate external and common iliac nodes;
  - superior/inferior pancreaticoduodenal nodes;
  - anterior and lateral groups of the superior diaphragmatic nodes → prepericardial / lateral pericardial;
  - superior anorectal → superior rectal nodes; retroaortic → postaortic nodes.
- **Vessels filed at their first TA98 node station:**
  - anterior/lateral cervical superficial and deep groups; facial nodes;
  - brachial, superficial and deep nodes of the upper limb;
  - anterior/posterior tibial nodes, popliteal deep/superficial nodes, inferior superficial inguinal nodes, deep inguinal nodes;
  - subscapular nodes (back); tracheobronchial nodes (cardiac trunks); brachiocephalic nodes (thymus); left/right lumbar nodes (ovary).
- **Cleaned up.** Coarse male "lymph node" exact terms were removed, and no part is left on a broad regional group. Final: 197 exact, 109 parent placements.
- **Viewer check.** Left lateral external iliac nodes → TA98 A13.3.04.012 under pelvic › parietal › external iliac nodes. Left deep popliteal lymphatic vessel → under popliteal nodes › deep nodes.

**Lymph nodes kept out of organs.**
- **The audit.** `bdl_organs.py` found 86 lymph parts more than 5% inside a reference body viscus. Examples:
  - paratracheal nodes 100% inside the trachea;
  - retrocaecal nodes inside the colon;
  - pararectal nodes inside the rectum;
  - parauterine nodes inside the uterus;
  - sigmoid nodes inside the small bowel.
- **The constraint.** The lymph solve now pushes nodes out of 26 visceral organs (`OUTORG`).
  - Node level: 3 mm margin, weight 30, applied after the vessel ties so it takes priority.
  - Vertex level: weight 3.
  - Exempt (intra-organ by anatomy): intrapulmonary and bronchopulmonary nodes, intraglandular parotid nodes, splenic and hepatic lymphatic vessels, spleen, thymus, tonsils.
- **Result.** 86 → 40 parts more than 5% inside, most of them the exempt ones.
  - Remaining: prepericardial nodes against the oversized reference body left lung; the right parauterine nodes (43%); a few mesenteric groups inside the reference body bowel loops.
- **Junction metric.** The raw lymph–vessel contact metric counts incidental crossings, e.g. external iliac lymphatics crossing the ovarian artery. The namesake check (`bdl_companion.py`) is the real one: 67 lymph groups keep their distance to their namesake vessel (median change 0.1 mm). The right common iliac nodes are now 6–7 mm farther (registered-collection 3 mm) where they are pushed out of the bowel.
- **Other metrics.** Folds 0.75%. Internal joints: 3 over 5 mm. Skin: the hips stay inside; the hand dorsal vessels remain 3.5 mm out, because the hand skin is too thin for the skin fix.

**Node-weighted organ exclusion and detailed TA descriptions.**
- **Fit.** The organ exclusion now pushes lymph nodes harder than lymphatic vessels (vertex weight 10 vs 3), so nodes leave the organs while vessels keep following their blood vessels. Right parauterine nodes inside the uterus: 43% → 26%. Prepericardial nodes inside the (oversized) left lung: 72% → 36%. Namesake-vessel companionship is unchanged (median 0.1 mm). Folds: 0.76%.
- **Descriptions.** `bd_describe.py` writes a TA-based description for every registered-collection part (lymph 306, nerves 1,485, vessels 1,758). The installer attaches it as `description`, which the viewer shows. Each description gives:
  - the TA98 ID, English and Latin terms, and FMA ID; or, where TA98 has no term, the TA98 group it is filed under;
  - the TA98 path;
  - the neighbouring structures it is continuous with in the registered-collection source, each with its TA98 ID;
  - for lymph and nerves, the blood vessels it accompanies.
- **Viewer check.** The right external iliac lymphatic vessels show their filing (A13.3.04.009, *nodi iliaci externi*), the six external/common iliac node subgroups they connect to (with TA98 IDs), and the external iliac artery and vein they run with.

**Hand lymphatics.**
- **The problem.** The dorsal lymphatic vessels of the hands cross the gaps between fingers. Pushing each point to its nearest skin pulled neighbours toward different fingers and stretched 1.9 mm edges to 9 mm.
- **The fix.** `bd_skinfix.py` now moves hand parts coherently toward the centre of that hand's metacarpals, and smooths the push field without pinning for limb parts.
- **Result.** The hands are at most 2.0 mm out (were 3.5 mm). Every lymphatic part is now at most 2.0 mm outside the skin.

**Lymphatic deep pass.**
- **Lymphoid organs** (`ORGFIT`). The registered-collection spleen, thymus and palatine tonsils are fitted by their own similarity ICP to the reference body organs, like the heart:
  - spleen: scale 1.32, two-sided 5 mm;
  - thymus: scale 1.65, about 3 mm;
  - palatine tonsils: about 0.5 mm.
  
  They are held fixed in the continuity solve, so their vessels and nodes attach to them. Before, the registered-collection spleen was 4.0 × 5.3 × 8.3 cm against the reference body 6.8 × 8.1 × 8.9 cm, and the thymus was 3 cm off.
- **Namesake ties** (`NAMETIE`). Lymph groups named after a blood vessel keep their registered-collection offset to its fitted position (node level, weight 8). This fixed the right common iliac drift (6–7 mm). Thymic lymphatics are excluded and follow the refitted thymus.
- **Fold polish** (`POLUNPIN`). The polish unpins folded faces plus two rings, so ARAP can repair dense abdominal networks. Folds: 0.78% → 0.57%.
- **Checks.**
  - Lymph joints: 3 over 5 mm, none over 10 mm.
  - Namesake companionship: median change 0.1 mm.
  - Skin: at most 2.0 mm outside everywhere.
  - Organ penetration: 39 parts more than 5% inside, mostly legitimate (intrapulmonary, hepatic, splenic).
  - Remaining overlaps: prepericardial nodes against the oversized reference body left lung (36%), right parauterine nodes and vessels (24–28%), right colic vessels in the reference body ileum (28%).

## Eyes (open) — registered-collection eyeball fitted to the reference body sockets, full TA98

- Scripts: `bdf-circ/scripts/bdf_eyes.py` (writes `eyes.pkl`, `eyes_report.json`), `bd_eye_ta.py`; installer flag `EYES=1` (run `bdf_eyes.py` on an atlas installed without EYES).
- Fit: per side, a sphere fit of the reference body eye exterior (cornea bulge excluded) and of the registered-collection sclera. Similarity: scale 0.955/0.961, axis turn 8.7°/10.0°, shift 5.6/6.1 mm. The registered-collection globe now sits in the reference body globe position, inside the open palpebral fissure of the reference body skin.
- What moves:
  - Globe parts move rigidly.
  - Orbital registered-collection structures of every system follow with a smooth falloff: optic nerve, ciliary nerves and ganglion, ophthalmic vessels, vorticose veins, lacrimal gland.
  - Structures on the globe surface (within 1 mm) move fully.
  - Eyelid and face structures in front of the globe stay with the lids (palpebral arches and veins).
- Open eyes:
  - The registered-collection palpebral conjunctiva is one sheet across the eye. Faces with no reference body skin in front of them are removed, tested straight on and at 15–20° oblique views (615/632 of 1180 faces), so the remainder lines the lids only.
  - The duplicate reference body eye shells ('Eye exterior source surface', 'Eye source cap'; parent-only TA) are hidden by default.
- Materials:
  - Opaque: white sclera, brown iris, black pupil, retina, choroid, ciliary parts.
  - Transparent: cornea, lens, vitreous, bulbar conjunctiva, aqueous humor.
- HRA retina and chamber parts that registered-collection lacks are fitted by the same globe similarity and shown:
  - macula lutea (now exact A15.2.04.021)
  - fovea centralis
  - optic disc
  - ora serrata
  - aqueous humor
  - scleral venous sinus
  - Macula, fovea and disc are snapped radially onto the registered-collection retina (median ≤0.5 mm).
- TA98 curation (`bd_eye_ta.py`):
  - optic nerve → A14.2.01.006 optic nerve [II]
  - supra-orbital lateral and medial branches → A14.2.01.022 and A14.2.01.023
  - anterior ethmoidal internal nasal branches, lateral and medial → A14.2.01.032 and A14.2.01.033
  - sensory root of the ciliary ganglion → A14.2.01.026
  - parasympathetic root → A14.2.01.010; it was wrongly the sublingual-ganglion root before
  - meningeal branch of the anterior ethmoidal nerve and middle superior alveolar artery → parent filing (TA98 has no own term for either)
- Eye and orbit coverage: every registered-collection eye part is exact TA98:
  - sclera, cornea, corneoscleral junction, iris, pupil, lens, zonular fibres
  - choroid, ciliary muscle, ciliary processes, retina, vitreous body, trabecular tissue
  - bulbar and palpebral conjunctiva, lacrimal caruncle, gland, canaliculi and sac
  - ciliary ganglion and its three roots, and the orbital nerves and vessels
- Eyelids, narrowed after review:
  - The skin opening is a UV cut-out (`skinEyes`) on a lid rim that curls into the orbit. Before: 28.5 × 12.6 mm, with the whole iris exposed top to bottom.
  - `bdf_eyes.py` moves the lid geometry so the cut-out rides on its UVs: the upper lid comes down 1.2 mm (`LIDUP`) and the lower lid up 0.4 mm (`LIDLO`) at the centre. The shift tapers to 0 at the canthi and 9 mm from the margin.
  - The lids rest 0.4 mm on the fitted eyeball. The eyelashes and the palpebral conjunctiva move with them.
  - The skin keeps its UVs; normals are welded across UV seams.
  - The conjunctiva visibility cut treats cut-out skin as see-through.
- Eyelid skin, final version (replaces the narrowing above). The web-model head had a 512 px face paint, and its eye opening was a UV cut-out on a rim that curls into the orbit, which left pale bands and a hard cut. Now:
  - `glb_head.py` reads the original reference body head skin from `Female_2024.glb`: 16 pieces, triangle strips decoded.
  - `face_bake.py` aligns it to the atlas head by similarity ICP (face residual median 0.8 mm). It bakes the 2048 px face texture into the atlas head UV layout (`head-open-eyes-v2.jpg`), colour-matched to the old head and body skin (gain 1.24/1.08/0.92) and feathered into the old texture away from the face.
  - `glb_eye.py` gives the original eyeballs (24.3 mm sclera, iris centre ±28.5/1616.2/79.8 mm). `bdf_eyes.py` now fits the registered-collection eyeball onto them, so the lids rest on the cornea.
  - The original high-resolution eyelid and periorbital patches (about 6.6k vertices per side) are added as 'Eyelids (left/right)' (TA98 A15.2.07.024), with the colour-matched 2048 px texture (`eyelids-v1.jpg`) and a polygon offset.
  - 443 atlas head faces under the patches are removed, and the surrounding atlas skin is snapped onto the original surface within 4 mm (seam).
  - The head material no longer uses `skinEyes`, so there is no cut-out and no repaint. Skin positions stay in world coordinates, because the neck tone blend uses vertex height.
  - Opening: original reference body, 21.6 × 8.0 mm (front coverage; the old opening was 28.5 × 12.6).
  - Each lash strand moves rigidly, root onto the nearest original lid-margin point. The palpebral conjunctiva follows the lids, then is cut to what is not occluded by skin or lids.
- Merged eyelids and lashes, final:
  - The original eyelid patches are merged into the head Skin part, so the Skin alone is complete. They use one 4096 × 2048 texture (`head-eyes-v3.jpg`: baked face on the left, original face on the right). 551 head faces under the patches are removed (including the rim curling into the orbit). Boundary faces that reach under a patch tuck 0.4 mm behind it, so there is no gap and no z-fighting.
  - Lashes: only the original reference body lash strands from the lid patches are kept (141/139 strands per eye, strip fragments welded). Strands that enter the fitted eyeball or root more than 1 mm from the lid surface are removed (19/12).
  - Both web-model lash sets (coarse sensory shards and the 30-strand integumentary set) are removed.

## Lymph — crease repair and TA98 vessel identities

- Creases: the old fold metric (normals against the registered collection source) also counts vessel segments the fit legitimately turned by more than 90°. True creases are adjacent faces turned back to back while flat in the source.
  - `bdl_untangle.py` repairs each crease cluster separately. The region becomes the crease-free transferred shape W plus a harmonic displacement from the region's boundary. It is accepted only if creases drop and nothing moves more than 3 mm.
  - Creased edges: 1,970 → 1,541. Fold metric: 0.57% → 0.53%.
  - Contacts and joints are unchanged in the audits, lymphoid organ fit is unchanged, and skin stays at most 2.0 mm outside.
  - The rest are in squeezed pelvic and mesenteric vessels, where the fit displacement itself is steep. Repairing those would move vessels more than 3 mm off their fitted course.
- TA98 (`bdl_ta98_vessels.py`, after `bdl_ta98_refine.py`):
  - Lymphatic vessels are filed under their TA98 vessel type instead of a node group: superficial lymph vessel A12.0.00.039 (26) and deep lymph vessel A12.0.00.040 (76). The TA98 node station each one drains to is kept, and the description says "Drains to: TA98 …".
  - Axillary lymphatic vessels → axillary lymphatic plexus A12.4.01.004 (exact).
  - Paramammary lymph nodes → paramammary nodes A13.3.02.002 (exact).
  - Lymph: 201 exact, 105 parent placements, of 306. The remaining parents are vessels (TA98 has no individual vessel terms) and two cervical node groups TA98 lacks (juxta-oesophageal cervical, sternocleidomastoid).

## Muscles — registered collection female muscular system

- Source: registered-collection `human_20_female_muscular_system`, 548 parts (230k vertices), via `bdm_extract.py` and the anchor transfer (`bdm_transfer.pkl`).
- Fit (`bdm_fit.py`, run in two stages):
  1. 310 registered-collection muscles get a reference body target:
     - same name;
     - registered-collection subdivision onto the reference body muscle containing it (one-way);
     - registered-collection whole muscle onto the union of its reference body parts (two-way).
     - Aliases cover reference body naming (internal pterygoideus = medial pterygoid, iliotibial band = iliotibial tract, and others), and a centroid gate rejects distant targets.
     - Fit: similarity ICP, then ARAP-regularised non-rigid ICP at the strongest strength that keeps creases under 0.5% of edges. Median surface distance 1.38 mm, p90 3.6 mm.
  2. The other 238 muscles follow the fitted displacement field through ARAP, so they keep their shape.
  3. Bone and skin: a batched ARAP solve with sticky targets. Inside bone 14.3% → 2.5% (median depth 0.13 mm, p90 0.7 mm; mostly the rotatores wedged between vertebrae). Outside skin → 0.01%.
  - Creases: 0.62% of faces.
- Native reference body and HRA muscular parts duplicated by registered-collection are hidden by default, by name or TA98 term, or by volume overlap (`bdm_dups.py`: at least 60% inside or within 4 mm). Ligaments, membranes, cartilage, retinacula and the trochlea stay.
- TA98: base seeded from the reviewed reference body labels (268), then `bdf_ta98.py` (A04/A05/A06/A09/A15), then curated rules (`bdm_ta98_refine.py`). Result: 473 exact and 75 parent placements of 548. The parents are mostly numbered interossei and lumbricals under their TA98 group terms.
- Distinct exact TA98 muscle terms in the female model: 252 (was about 191).
- Descriptions: TA98 identity and path, muscles in contact, nearby vessels and nerves.
- Muscles, refinement pass:
  - registered-collection muscles are subdivision-surface cages (median edge 4.6 mm, up to 27 mm; registered collection smooths them at render time). `bdm_subdiv.py` Loop-subdivides them to a median edge of 4 mm or less, capped at 6,000 triangles per muscle: 373k → 849k triangles.
  - Subdivision-shared reference body targets: registered-collection subdivisions sharing one reference body muscle (for example the medial pterygoid heads) are fitted together as a union, keeping their registered-collection arrangement. Stage 1 now runs in 3 shards.
  - Fit to the same-named reference body muscles: median 0.81 mm (p90 2.2 mm), was 1.38/3.6.
  - Inside bone: 1.63% of vertices (median depth 0.15 mm, p90 0.6 mm). Outside skin: 0.
  - `bdm_separate.py`: in each interpenetrating pair, the muscle more inside the other moves out (ARAP, sticky targets, bone still a constraint). Pairs that already overlap in the registered-collection source are kept. Interpenetration over 1 mm: 1.16% of vertices, below BD's own 1.76%.
  - Crease repair (`bdl_untangle.py` with `TRF=bdm_transfer.pkl`): 7,253 → 5,675 creased edges (0.66% of faces).
- Muscles, malformation fix (after review: "abdominal muscles malformed"):
  - Cause: 237 registered-collection muscles had no reference body target and only followed the neighbours' displacement field. reference body models some muscles as one mesh across both sides (rectus abdominis, external oblique, transversus abdominis, trapezius, diaphragm), names some with digits ("Trapezius1"), and large muscles failed a fixed 40 mm centroid gate.
  - Matching fixes:
    - unsided reference body meshes are cropped to the registered-collection muscle's side;
    - gates scale with muscle size;
    - digits and letter suffixes are stripped;
    - aliases added (transversus abdominis, inner intercostals, peroneus for fibularis, subclavian, interosseus);
    - iliotibial tract and tensor fasciae latae are fitted jointly.
    - Result: 395 of 548 targeted.
  - Fit: normal-aware point-to-plane non-rigid ICP (thin sheets keep both faces apart); local kinks relaxed. Median 0.76 mm, p90 1.8 mm.
  - Field-following muscles: crease-guarded ARAP, with a rigid fallback.
  - Bone-push creases repaired from the pre-push reference (`bdl_untangle.py` with `REFF=bdm_fit_stage2.pkl`), then `bdm_smooth.py` (local Taubin). `bdm_fixinv.py` replaces inside-out or collapsed muscles (22) by their similarity-aligned registered-collection shape.
  - Hybrid display (`bdm_usez.py`):
    - 319 registered-collection muscles that are one-to-one the same muscle as a native reference body muscle are shown by the native (textured, well-formed) mesh, which carries the registered-collection TA98 term where it had none.
    - registered-collection subdivisions and groups and the registered-collection-only muscles are shown as fitted registered-collection parts.
    - Group concepts keep every term: the registered-collection whole term over native heads (8), and native whole terms over registered-collection subdivisions (26).
    - 251 distinct exact TA98 muscle terms over visible muscle parts.
  - registered-collection muscle material tuned to the rendered tone of the reference body muscle textures.
- Muscle colour (living tissue):
  - Every textured muscle material (reference body textures, and registered-collection parts through a triplanar procedural fibre texture `muscle-fibre-v1.jpg` built at the average reference body muscle albedo) gets the gain `skinToneGain` [1.0, 0.74, 1.22] (`MUSCLEGAIN`), a pink-red living tone. Installer step under MUSCLES.
  - registered-collection parts and native parts now render through the same texture path, so they match from every view.
  - registered-collection heads and parts with exactly one reference body counterpart (for example the gastrocnemius heads) are shown by the native mesh: 333 registered-collection muscles in all.
- Rectus abdominis as a group (`bdm_rectus.py`, installer `extra_parts.pkl`):
  - The native reference body rectus is one mesh across both sides, with its anatomy painted on the texture. It is upsampled linearly, and each face is classified by its texture: bellies are labelled in texture space at full resolution, separated by the dilated pale bands.
  - It is split into left and right bellies 1 (superior) to 4 (inferior), left and right tendinous intersections, linea alba (A04.5.01.022) and left and right linea semilunaris (A04.5.01.025). Same material and UVs; normals welded across UV seams.
  - "Rectus abdominis" (A04.5.01.001) is now a group: its dropdown contains bellies 1–4 and the tendinous intersections, each pairing left and right. Group concepts: left and right rectus abdominis, and tendinous intersections (exact A04.5.01.002).
- Left/right split of the remaining single-mesh paired muscles (`bdm_split.py`, `extra_parts_split.pkl`):
  - Split, as shown native meshes (faces by midline side; same material and UVs): serratus posterior superior and inferior, coccygeus, iliococcygeus, transversus abdominis, external oblique, puborectalis, genioglossus, mylohyoid, depressor septi nasi, procerus and trapezius.
  - Each muscle is a group concept with its TA98 term (left and right in its dropdown). Each half is selectable, with its opposite side.
  - The diaphragm and orbicularis oris stay whole, as unpaired muscles. The auricularis superior and buccinator meshes cover one side only.
  - Installer: any number of split bundles. 254 distinct exact TA98 muscle terms over visible parts.
- Muscle labels are TA98 (installer, MUSCLES):
  - Every muscular part with a TA98 term is labelled with the TA98 English name plus side. Subdivision terms take their TA98 parent ("Long head of triceps brachii", "Upper head of lateral pterygoid"). Names TA98 uses twice take hand or foot from the TA98 code range ("Lumbricals of hand").
  - Parts that TA98 doesn't name show the TA98 term they are filed under, then the source qualifier ("Dorsal interossei of hand, first", "Longus colli, inferior oblique part"). The source name stays in `sourceName`.
  - The rectus split now has TA98 structures only: rectus abdominis (left/right), tendinous intersections (left/right, part of rectus abdominis), linea alba and linea semilunaris (left/right). No invented "bellies".
  - Curated: adductor hallucis oblique and transverse heads (were wrongly filed to adductor pollicis); lateral pterygoid upper and lower heads; cricothyroid straight and oblique parts; external anal sphincter subcutaneous, superficial and deep parts. Native parts: masseter superficial and deep parts, trochlea of superior oblique, patellar ligament, nasal cartilages, digastric fibrous loop (filed under digastric).
  - Every visible muscular part has a TA98 term. 277 distinct exact terms (264 muscle terms).

## Quick-win TA98 curation of native parts (`native_ta.py` → `native_ta.json`, installer)
- 652 native reference body/HRA parts (skeletal, connective, digestive, integumentary, respiratory) given a verified TA98 id and a TA98 label: 476 exact, 176 filed under their TA98 entity.
  - Teeth: typed by position (incisor, canine, premolar, molar, third molar), labelled with the TA98 type plus the FDI number ("Incisor tooth [FDI 11] (upper right)").
  - Intervertebral discs, with levels ("C V - C VI").
  - Phalanges of hand and foot by row; ribs (first, second, rib [n]); costal cartilages; manubrium; cranium.
  - Vertebrae C3–L5, metacarpals and metatarsals are filed under their TA98 groups (TA98 has no individual terms).
  - Ligaments by name match or curation, with joint capsules as TA98 "joint capsule" plus the TA98 joint and digit. Bursae: subcutaneous olecranon, subcutaneous acromial, subtendinous of infraspinatus, superior and inferior biceps femoris. Located by position: extensor retinaculum, ligament of head of femur, fibular collateral ligament.
  - Colon, caecum, appendix, rectum and oesophagus wall layers: serosa, subserosa, submucosa, muscularis mucosae, circular and longitudinal layers, mesocolic/omental/free taenia, adventitia. Lamina propria and epithelium are filed under mucosa.
  - Breast: areolar tubercles, lactiferous duct.
- Visible parts without any TA98 term: 1 (an unidentified "Fan ligament"). Exact terms: joints/ligaments 64 → 93, bones 39 → 48 (most bone instances are class terms or groups), digestive 68 → 82. 1,700 distinct exact TA98 terms over visible parts.

## registered collection female connective tissue and fascia

- Source: registered collection `human_20_female` connective_tissue and fascia systems. There are 1773 parts in all. 1041 are installed after removing articular cartilages (TA98 has no term for them), parts with no TA98 term, and native reference body duplicates (same exact term and side within 3 cm).
- TA98 labelling:
  - Automatic authority match, then `bdc_ta98_refine.py`.
  - Curated rules cover joint capsules, vertebral, costal and costotransverse ligaments, palmar and plantar plates, digital pulleys (anular and cruciform parts of the fibrous sheath), collateral ligaments by joint, the tibial and fibular collateral ligaments, and the check ligament of lateral rectus.
  - When one TA98 class term covers many structures, the label is the TA term followed by the source location (for example, "Joint capsule, interphalangeal proximal of fifth toe (right)").
  - When the same English term appears in different TA98 entries, the TA98 parent qualifies it (for example, "Articular disc of temporomandibular joint").
- Fit:
  - Transferred by the bone/organ anchor field.
  - Then a bone/skin ARAP sticky solve (`bdc_fit.py`). Faces inside bone went from 21.6% to 3.8%; faces outside the skin went from 5.4% to 0.1%.
- Materials by tissue: ligament, cartilage, fascia (translucent), serosa (pleura, peritoneum, pericardium; translucent), meninges.
- Viewer: ligaments appear under Skeletal › Joints and ligaments, fasciae under Muscular, pleura under Respiratory, peritoneum under Digestive, meninges under Nervous. Thoracic ligaments are placed in the Torso region.
- Coverage of exact TA98 terms after this pass: A03 136/383 (was 93), A04 291/654 (was 242), A05 112, A06 87, A10 18, A14 349.
- Geometric duplicate check (`bdc_dedup.py`, run after a full install with `CDEDUP_OFF=1`):
  - A candidate pair has the same exact TA98 term or base name, and the same side.
  - It counts as a duplicate when ≥60% of one surface lies within 5 mm of the other, or when the names match and the bounding boxes overlap ≥60% of the smaller box.
  - Results: 148 registered-collection parts dropped (the native is kept) and 49 numbered native names (e.g. "[III]") replaced by the joint name.
  - The native interspinous ligament mesh is hidden; the registered-collection per-level segments cover it.
- Native mislabel: the reference body "Calcaneofibular ligament" 460/456 is a 27 cm muscle-textured strip running from the heel up the calf. It is the calcaneal tendon, which the registered-collection muscles already provide, so it is now hidden. The real calcaneofibular ligament is the small native part at the lateral malleolus.
- The native patellar ligament and interosseous membranes moved from the muscular system to connective, with an untinted material. They now sit under Joints and ligaments.

## Labels (TA98 naming of all registered-collection parts and native numbered pieces)

- registered-collection vessels, nerves, lymph and heart (3380 parts) are now labelled from TA98, as the muscles and connective tissue already were:
  - The label is the TA98 English name plus the side.
  - An English name shared by several TA98 entries gets "of <TA98 parent>" (e.g. "Ureteric branches of renal artery").
  - One TA98 term covering several structures, or a parent filing, gets "<TA98 term>, <qualifier>". The qualifier is the spinal segment ("Anterior root, L3"), else the source words not in the TA98 label, with a series number as an ordinal ("Posterior intercostal arteries, fifth").
  - Labels that still collide take the longer qualifier.
- `bd_ta98_upgrade.py`: 28 registered-collection names that were only filed under a parent now have an exact TA98 term. The TA98 name must appear as a contiguous phrase in the source name, the vessel/nerve kind must match, and the term must lie under the filing or name its whole TA98 parent (e.g. pharyngeal branch of vagus nerve, renal ganglia, accessory saphenous vein, rami communicantes).
- `native_names.py` names 418 native reference body pieces that only carried running numbers, using the bones they touch:
  - Metacarpals, metatarsals, ribs and vertebrae with ordinals and TA numerals ("Third metacarpal [III]", "Rib, sixth [VI]", "Atlas [C I]", "Thoracic vertebra [T XI]").
  - Phalanges get their digit ("Proximal phalanx of hand, ring finger").
  - Ligaments and capsules get their digit and joint level ("Joint capsule, interphalangeal joints of foot, fourth toe, distal"), or else their two bones ("Dorsal intercarpal ligaments, trapezoid–capitate").
  - Ethmoidal cells are anterior, middle or posterior by position.
  - "<group>: <structure>" labels become the TA98 structure, or "<group>, <structure>".
- Final clean-up: Latin wall layers become English ("Longitudinal layer of muscular layer of oesophagus"), capsules read "Joint capsule, <joint>", and descriptions name related parts by their current labels.
- Order: install with `CDEDUP_OFF=1`, then `bdc_dedup.py`, `native_names.py`, `bd_ta98_upgrade.py`, then the final install (in `bdf_chain.sh`).
- Lymphatic labels:
  - Generic node groups carry their TA98 group up the chain, stopping at regional classes: "Apical nodes of axillary lymph nodes", "Superior deep nodes of lateral cervical nodes", "Proximal node of deep inguinal nodes", "Paratracheal nodes of deep nodes of anterior cervical nodes".
  - Lymph vessel qualifiers drop the repeated "lymphatic": "Deep lymph vessel, anterior tibial".
  - Description relabelling is case-insensitive, so "Accompanies:" lists current labels too.

## Brain (HRA Allen atlas fitted into the registered collection female brain)

- Base: the registered collection female brain, already fitted to the reference body skull. Added detail: the HRA Allen brain (141 structures per side), which was hidden and sat 6–9 cm low and posterior in the HRA frame with left and right mirrored.
- `brain_reg.py`:
  - Mirror into this frame.
  - Affine fit from 63 matched structure centroids: ventricles, basal nuclei, thalamus, hippocampus, amygdala, colliculi, commissures, brainstem, cerebellum and the registered-collection gyri. Mean error went from 83 mm to 5.7 mm; volume scale 0.60.
  - Three rounds of thin-plate warping (centroids, plus the Allen outer cortex pulled to the registered-collection surface). Median surface gap 2 mm.
- `brain_fit.py`:
  - Allen subdivisions are fitted as sets into their registered-collection structure (affine ICP, boundary snap, containment) and fill it to 0.6–2.8 mm median. The subdivided structures are thalamus (14 nuclei), caudate (head/body/tail), hippocampus, amygdala (10), hypothalamus (regions), lateral ventricle (horns, central part, atrium), cerebellum, pons, medulla and midbrain.
  - Allen-only interior structures are pushed out of the registered-collection structures they cross. These are accumbens, claustrum, internal pallidum, red nucleus, substantia nigra, subthalamic nucleus, zona incerta, septal nuclei, bed nucleus of the stria terminalis, basal forebrain, cerebellar nuclei, cerebellar peduncles, inferior olive and the mammillothalamic tract (which enters the thalamus).
  - Where registered-collection has the structure itself (putamen, lateral pallidum, colliculi, geniculate bodies, pineal, fornix, commissures, third and fourth ventricles, aqueduct, tracts), registered-collection is kept.
  - Left out because their HRA position contradicts their name: frontal pole (behind the putamen), posteroventral putamen (in the frontal lobe), frontomarginal gyrus.
  - Known overlaps: claustrum against the registered-collection insula and bed nucleus against the internal capsule (17–28%); registered-collection leaves no gap there.
- `brain_build.py`:
  - Each registered-collection container becomes a group concept over its Allen pieces.
  - The registered-collection hemisphere surface not lying under a registered-collection gyrus is cut into the Allen parcels registered-collection lacks: precuneus, superior parietal lobule, paracentral gyri, opercula, transverse temporal gyri, temporal pole, straight gyrus.
  - registered-collection occipital lobe, insula, orbital gyri and inferior frontal gyrus are cut into their Allen subdivisions: cuneus, lingual, occipital gyri, occipital pole, long and short insular gyri, orbital gyri, opercular and triangular parts.
  - TA98 mapping is curated in `allen_ta98.json`.
- Visible exact A14.1 terms: 49 → 131.
- Order: install without `BRAIN`, then `brain_reg.py`, `brain_fit.py`, `brain_build.py`, then install with `BRAIN=1`.
- Brain appearance: registered-collection CNS meshes carry the registered collection nervous-system colour atlas (`bdn_uv.py`).
  - BDA v18 stream order is positions, UVs, normals.
  - The decoded atlas comes from the registered collection brain capture (male state). Female UVs coincide except one block (u ≥ 0.75, 0.5 ≤ v < 0.75), which is shifted by (−0.75, +0.25).
  - Meshes sampling mostly atlas background (eye parts, choroid plexus) stay untextured.
  - Allen pieces inside a textured registered-collection structure take its texture, one sample per face.
- Viewer tree: several qualified entries on one exact TA98 term are grouped under a node named after that term (e.g. Orbital gyri › medial, lateral, anterior intermediate, posterior intermediate).
- Brain surfaces: each registered-collection structure that Allen subdivides keeps its own closed, textured registered-collection surface, cut into its TA98 subdivisions by the nearest fitted Allen structure, as for the cortex. Allen structures that reach no part of the surface stay as volumes inside, coloured with the mean atlas colour of their structure. This replaced Allen meshes standing in for the registered-collection surface, which left holes in the medulla and on the cerebellar edges.
- TA98 intermediate groups (brain, `BRAIN=1`): every TA98 ancestor (A14.1) with at least two modelled children becomes a group concept with that exact TA98 term, per side and for both sides (114 groups). Examples: ventrobasal complex, ventral/medial/intralaminar/dorsal nuclei of thalamus, grey substance of thalamus, lentiform nucleus, basal nuclei, insular gyri, cerebellar peduncles, tectum and tegmentum of midbrain, diencephalon … brain. Labels and the inspector open level by level ("Part of", "Contains"). Generic TA98 names take their parent ("Grey substance of thalamus").

## TA98 intermediate groups, all systems

- The installer (bdf-circ/scripts/bdf_install.py) makes one group concept for every TA98 ancestor that has at least 2 modelled children. Each group is made per side and for both sides together. TA98 A01 body regions are left out.
- Current count: 1182 groups (for example: Elbow joint › Synovial joints of free upper limb › Joints of upper limb; Pancreatic nodes › Visceral lymph nodes of abdominal lymph nodes › Abdominal lymph nodes).
- Names that TA98 reuses for several ids, and bare names (muscles, fascia, phalanges, trunks, ligaments, grey/white/basal substance), take their parent's name: "Posterior rami of cervical nerves [C1-C8]", "Cerebral part of internal carotid artery".
- A group whose qualified name would equal its parent's name is skipped.
- The TA98 parent table lists the sacral anterior rami under the posterior rami. The group pass files them under the sacral nerves instead.

## open anatomy structures fitted into the female model (za_cns.py, ZA=1)

Run order:
1. Install with BRAIN=1 and without ZA.
2. Run scripts/za_cns.py. It reads the installed atlas and refuses to run on a ZA install.
3. Install again with ZA=1 BRAIN=1 …

What it adds (open anatomy, CC BY-SA 4.0, geometry from ta98-male.json; all of it unreviewed):

- **Brainstem nuclei, 22 parts, TA98 exact.** open anatomy brainstem registered onto the registered-collection brainstem (centroid/spread start, damped affine ICP, thin-plate snap). The registered-collection oculomotor and trochlear nuclei are landmarks (1.3–2.1 mm after the fit). Nuclei are kept inside the registered-collection brainstem and get one loop subdivision. Covered:
  - accessory oculomotor, abducens, motor facial;
  - superior/inferior salivatory, ambiguus, hypoglossal, posterior vagal;
  - vestibular, anterior/posterior cochlear.
- **Cerebellar lobules, 27 parts, TA98 exact.** The textured registered-collection cerebellum surface (upsampled) is cut by the nearest registered open anatomy lobule surface (fit median 0.68 mm). The cut is majority-smoothed, and fragments under 15% of a lobule's largest piece are merged into their neighbour. This replaces the Allen hemisphere/paravermis/vermis cut.
- **Spinal cord, 55 parts.** The open anatomy cord section (one 5 mm slab, right half, mirrored) is swept along the registered-collection cord and scaled per level to the registered-collection section. Tracts, grey columns, intermediolateral/intermediomedial nuclei, nucleus proprius, spinal reticular formation and the central canal follow TA98 levels from the HRA cord segments:
  - cuneate fasciculus above T6;
  - posterior spinocerebellar tract above L2;
  - intermediolateral nucleus T1–L2.

  The section shape is the same at every level, so this is schematic. Materials: ascending blue, descending red, propriospinal pale, grey columns brown, autonomic green.
- **Ear.**
  - The tympanic membrane is carried by the open anatomy ossicles registered onto the registered-collection ossicles.
  - The registered-collection bony labyrinth is cut into cochlea (by the registered open anatomy cochlea), vestibule and semicircular canals. The canals are the thin parts of the non-cochlear labyrinth, found from local thickness.

## Skull bones, sutures, liver segments, external genitalia, nasal mucosa (bdx_*.py, BDX=1)

Chain:
1. `scripts/bdx_extract.py`: decodes the registered-collection female source meshes into bdx_raw.pkl.
2. `RAWF=bdx_raw.pkl SYSF=extra OUTF=bdx_transfer.pkl bdf_transfer.py`: moves them to the reference body female.
3. `scripts/bdx_build.py`: builds the parts (bdx_parts.pkl).
4. Install with `BDX=1`.

It reads the reference body source atlas (ta98-female-runtime.json), so it does not depend on the installed state. Shared fitting helpers are in `scripts/fitlib.py`.

- **Cranium.** The textured reference body cranium (3 source pieces, merged) is cut into 21 bones by the nearest registered-collection female skull bone. The registered-collection skull is refined onto the cranium with a median distance of 2.0 mm. The registered-collection maxilla and palatine are single meshes, so they are split at the midline. Adjacency uses position-welded vertices (the reference body mesh repeats vertices along UV seams).
- **Sutures, 47 parts.** Smoothed ribbons (1.5 mm wide, lifted 0.25 mm) along the boundary between the two bones of each TA98 suture: coronal, sagittal, lambdoid, squamous, sphenoid, facial and palatine sutures. The cranial-relations limitation text about the combined cranium no longer applies; the TMJ/C0–C1 relations keyed on the old cranium ids are not shown.
- **Liver.** Couinaud segments I–VIII are computed on the registered-collection liver with its own vessels:
  - vertical planes through the middle and right hepatic veins;
  - the umbilical fissure between the segment III and IV branches;
  - the left hepatic vein plane between segments II and III;
  - the horizontal portal plane;
  - segment I as the caudate wedge behind the hilum.

  registered-collection shares: I 2%, II 17%, III 13%, IV 21%, V 13%, VI 14%, VII 11%, VIII 9%. The labels are then carried to the reference body liver by nearest face (registered-collection liver refined onto it, median 4 mm). The reference body left lobe takes a larger share.
- **External genitalia.** registered-collection glans and crura of the clitoris, bulbs of the vestibule and labia minora (split at the midline), moved by the BD→reference body vagina fit (median 2.2 mm). The bulbs and labia sit 3–5 mm from the vaginal opening; the glans is 16 mm under the skin.
- **Nasal mucosa.** registered-collection mucosa moved with the skull warp.
- **Not built.** Mediastinal divisions, peritoneal recesses and fossae, mons pubis and labia majora are spaces or skin regions with no source mesh, so no geometry was invented for them.

## Gap fill (gapfill.py, GAP=1) and the single folder

**Gap fill.** Female meshes are cut by registered source parts; new parts are added only for TA98 terms the model lacks.
- **Lungs:** cut into lobes by the registered collection (BD) female lobes, then into bronchopulmonary segments by the open anatomy segments inside each lobe. The left lung has 8 segments (open anatomy has no left S VII).
- **Organs cut by registered-collection parts:**
  - sternum into manubrium, body and xiphoid;
  - pharynx into naso-, oro- and laryngopharynx;
  - the duct mesh into bile, pancreatic and accessory pancreatic ducts;
  - parathyroids into the 4 glands.
- **Organs cut by open anatomy parts:**
  - thyroid into its 2 lobes, isthmus and pyramidal lobe;
  - hyoid into body and greater/lesser horns;
  - gallbladder into fundus, body, infundibulum and neck;
  - the registered-collection heart valves into 11 cusps.
- **registered-collection additions:**
  - larynx: vocal and vestibular folds, laryngeal vestibule, cuneiform and triticeal cartilages, trachealis (refined on the thyroid cartilage);
  - pericardium: fibrous and serous layers;
  - fat: buccal fat pad, retrobulbar fat, perinephric and paranephric fat;
  - eyebrows and tarsal glands;
  - ducts: hepatic ducts and their branches, cystic duct, parotid and submandibular ducts.
- **open anatomy additions (66 parts):** bursae, intermuscular septa, intercostal and interosseous membranes, fibular head ligaments, the cauda equina, the nasolacrimal ducts, the infrapatellar fat pad, the interpubic disc and others.
  - Placed by blended similarities of the 29 nearest matched bones (median bone fit 2.4 mm).
  - The tarsi are fitted onto the registered-collection tarsal glands of the same lid.
- **Coverage after this pass:** 2279 distinct TA98 structures, 2711 with the TA98 group levels, in 8046 parts. The largest open areas are bone surface landmarks (A02), CNS nuclei and tracts above the medulla (A14), body regions (A01) and sense-organ details (A15).

**Single folder.** `.local-models/female-ta98-complete/` holds:
- `atlas.json`;
- `geometry/` (222 used chunks, hard links);
- `textures/` (58 files);
- `manifest.json`;
- `README.md`.

Open it with `/?model=local-female-ta98&atlasSource=/local-models/female-ta98-complete/atlas.json`. Build everything with `zsh .local-models/geometry-workbench/bdf-circ/scripts/build_all.sh`. It runs:
1. bdf_chain.sh, ending with the base install;
2. bdx_extract and transfer;
3. za_cns, bdx_build and gapfill;
4. the final install;
5. package_model.

Steps 3 and 4 depend on the base install, so the order matters.

## Bone parts and landmarks (bone_landmarks.py, BONES=1)

- **Bones covered: 85.** Femur, tibia, patella, humerus, radius, ulna, clavicle, scapula, hip bone, talus and calcaneus on both sides; sacrum; atlas, axis and C3–L5; 24 ribs; mandible; and the cut skull bones (frontal, parietal, temporal, sphenoid, occipital, maxilla, palatine, zygomatic). Each textured female bone is cut into its TA98 parts, processes, surfaces, borders, lines and crests: 793 parts, 267 distinct terms.
- **Method.** No labelled template source exists, so position rules are used in each bone's frame:
  - long bones: the level along the shaft axis, and position within each cross-section in local radii (lateral/anterior);
  - other bones: box fractions, distance from the midline plane and surface normals;
  - `spot()` patches around extreme points for tubercles, spines and epicondyles.
- **Naming.** Rules name TA98 terms inside the bone's own TA98 subtree, so ids cannot drift. Vertebra and rib parts carry their bone ("Spinous process of thoracic vertebra [T VI]"). The rest of a bone keeps its parent term (for example "Shaft of femur").
- **Not cut.** Foramina, canals and internal grooves have no surface of their own here.
- **Review.** `BONECOLOR=1` builds a colour-per-part review version. A point-splat contact sheet (front, both sides, back) was used to check the femur, tibia, humerus, radius, ulna, clavicle, scapula, hip bone, sacrum, L3, T6, a rib, the mandible, the temporal, sphenoid, occipital and frontal bones, the maxilla and the talus. Boundaries are approximate and unreviewed.
- **Installer.** A later bundle may now replace a part made by an earlier bundle (the skull bones from bdx_build). A group with `concept: False` reuses the replaced part's own concept.
- **Coverage.** Skeletal chapter (A02): 381/981 terms (was 103). Model: 2531 distinct TA98 structures, 2989 with group levels, 8754 parts.

## Skin containment (skin_contain.py, after install)

- **Removed:** the registered-collection eyebrows, a thick second pair over the skin's own brows.
- **What it does:** after the install, every visible non-skin part that crosses the skin, or lies within 1.5 mm under it (and so shows through), is placed 1.5 mm under the skin. The correction spreads smoothly over the rest of the part. Only each part's own vertices count, since some parts share a vertex buffer.
  - Installer-written parts are rewritten in place.
  - reference body source parts get a corrected copy appended to the model's own geometry.bin; the source files stay untouched.
- **Excluded:**
  - skin, hair, lashes, nails, nipples and areolae;
  - eyes, teeth, lips and mucosa;
  - the openings: anus, urethra, vagina, nostrils, external acoustic meatus.
- **Ear cartilage:** the auricular cartilage is first fitted onto the skin of the auricle (damped affine ICP), then placed inside it.
- **Runs:** run it 2–3 times, until it converges (768, then 223 parts moved). In thin places (fingers, ear, nose) a part cannot be 1.5 mm from both skin surfaces, so those parts keep moving slightly on every pass.
