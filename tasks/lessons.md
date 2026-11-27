
## Female circulation retarget (2026-11)
- Read instructions literally: "use the registered collection heart" means the BD heart, not the female reference body heart.
- When importing a replacement set (BD vessels), hide the native set it replaces in the same install step, then count visible duplicates (by normalised name) before saying done.
- Always end with the viewer link; keep replies short.
- Viewer links must show only the system under review: append `&layers=cardiac,arterial,venous` for circulation (`layers=` takes system ids).
- "Interrupted vessels" means every visible branch must physically meet its parent. Check source-contact joints, free component ends, and native-to-BD attachments (female-specific vessels to their real BD trunks), not just BD↔BD contacts. Native leftovers with no BD connection read as breaks, so hide them.
- Before retargeting a male source onto the female, check for a native female source. The registered collection female capture (fetch/content/registered collection/female-complete, packed in fetch/viewer/data/female) existed and has the complete female vasculature.
- When a transfer blends per-anchor rotations, thin tubes made of long triangles fold. Average the displacement over a few-mm source-space ball, so each cross-section moves as one.
- A catalogue install must also add concepts for the new parts. Viewer search lists only concepts; concepts whose elements are all hidden are dropped, so replaced parts become unsearchable.
- Use exact parent/ancestor equality for terminology matching. Substring matching mapped lumbar branches to intercostal terms and "callosomarginal" to "marginal artery".

## Eyes: judge the open-eye look against real proportions, in millimetres
- Pattern: I declared the eyes good when the iris became visible. But the skin opening was 28.5 × 12.6 mm and showed the whole iris, and the user flagged that the openings were too large.
- Rule: when fitting eyes or eyelids, measure the palpebral fissure in mm against the iris. The target is about 28 × 10 mm, with the upper lid covering 1–2 mm of the iris and the lower lid at the lower limbus. Check it straight on and from an oblique view before saying it's done.
- Rule: skin openings that are shader cut-outs in UV space move with the mesh. Narrow them by moving lid geometry, and keep the UVs. Any visibility test must treat cut-out skin as see-through.
- Before reworking a low-resolution web-model surface, look for the original high-resolution source (here `.local-models/models/Female_2024.glb` + `textures1/` 2048 px). Transfer its geometry and texture rather than tuning cut-outs on a 512 px paint.
- Skin materials with neck-tone blending read vertex height: keep skin positions in world coordinates (no `geometryOrigin` recentring).
- A part the user may isolate must look right on its own: merge replacement geometry into that part (shared texture atlas) instead of adding an overlay part.
- Lashes and other strand geometry: judge each strand as a whole (weld triangle-strip fragments first) and test it against the eyeball; never move strands by per-vertex fields or root slides. Both made them stray.

## Imported muscles looked malformed
- Pattern: I reported good aggregate fit numbers, but 237 muscles had silently received no fit target and were only carried along by neighbours; the user saw a lumpy abdomen.
- Rule: after any matching step, list what did not match, especially large structures, before reporting quality. Aggregate distance over matched items hides the unmatched ones.
- Rule: for thin sheets, use normal-aware (point-to-plane, same-facing) correspondences. Check for inside-out meshes (face-normal majority against the rigidly aligned source), not only creases.
- Rule: where an imported structure is one-to-one the same as a good native mesh, show the native mesh and keep the term. Re-meshing an equivalent adds artefacts, not information.

## Never invent anatomical labels
- Pattern: I labelled split parts "Rectus abdominis, belly 1–4". TA98 has no such terms; the user objected ("use TA!!!").
- Rule: every anatomical label comes from TA98: the English term plus side. Subdivision terms take their TA98 parent ("Long head of triceps brachii"). A structure TA98 doesn't name is labelled by the TA98 term it is filed under, plus the source qualifier. Never split geometry into entities that TA98 doesn't name.
- Re-run duplicate selection after any TA98 refile: a wrong automatic filing (medial/lateral collateral → "collateral ligaments") hid real native duplicates (tibial/fibular collateral) until refiled.
- Duplicate checks must be geometric (point-to-surface on dense samples), not name/centroid only, and must run on a full (un-deduplicated) install. Also audit native meshes for implausible size/material vs. their label (reference body "calcaneofibular ligament" was the calcaneal tendon).
- Label passes must cover every imported system, not only the one being worked on; audit all visible labels (digits, brackets, colons, source names in descriptions) after each pass.
- Scripts that derive names from installed labels must read a clean install (CDEDUP_OFF=1), never their own output.
- Check imported atlas labels against position before trusting them (HRA Allen "frontal pole" lay behind the putamen); exclude contradicting labels instead of mapping them.
- Verify decoded UVs by triangle coherence (UV edge vs 3D edge), not by bounding boxes: the BDA UV stream sits between positions and normals, and bounding boxes looked right with the wrong offset.
- geometry_tools runtime (run.sh) caps a run at 300 s; asking for more fails at once with a ValueError, not with a timeout.
- TA98 English names repeat across ids (posterior rami, basal part, deep nodes). Any generated label must be qualified by its parent when its name collides; check by dry-running the naming against the installed atlas before reinstalling.
- Registration between two body models (open anatomy male to BD female) needs an initial centroid + spread alignment; damped ICP alone from a 10 cm offset converged to a 13 mm misfit.
- Labelling a target surface from source parts: use the distance to each part's surface, not the nearest source vertex (vertex density biases toward dense meshes; all vermis lobules vanished).
- BD meshes can be open (spinal cord tube): winding-number inside tests then report 0; check extents per section instead.
- Scripts that read the installed atlas and replace parts must refuse to run on their own output (the labyrinth disappeared after reading a ZA install).
- reference body meshes repeat vertices along UV seams (cranium: 68% of edges open). Weld by position before any adjacency work (smoothing, islands, boundaries).
- Region labels on a coarse mesh: build boundary features (sutures) as ribbons along the welded boundary edges, not as faces within a distance (large faces are skipped or make jagged strips).
- Territory cuts from short vessel stubs: nearest-seed Voronoi favours the central stubs. Use anatomical planes (hepatic veins, portal plane) and compute in the source frame where vessels and organ agree.
- zsh does not word-split an unquoted `$E`: `env $E cmd` passes one variable with spaces, and the install silently dropped the nerve, lymph and connective layers. Write the variables out inline, or use `${=E}`.
- Scripts that cut parts must read the base install, not their own output. Guard with an assert, and keep the order in build_all.sh.
- Packaging: rewrite every `/local-models/` string in the atlas (chunks, maps, height maps, hair fields), not only material maps; then check that no external URL is left.
- Region cuts that share one texture are invisible in the viewer. Review them with a colour-per-part build plus an offline point-splat contact sheet (numpy + sharp); the viewer search and UI are too slow for that.
- Select source meshes by their stable source ids (reference body LOCAL ids), not by installed names: an install that replaces them makes name lookups silently miss bones (85 → 39).
- Being inside the skin is not enough: parts within about 1 mm under it show through in the viewer. Keep a clearance of 1.5 mm.
- reference body parts can share one vertex buffer: geometry edits must use only the vertices the part's own faces reference.
- Saved view links encode visibility by part order: removing parts shifts what an old link hides.
- Checking only vertices against the skin is not enough: large triangles cut through concave folds (gluteal fold, ear concha) between inside vertices. Test points inside each face too.
- Name filters must match whole words: 'retina' caught every retinaculum and 'anal' caught every canal, so these parts were silently never processed.
- Never call a fix done from numbers alone: check it in the viewer at depth 0, close up, with a fresh load. The viewer caches chunks in IndexedDB by URL.
- Closing the skin with the sclera for the inside test makes the eyeball interior count as "outside". Exclude the eye's internal parts.
- The camera aspect belongs to the canvas: a saved pose (view links) must never restore it, or the body stretches in a window of another shape.

## Male TA98: parts poking through the skin
- The drawn male skin is the default-visible ZA skin patches (single layer, ~89 L), not FJ2810. FJ2810 is hidden and is a 2 mm double-walled shell, so inside tests against it are wrong.
- Inside/outside tests on the patch union: nearest-face sign near the surface, and a point counts as "outside" only if the 6-ray parity vote / winding number agrees. Either test alone misfires at seams, ear canals and eye openings.
- Parts that are mostly outside by >10 mm are misregistered, not poking (DMTM orbicularis oculi from a Zygote source ~55 mm too high, a lymph node at the world origin): hide them, don't deform them.
- Pipeline: scripts/stage-male-ta98-skin-clearance.py (dmtm-v8 → dmtm-v8-skin-clearance-v1) → npm run prepare:male-ta98 → server/verify-male-ta98.mjs → rsync outputs/male-ta98-pages/models/male-ta98 to public/models/male-ta98. Audit with scripts/audit-male-ta98-skin-protrusion.py.

## Male TA98: duplicate structures
- "Two spinal cords" was one misregistered LOCAL cord running behind the ZA spinal dura; the user wants the duplicate removed: hide it (DUPLICATES in `scripts/stage-male-ta98-skin-clearance.py`), never refit or relabel it.
- Check `suppressed` as well as `defaultHidden` before calling a part visible; the viewer hides suppressed parts.
- Not duplicates: wall layers stacked <1 mm apart, DMTA regions (0.15 mm overlay shells on their bone), artery/vein pairs on one course, shared L/R sutures. Real duplicates were O3M copies of ZA ligaments and eye shells; hide them via `DUPLICATES` in the stage script. Audit: `scripts/audit-male-ta98-duplicates.py`.
- The ATTRIBUTION.md `atlas-source-references` table is not required; when a build input like that is missing, don't block on it, make the build tolerate it.

## Male TA98: malformed vessels
- Topology/gap audits miss shape defects. BDM vessels are 6-sided tubes that sometimes collapse into flat blades at bends; measure ring cross-sections (geodesic level sets, width vs depth), not just manifoldness.
- Find defects yourself: render every structure (contact sheets, then true-scale renders of suspects) and review all of it. Never ask the user for examples.
- Rounding a stretch needs a seam-cut harmonic angle (plain (cos,sin) interpolation collapses to one side) and a move guard (no vertex moves further than the flattened width). Side-mouth rings (where branches join) stay fixed.
- A second defect class: stray vertices (single vertices or a whole end ring sitting 20-50 mm away, joined only by long edges) make wedge fins and spikes; found around both kidneys. `fix_strays` in `scripts/vessel_sections.py` returns them before rounding stretches.
- One-sided fins/bulges can be only ~3:1 wide:deep; use FLAT_ASPECT 2 and measure the guard by the ring's real extent, not 2x its standard deviation.
- Verify vessel fixes with offline renders (scratch venv with matplotlib) when the viewer state is hard to reproduce.
- Find where a malformation came from before patching it: walk the part back through the
  atlas history. The renal/hepatic blades were round until `thoracic-skeletal-refit-v17`
  squashed them while moving them 40-50 mm; 799 vessels carry such squash. Restoring the
  section from the installed shape (same topology) along the current course fixes the class.
- Restore only the offset across the section; the position along the vessel belongs to
  the current course, or refitted (stretched) vessels shrink back.
- Pin open ends (where other vessels join) and fade any reshaping in from them.
- Check a repair's own output for new defects: on coarse meshes single-vertex "rings" made
  the median ring width 0, every ring looked flat and a dural sinus was rebuilt as a thread.
  Statistics must ignore degenerate samples, and a rebuild must refuse a zero radius.
- Dural sinuses and other naturally flat channels are not flattened tubes.
- Every automatic reshape needs an acceptance test against ground truth, applied per part:
  `section_error` (local section radius vs the installed shape) flagged all four visual
  regressions of a whole-vessel ARAP re-pose and all improvements. Candidates that do not
  lower it are rejected; triangle flips are backed off with `keep_orientation`.
- Re-render before/after for every fix round; a fix that helps the target can damage others.
- Before any constrained solve (ARAP, Laplacian), check that every connected piece keeps at least one held vertex; a piece that is wholly free has no anchor and can drift arbitrarily (a spinal branch piece moved 1.2 m). Section-shape metrics cannot see such a jump.
- Measure local rotation over a distance tied to the mesh spacing (longest incident edge), not only the radius: on coarse thin vessels a ball of a few radii holds a single ring, and folds spread over several edges never exceed a per-edge threshold.
- When blending a push field, give each offender its full push and fade it smoothly (weighted mean × max kernel weight), with a fade length tied to the push depth. Diluting the push leaves residues that a final per-vertex snap turns into steps and hooks.
- Never reshape one vessel at a time without holding its joins. A join is an open end meeting another vessel, or another vessel's open end meeting this one; parallel contact is not a join and holding it blocks repairs. Clear all visible vessels in one shared push field so joined vessels move together. After any vessel change, scan join gaps against the input (scratchpad gaps.py): per-vessel steps had opened 753 joins.

## Smooth continuity at vessel joins (2026-11-26)
- A shared, coincident ring is not enough: the user sees a step when one vessel is much wider just past the ring (axillary 4 mm ring → basilic bulb 6.5 mm). Check the radius profile across each end-to-end join, not only that ends meet.
- Fix by widening the narrow side to a 15° envelope of both surfaces (30° was too steep: the ledge stayed visible), per angular sector (keeps oval / one-sided confluences); never squash the wide side or push along the axis (both looked worse or did nothing).
- Gate hard and render: mouths in a wall, forks (another like-sized ring of the same system within 3 radii — on ANY vessel, incl. the same one), and artery↔vein coincidences must be skipped, or aortic branch bases and the SVC/aortic bifurcation grow fins. Move nested (double) walls by the outer wall's amount, and give every copy of a shared point one common move, repeated along chains; then check flipped triangles = ~0 and separated coincident pairs = 0.

### Join seams, round 2 (2026-11-26)
- Source end rings of joined vessels are 0.001–0.1 mm apart, not welded: pair ring points by proximity (5% of ring radius), never by exact coordinates, or the two sides move apart and leave a step.
- A move received at a shared point must not fade into vertices the taper already sized (the wide side), or it inflates them into a shoulder; give every copy of a point the largest move, not the mean.
- Shading: copies of a ring point in two meshes must carry one shared normal; measure cross-vessel normal mismatch (was a median of 19°).
- Check in the viewer at close zoom with the neighbouring vessels shown; isolating a pair makes inflow openings look like notches.
- Branch mouths: a cluster of ring-point copies (branch end plus two aorta copies) gets ONE normal, the sum of those facing the cluster's way; per-copy averages differ and leave a seam (707 → 97 mismatched points).

### Lone-ring vessels (2026-11-26)
- Coarse thin vessels (rings ~6 mm apart) put each ring alone in its geodesic ball, so `restore_sections` treated it as a flat sheet and skipped it; the re-pose's anisotropic stretch (round 0.5×0.5 → 0.95×0.52) survived as flat blades. A flat single ring is now handled: in-plane squash ratio, axis = plane normal. Section error of left posterior intercostal vein 10: 1.70 → 0.81; restored vessels 291 → 509.
- Always open the user's pasted link with the exact hash (camera included) before theorising; my offline renders missed that it read as a blade from their angle.

### Sweep, don't wait for links (2026-11-26)
- The user should never have to hand me examples. After a fix, sweep every vessel myself: `scratchpad/sweep.py` renders 48 shaded thumbnails per sheet, and `rsheet.py` renders before/after for whatever a repair step touched. Judge the before/after sheet of each repair, not just numbers.
- `reshape_small` (whole-vessel ARAP from the source shape, joined ends held) removes node fans, wedges and corners that no per-ring step can; it was capped at 150 vertices. Cap now 4000 with a 20 mm max-move guard: 70 vessels re-posed, joins still closed.
- Thumbnails scaled to the whole vessel hide fine node-level defects; to find those, test candidate defect classes against a known bad case first (6a0 was a good one) before trusting a detector.

### Coarse source vessels (2026-11-26)
- Flat-shaded sweep sheets make every coarse hexagonal tube look faceted; judge from the viewer's own normals (`scratchpad/cs/`, vertex count < 100) and then in the viewer.
- The real remaining defects are source meshes of 12–70 vertices: pinched or twisted waists (left common facial vein), collapsed sheets (right truncus anterior), spiked or toothed ends (superior gluteal deep branches, thoracoacromial artery), lumps (coeliac trunk). Most are branched pieces (3+ boundary loops, Euler number < 0).
- Tried and dropped: Taubin smoothing (a 22–70 vertex mesh slides along the tube, moves > 1 radius, roughness barely falls) and a whole-tube rebuild between the two end rings (only plain two-ring tubes qualify, gains are marginal, some come out twisted). A real fix needs a remesh, which changes vertex and index counts and so breaks the "source arrays retained / topology unchanged" claims in the atlas labels; do not do it silently.
- `gaps.py` against the current build shows 3 opened joins; the 753 in an old output file were from a stale build. Re-run before believing an old log.

### Remeshing coarse vessels (2026-11-26, approved by the user)
- `scripts/remesh-male-ta98-coarse-vessels.py` (+ `vessel_remesh.py`) runs BEFORE the stage script: the stage now reads `dmtm-v8-remeshed-v1`. 138 crumpled coarse vessels (≤250 vertices, no UVs) are regenerated as round tubes (trunk = two widest ends, other ends as branches into the trunk axis, closed by a fan); end-ring vertices are kept exactly, so joined ends stay closed (checked: 2100 → 2103 joined boundary vertices).
- Selection: source surface p90 distance from its rebuilt tube between 0.3 and 2.5 section radii (lumpy but still followable); beyond that the rebuild does not follow the source (long tortuous vessels, 490 of 718), so it stays. Pieces of end rings only (no inner vertices) always qualify: their score is ~0 though they are a skew ruled strip (right truncus anterior).
- The atlas now says so: `vesselRemesh` block, `topologyUvsIdentitiesUnchanged` False, per-part `vesselRemesh` and provenance text.
- After any change: re-render before/after, check joined boundary vertices, and open the vessels in the viewer.

### Remesh follow-ups (2026-11-26, user found a broken facial vein)
- Check what the REMESHED vessels do to their neighbours, not only their own ends: an end ring of another vessel lying on the remeshed vessel's WALL must stay on it (external carotid → facial artery lifted 1.2–2.4 mm). `remesh-…py` now collects those contact points, pulls the new surface onto them, and keeps the source mesh if it cannot (10 vessels kept).
- Winding: per-face inward/outward tests flip faces next to skewed end rings, leaving inconsistent winding (80 of 128 vessels). Make faces consistent by walking edges, then flip each piece as a whole. Verify: no directed edge used twice.
- Trunk of a Y: the pair of ends furthest apart (weighted by width), not the two widest — the common facial vein stepped sideways where the facial end entered.
- Whole-model checks worth re-running after any geometry change: scratchpad `contacts.py` (end of one vessel on another's wall: still touching?) and `rings.py` (coincident end rings: still together, radii similar?). Its small ring "gaps" for cremasteric/testicular are merged-loop artefacts; per-point gaps there are ≤0.17 mm.
- Tried and reverted: a final "re-seat" of coincident source points (0.3 mm clusters chain along dense rings; lifted contacts 51 → 374 pairs). One real single-vertex lift remains: left perineal vein ↔ left posterior scrotal vein (4.8 mm), caused by the clearance push moving the two copies in opposite directions.
- The dev server (vite) does not notice replaced files in public/ after the rsync: restart it (preview_stop / preview_start) or the viewer reports "An anatomy file was incomplete".

## Tree-shaped coarse vessels (submental / inferior labial veins)
- Whole-vessel tree remesh (trunk + attached branches) scores p90 11–29 section radii against these veins and drops or mis-attaches branches (the source branch graph has loops/ends the skeleton cannot place). It is NOT usable as-is: do not raise MAX_VERTICES above 250 for it.
- A radial-regularise pass (pull vertices to the local median radius around the local axis) moves <0.9 mm and changes nothing visible.
- In the viewer the branches of these veins are clean tubes; the remaining defects are the sparse 54-face stub at the trunk end ring (41,1618,37) mm and the crumpled thin tips. A fix needs a LOCAL rebuild of just that stub, stitched to the existing mesh, not a whole-vessel remesh.
- Fix that worked: `rebuild_end_stubs` (vessel_remesh.py) cuts the stretch ~10 end radii back from each wide open end, remeshes it as a small vessel, and stitches it back (ring and outside vertices untouched). Applied only to vessels within 80 mm of (20,1615,60) mm (face/neck); roughness cannot tell defect from coarse hex tubes, so a body-wide run would hit 485 vessels. Lifted contacts after: 48 (was 51).

## TA98 label counts
- A part is labelled if it has `ta98Term` OR its concept has an entry in `byConceptMatch`/`byConcept` of app/data/ta98-metadata.json. Counting only `ta98Term` (or only atlas bindings) understated labels twice; always use the viewer's own resolution (terminologyForConcept).
- Running build-ta98-metadata.py would rewrite ~5,000 concepts (it sees the female concepts too), so new reviewed rules go in via scripts/apply-ta98-reviewed-nervous-identities.py, which appends to scripts/ta98-reviewed-mappings.json and patches the metadata in place.
- `server/ta98-catalogue-tags.test.py` fails on TA98DIAPHRAGM:A04.4.02.001 (parent vs exact) independent of the nervous labels.

## Muscular labels (male TA98)
- Only 39 visible muscular parts lacked a label; 13 got reviewed rules via scripts/apply-ta98-reviewed-muscular-identities.py. Extensor hoods and the tongue "transverse-vertical-shared" carrier have no TA98 term, so stay unlabelled.
- Check hand vs foot before mapping a shared name (abductor digiti minimi has A04.6.02.062 and A04.7.02.063).

## Connective labels (male TA98)
- 307 unlabelled connective parts: 259 concepts labelled via scripts/apply-ta98-reviewed-connective-identities.py (54 exact, 76 parent rules). Source names are often truncated ("ligame", "join", "shea"), so list both truncated and full spellings as keys. Articular cartilages, finger lateral bands, serratus bursae and frenula capsulae have no TA98 term and stay unlabelled.

## Skeletal labels (male TA98)
- Only 38 skeletal concepts lacked a label: the 16 teeth (canine exact A05.1.03.005; incisors, premolars, molars parent since TA98 does not number them) via scripts/apply-ta98-reviewed-skeletal-identities.py. Nasal septal lateral process, unresolved distal surfaces and the unclassified patellar border stay unlabelled.

## Digestive labels (male TA98)
- 51 unlabelled digestive parts, all wall layers: 46 labelled by scripts/apply-ta98-reviewed-digestive-identities.py. Epithelium and lamina propria have no TA98 term, so they sit under the organ's mucosa (parent); muscularis mucosae and submucosa in small-intestine folds are exact. Organ-unspecified strata, epithelium and "external connective coat" stay unlabelled.
- BSD sed on macOS rejects some patterns used with `&&`; use python for file edits.

## Integumentary labels (male TA98)
- 28 unlabelled skin names (56 concepts) were surface regions, nails and eyebrow hair: labelled via scripts/apply-ta98-reviewed-integumentary-identities.py using A01.2 regions and A16 terms. Only the radial foveola has no TA98 term.

## Remaining systems (male TA98)
- scripts/apply-ta98-reviewed-other-identities.py: lymphatic nodes, fascia, urinary adventitia, suprarenal capsule (15 rules). scripts/apply-ta98-reviewed-vessel-identities.py: 37 vessel names whose TA98 term is the named parent vessel (branch/division → parent) or the plural of a TA98 vein. Left unlabelled: attachments (not TA98 terms), vestibular nuclei and retropharyngeal nodes (two TA98 homes), intercostal collaterals and cutaneous branches, segmental hepatic branches of unnamed side.
- Vessel names embed left/right inside the name ("left medial circumflex femoral veins"); normalise before matching.

## Vessel branches without a unique TA98 parent
- scripts/apply-ta98-reviewed-vessel-branches.py: 107 numbered/side-neutral branch names (intercostal collaterals, cutaneous and dorsal branches, lumbar artery branches, perforator veins, segmental hepatic veins) placed as parent under the TA98 vessel family (204 concepts). Left unlabelled when TA98 has several equally valid parents (lumbar veins A12.3.07.012 vs .09.003, segmental medullary arteries, gastric branches without a side) or no matching term.

## Remaining vessels
- scripts/apply-ta98-reviewed-vessel-remaining.py: 63 more parent placements (109 concepts). Where TA98 has several homes, used their lowest common ancestor (gastric branches → coeliac trunk, straight artery → superior mesenteric artery, segmental pulmonary veins → pulmonary veins). Lumbar-vein branches placed under the IVC lumbar veins (A12.3.09.003). Radicular, segmental medullary and spinal arterial branches stay unlabelled (common ancestor is only the aorta).
- Last 9 vessels: radicular arteries, segmental medullary arteries and spinal arterial branches span the whole cord, so they sit under the aorta (parent, no code); para-arterial veins of the arm sit under vena comitans A12.0.00.031. Every visible vessel now resolves.

## Attachments
- Origins/insertions are not TA98 terms. scripts/apply-ta98-reviewed-attachment-identities.py places 699 attachment concepts as parents under their muscle's TA98 term (found via the labelled muscle of the same name, then 13 manual names). Suffixes: .ol/.or/.el/.er and numbered .o1l/.e12r variants. Common extensor/flexor tendon and "trochanteric insertion" stay unlabelled.
- Shared attachments: common extensor/flexor tendon sit under the posterior/anterior compartment of forearm (A04.6.01.007/.004); trochanteric insertion under greater trochanter A02.5.04.005. All visible attachments now resolve.

## Fused end rings (kidney segmental arteries)
- A vessel whose open end loop is a "bowtie" (vertices at two sites 20 mm apart, 2 long edges) is two mouths fused into one; the remesh kept it as a spiky 11 mm ring on a 0.75 mm artery. `fused_end_tree` in scripts/vessel_remesh.py cuts such loops at long edges, makes each arc an end, and joins all ends to their geometric-median junction with tapered tubes (4 vessels: superior, anterior superior, inferior left; inferior right). Contacts may move up to FUSED_TOUCH_MM=1.7 for these because they lay on the old fan.
- Do not name a helper `_tube` in vessel_remesh.py: it already exists and a duplicate silently broke every other remesh (TypeError → "skipped"). Check the build report for "skipped" after any change.
- Check the report for "skipped" lines and compare vessel counts before publishing.
- Do NOT run pull_to_contacts on fused-end trees: with 17-32 contacts on 0.5 mm tubes it inflated them to 3-6 mm cups/funnels. Junction = centre of the fused loop; contacts may sit up to FUSED_TOUCH_MM=4 mm off the new tubes (they touched the removed fan). Always screenshot the result after a geometry change.

- Continuity for rebuilt kidney vessel trees: every open end must be seated, not just near. Clean ends keep their exact source ring; fused/flared ends are run onto the nearest THICK neighbouring wall (max end-loop radius >= 2x the tube) and end 0.5 mm inside it. `fused_end_tree` also treats one huge flared mouth (radius > 5 mm, > 4x the thinnest loop) as an end to seat. The driver leak-checks every end ring (<= 1 mm from a neighbour surface) and keeps the source otherwise.
- `pull_to_contacts` bloat guard: if the pull moves any vertex > 1.5 section radii, keep the plain tube (looser 0.8 mm contact limit) - otherwise thin tubes turn into sheets/cups.
- Ring sizes differ between ends (6 vs 8 sides): never `reshape(ends, -1, 3)`; use `ringSizes`. Report exceptions from the fused path instead of swallowing them - a silent `except` hid this for a build.
- Killing :3016 does not restart it: call preview_stop then preview_start.
