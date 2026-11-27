# Oesophageal dimensions and placement repair

## Status

A local candidate now follows the detailed male's native hollow-wall profile and clears every tested skeletal guide. It is **not installed or anatomically accepted**: soft-tissue packing remains incomplete. No new TA98 terms are counted by this repair.

## Inventory and assembly

The complete 17-term inventory remains in `reports/ta98-male/oesophagus-parts-plan.json`. This repair changes only the existing whole oesophageal wall; it does not establish cervical/thoracic/abdominal boundaries, histological layers, microscopic glands or constricting muscles. The whole gastric assembly includes the cardiac region and all existing stomach/duodenal parts. All lower interfaces and their coordinates remain unchanged.

## Reference constraints first

Primary guides: the original detailed-male `ZA:` skeleton and fascia. The full source inventory has 633 parts; 372 nearby guides were staged and proved to have identical world vertices/triangles in the master. Source offsets are applied once. Bones and fascial surfaces supply the frame and represented compartment boundaries. Open fascia and gaps between ribs cannot establish a complete watertight available-space volume by themselves.

The earlier reconstructed wall crossed C6–T7 vertebrae/cartilage. The native detailed-male oesophagus crossed none of those bones. Its centreline differed by less than a millimetre, but the reconstructed outer section area was approximately 1.5–1.9 times the native area: an oversized circular profile had replaced the flattened native profile. Moving the centreline alone was insufficient.

## Geometry changed

Final staged candidate: `.local-models/ta98-review/gastric-reference-corridor-profile-thoracic-end-candidate.json`.

- 42 exact-edge source/current section samples, with 360 radial directions per sample.
- Inner and outer contours checked for single radial boundaries and positive wall intervals.
- A common radial map preserves the original face topology, ownership, corner UVs and texture materials while fitting the native outer wall and lumen.
- Quintic support begins at Y=1.252 m and reaches full fit at Y=1.270 m. Geometry at/below the support start, including hiatus/cardia and every stomach/duodenal part, remains bit exact.
- A source section just below the existing proximal cut constrains its dimensions; no extrapolated coarse proximal profile is accepted.
- 43,273 oesophageal vertices changed; maximum displacement 4.817 mm. No added or removed faces; whole assembly remains 408,500 triangles.

The intermediate `profile-thoracic` variant had 40 C6 crossing pairs and was rejected. Adding the measured proximal-end section removes those crossings. Earlier variants and diagnostic reports remain available; do not confuse their hashes with the final candidate.

## Verification and remaining defects

The final candidate passed independent attribute/face-ownership/lower-interface checks and has one closed genus-zero component, no open/nonmanifold edges, invalid links, duplicate/degenerate faces, winding defects or self-crossings. All 372 skeletal/fascial reference meshes were tested: no bone crossings remain; 288 diaphragmatic-fascia and 771 transversalis-fascia crossing pairs remain unchanged.

Adjacent-source checks are mandatory after the skeletal/fascial fit. They found unresolved left-atrium, tracheal/bronchial, thoracic-aortic and pulmonary-vein contacts. The native source has several of these defects itself. In particular, left-atrium crossing pairs changed from 2,546 to 2,931, so the candidate is not accepted. Triangle-pair counts describe surface intersections; they are not penetration volume or depth, and counts across differently tessellated native/reconstructed models are not equivalent quality scores.

`outputs/ta98-oesophageal-compartment-thoracic-end.png` shows twelve actual world-coordinate sections through the spine, fascia, left atrium, airway, aorta and pulmonary vein. It exposes the available corridor and remaining overlaps without inventing solid compartments from open sheets. The next fit must resolve those contacts within the skeletal/fascial constraints, then pass full tissue containment, clearance and 3D wall/lumen checks. Do not accept a shift that merely transfers an overlap to another neighbor.

The three textured reference views and both twelve-plane section figures were visually inspected. At Y=1.275–1.311 m, atrial traces still cross the oesophageal wall; aortic and airway interfaces also remain unresolved. These are explicit failure conditions, despite the successful skeletal and closed-surface checks.

## Anatomical references actually consulted

- [Rhyne et al., cervical CT mapping](https://pmc.ncbi.nlm.nih.gov/articles/PMC2200775/): article text reviewed, including methods and level-dependent esophagus–vertebral/longus-colli measurements. This supports close, variable spatial relationships, not fixed atlas dimensions. Figure requests failed; do not claim those CT images were inspected.
- [Oesophageal anatomical relationship review](https://www.ncbi.nlm.nih.gov/books/NBK482513/): relationship to the vertebral column, airway and heart. This provides relational guidance, not source-specific metric dimensions.
- [New insights into surgical anatomy](https://pmc.ncbi.nlm.nih.gov/articles/PMC5538986/): search excerpt described the peri-esophageal compartment; full-page retrieval was challenged. Do not treat unviewed figures or the inaccessible full article as visual evidence.

## Reproducible scripts

```sh
python3 scripts/prepare-ta98-placement-reference.py
python3 scripts/fit-ta98-oesophagus-reference-corridor.py --fit-profile --start-y 1.252 --full-y 1.270 --tag profile-thoracic-end
/Applications/Blender.app/Contents/MacOS/Blender --background --python scripts/validate-ta98-reference-corridor.py -- --candidate .local-models/ta98-review/gastric-reference-corridor-profile-thoracic-end-candidate.json --report reports/ta98-male/oesophagus-reference-corridor-profile-thoracic-end-validation.json
/Applications/Blender.app/Contents/MacOS/Blender --background --python scripts/check-ta98-placement-reference.py -- --candidate .local-models/ta98-review/gastric-reference-corridor-profile-thoracic-end-candidate.json --report reports/ta98-male/oesophagus-reference-corridor-profile-thoracic-end-contacts.json
/Applications/Blender.app/Contents/MacOS/Blender --background --python scripts/check-ta98-oesophagus-adjacent.py -- --candidate gastric-reference-corridor-profile-thoracic-end-candidate.json --report reports/ta98-male/oesophagus-reference-corridor-profile-thoracic-end-adjacent.json
python3 scripts/review-ta98-oesophagus-reference-sections.py --candidate gastric-reference-corridor-profile-thoracic-end-candidate.json --prefix ta98-oesophagus-thoracic-end
python3 scripts/review-ta98-oesophageal-compartment.py --candidate .local-models/ta98-review/gastric-reference-corridor-profile-thoracic-end-candidate.json --prefix ta98-oesophageal-compartment-thoracic-end
```

Use the configured scientific Python runtime when the shell Python lacks NumPy/Pillow. Preserve all staged candidates locally; source redistribution rights remain unverified. Match report/candidate SHA-256 before relying on findings. Textured Blender review: `outputs/ta98-oesophagus-thoracic-end.blend`; normal vectors are recalculated coherently across exact shared coordinates for rendering, not claimed to prove histological layers.

## Adjacent-tissue fit: rejected trial and guarded search

The `gastric-reference-corridor-soft-candidate.json` trial translated both wall contours by one common course field. Actual native bone/fascia traces were tested first, then atrial, airway, aortic and venous traces. It preserved faces, UVs, materials and lower junctions, remained one closed genus-zero component without self-intersections, and retained zero bone crossings across all 372 guides. Nevertheless, actual section review rejected it: the path entered the atrial chamber envelope. Geometry gates and reduced boundary intersections do not establish correct packing.

The revised `scripts/fit-ta98-oesophageal-soft-corridor.py` enforces a posterior atrial condition from the measured lowest-Z surface traces, alongside sampled bone exclusion. It uses a conservative outer convex hull and a 0.2 mm provisional numerical search margin. The margin is not a clinical measurement. The +/-10 mm, 0.5 mm translation grid found **zero independently bone-clear and posterior-centroid options** at Y=1.292173, 1.294173, 1.296173 and 1.298173 m. The guarded script wrote `reports/ta98-male/oesophagus-soft-posterior-10mm-search.json` and refused to stage a candidate. This is evidence about that bounded search and source geometry, not proof that correct anatomical fitting is impossible.

The next repair target is the adjacent atrial source/interface assembly, with skeleton/fascia unchanged. `scripts/prepare-ta98-atrial-packing.py` now stages the untouched textured source with its original corner UVs and records all 113 heart entries and all five left-atrium descendants before alteration. `reports/ta98-male/left-atrium-packing-plan.json` records two open boundary cycles (90 and 31 edges), source Euler characteristic -6, exact source provenance and nearest adjoining surfaces. The 31-edge cycle is within 0.177 mm of the source left ventricle; this is proximity evidence, not an accepted anatomical port label. The larger cycle remains unresolved. Neither loop was capped, and no atrial geometry was deformed. Open source surfaces do not certify a watertight available volume.

Commands:

```sh
python3 scripts/fit-ta98-oesophageal-soft-corridor.py --bound-mm 10
# Expected guarded refusal for the recorded source: report written, no candidate staged.
/Applications/Blender.app/Contents/MacOS/Blender --background --python scripts/prepare-ta98-atrial-packing.py
```

Installed master remains unchanged. No new TA98 coverage is claimed by these fits or staged source copies.


## Atrial source repair before another placement fit

The posterior-column trial (`b522f4284b2ed56133454a6bc41d03909c7e52d090e8f6052e8c89e7051c8c30`) is rejected. Matching oblique/section images show corrugation and continued crossing. The independent float64 narrow-phase checker reports seven self-crossings; atrial/oesophageal and neighboring cardiac contacts remain unresolved. The installed model is unchanged.

A separate source-only diagonal repair now removes three nearly collinear triangles by flipping their internal diagonals with their adjoining triangles. No vertex coordinates move, and source UVs, materials, 121 opening edges, component count and Euler characteristic remain exact. Local nonplanarity is 37–68 nanometres, within the measured packed float32 coordinate rounding scale. The repaired candidate SHA-256 is `aeb23bbbb9c14180cc2cf09564521cbb0672533ed4d05ff665644a6e19050d5a`; 5,105 triangles remain. Independent Blender verification finds zero self-crossings and no new edge/link/winding defects. Posterior, anterior and oblique textured renders were inspected. This is a local triangulation repair, not a completed wall, resolved packing or new TA98 term; it remains staged locally.

Scripts: `repair-ta98-atrial-slivers.py`, `validate-ta98-atrial-sliver-repair.py`. Reports: `left-atrium-sliver-repair-plan.json`, `left-atrium-sliver-repair-validation.json` under `reports/ta98-male/`. Run Blender validators with `--python-exit-code 1` so assertions produce a failed process status.


## Constrained atrial depth-compression trials

The source-only sliver repair supports conforming refinement without self-crossings. A subsequent depth-compression trial holds the measured anterior envelope fixed, preserves original opening edges and limits posterior motion to 85% of measured source depth. The trial (`5968aee7ec11f65968ca7a97245bc991ca0488f8bccd3099975c95c5d3e346a0`) has zero self-crossings, but its oesophageal overlap and visible corrugation remain; it is rejected for installation.

A separate attachment preparation now measures actual point-to-triangle distances, reducing venous anchor samples from all source vein vertices to 85/73/97/50 nearby unique vertices for the left superior/inferior and right superior/inferior pulmonary veins. This avoids treating distant vein projections as attachment evidence. Samples remain provisional and do not establish full ports. The resulting trial (`ba9fc1b405d843bc3a37338626d9126110813a391257893cd60a2331cea40143`) produces 24 self-crossings and worsens right inferior pulmonary vein contacts, so it is also rejected. The validator correctly returned exit code 1 using `--python-exit-code 1`. Matching textured oblique images and twelve sections were inspected.

The native oesophagus itself has source triangle crossings with the left atrium, airway, aorta and right inferior pulmonary vein under the independent float64 narrow-phase check; matching native coordinates alone cannot certify packing. See `oesophagus-native-and-fitted-adjacent-float64.json` (counts depend on tessellation and do not measure penetration volume). Upcoming fitting must address the actual shared cardiac/venous interface and smooth discrete surfaces together, rather than relax the source or candidate gates. No new material TA98 term, installed geometry or completed packing is claimed here.

## Current male distal junction plan and full terminology inventory

Before replacing the distal junction, bind the 17 oesophageal entries to the repository TA98 authority and inventory the current catalogue. The whole assembly is `ZA:Oesophagus` plus `TA98REGION:A05.4.01.001:gastric-junction`; the terminal junction is not a distinct TA98 organ. The new paired collar must share its lower outer/luminal coordinates with the gastric cardia and its upper coordinates with the retained oesophagus. Leave the rest of the oesophagus and the accepted female stomach unchanged. Preserve both openings; bridge wall thickness only. No new histological coverage is established by a generic paired collar.

| TA98 | Term | Parent | Representation | Current male source / missing status |
|---|---|---|---|---|
| A05.4.01.001 | oesophagus | A05.0.00.000 | whole or regional tissue | Missing independent representation |
| A05.4.01.002 | cervical part | A05.4.01.001 | whole or regional tissue | Missing independent representation |
| A05.4.01.003 | thoracic part | A05.4.01.001 | whole or regional tissue | Missing independent representation |
| A05.4.01.004 | thoracic constriction | A05.4.01.003 | localized narrowing of existing wall, not a solid plug | Missing independent representation |
| A05.4.01.005 | diaphragmatic constriction | A05.4.01.003 | localized narrowing of existing wall, not a solid plug | Missing independent representation |
| A05.4.01.006 | abdominal part | A05.4.01.001 | whole or regional tissue | Missing independent representation |
| A05.4.01.007 | serosa | A05.4.01.001 | independent histological wall layer | Missing independent representation |
| A05.4.01.008 | subserosa | A05.4.01.001 | independent histological wall layer | Missing independent representation |
| A05.4.01.009 | adventitia | A05.4.01.001 | independent histological wall layer | Missing independent representation |
| A05.4.01.010 | muscular layer | A05.4.01.001 | independent histological wall layer | Missing independent representation |
| A05.4.01.011 | crico-oesophageal tendon | A05.4.01.010 | independent supporting muscle or tendon | Missing independent representation |
| A05.4.01.012 | broncho-oesophageus | A05.4.01.010 | independent supporting muscle or tendon | Missing independent representation |
| A05.4.01.013 | pleuro-oesophageus | A05.4.01.010 | independent supporting muscle or tendon | Missing independent representation |
| A05.4.01.014 | submucosa | A05.4.01.001 | independent histological wall layer | Missing independent representation |
| A05.4.01.015 | mucosa | A05.4.01.001 | independent histological wall layer | Missing independent representation |
| A05.4.01.016 | muscularis mucosae | A05.4.01.015 | independent histological wall layer | Missing independent representation |
| A05.4.01.017 | oesophageal glands | A05.4.01.001 | microscopic glands | Missing independent representation |

A05.4.01.004/.005 are local constrictions, not lumen plugs. A05.4.01.017 remains microscopic. Current body/terminal surfaces alone do not prove the seven histological layers or their microscopic contents. Oesophageal packing defects recorded earlier remain unresolved and must not be erased by this junction work.

The initial straight and curved collars to the existing Y=1.30958962 m oesophageal terminal contour are rejected by strict finite-wall and crossing gates. The gastric rim is about 13 mm in radius and its oblique plane would require too tight a bend within the short available span. The next local trial retains the oesophagus above Y=1.32958962 m and clips its distal surface at that plane using exact shared-edge interpolation, without a lumen cap. A longer paired collar then connects this retained outer/luminal pair to the unchanged unwarped male gastric donor. This is a 20 mm extension of the reconstruction support, not a clinical length assertion. It does not move retained oesophagus, native neighbors or the accepted female anatomy. It needs source attribute transport, whole-wall seam, continuity, packing and actual render checks before any installation. Cardiac region ownership remains gastric; the reconstructed continuation belongs to the oesophagus.

2026-10-05 male gastric transfer preparation: source-face corner UVs, packed native normals, original material/provenance and exact source-face ancestry now survive terminal seam subdivision in a local runtime-format candidate. The shared inner/outer seam normal overrides replay independently; native normals away from the seam are preserved. No oesophageal placement or accepted anatomy changed; local candidate remains uninstalled pending coupled gastric compartment fit. See gastric-horizontal-native-attributes-validation.json.

## Fuller distal male course, local proposal

The TA98 inventory and assembly interfaces remain those in TA98-RECTUM-OESOPHAGUS-PLAN.md; this adds no new anatomical term. `scripts/stage-male-oesophageal-area-course.py` carries all eight shared boundaries together, preserves both source joins and the complete native upper fields, and retains the original carrier, faces and wall identities. Larger/refined trials failed packing and are not promoted. The retained course has a straight oblique passage segment and a smooth connection to the stomach. The upper-join lumen reference is 79.56 mm². Actual triangle-plane axial cuts sampled every 2 mm over the revised approach measure 57.69–107.75 mm², compared with a 12.16 mm² minimum in the previous narrow proposal. These discrete samples do not establish a continuous area bound.

Independent source-field replay, eight boundary self checks, seven finite cell topology/self checks and the actual joined stomach/oesophagus self check pass. The source diaphragm has zero outer/lumen crossings in the full changed-extent screen. Diaphragmatic fascia still has 217 exterior contacts; its lumen contacts are zero. Liver/caudate, omentum, upper lung/pleura, aortic and a caudate portal branch contacts remain. This candidate is not anatomically accepted or installed. All source context and both ordinary model catalogues stay unchanged. The local comparison is `/oesophageal-junction-review.html?sex=male&profile=area`; its optional contexts expose the remaining obstacles. Evidence is hash-bound in `reports/ta98-male/oesophageal-hiatus-area-preserved-{candidate,axial-sections,packing,self-validation,joined-validation}.json`.

## Rounded male course, 2026-10-09 local review

The organ and interface inventory remains TA98-RECTUM-OESOPHAGUS-PLAN.md; no new solid, TA98 identity, lumen cap or microscopic feature is introduced. The earlier angular proposal is retained as source `442873b6a6eb7b9d29358951de27333ce41ea898c7006e17b9ea0a31dc474ba7`. A single nonnegative clamped cubic basis carries all eight paired connector boundaries together, with both physiological joins fixed. Forty-five longitudinal samples approximate the continuous field. A small common transverse course field smooths the native upper portion; its sampled maximum displacement is 1.481 mm. The actual referenced carrier has 59,895 controls; unused original controls remain only as ancestry records and are omitted from the review renderer. No UVs or images are added. This is a provisional course, not an anatomical shape acceptance.

Independent source/basis/spatial-field replay is exact after Float32 packing. All seven local finite prism maxima remain strictly negative; all eight surfaces and seven finite wall cells pass Blender self-intersection checks, and the seven cells retain closed genus-one topology. The literal joined gastric/oesophageal macro surface also has zero detected self-intersections, one closed genus-one component and exact shared gastric coordinates. Actual triangle-plane cuts sampled every 1 mm over Y=1.273–1.424 m have a minimum lumen area of 62.746 mm², compared with 57.082 mm² in the angular source over the same samples. These samples do not establish a continuous bound. Focused whole-inner-surface checks show no diaphragm, fascia, either aortic source, right-lung or pleural crossings. Wider atlas packing still detects liver/caudate, heart, vessels/nerves and proximal source contacts; it is not accepted or installed.

The lesser omentum's oesophageal attachment is refitted to the actual new exterior using only controls 79–81, with 0.3 mm clearance; its other 86 controls, 139 source faces and remote attachments remain literal. Every oesophageal boundary is clear of the fitted sheet. There are no new neighbor or joined-macro crossing pairs. Remote source gastric contacts remain (430 pairs in the current joined-macro screen), so this does not accept the whole lesser omentum. The source and fitted geometry were rendered in the existing Blender session without reopening or saving the user's file. The ordinary male/female models, skeleton, source contexts and camera/explosion behavior remain unchanged.

Evidence: `reports/ta98-male/oesophageal-rounded-course-{candidate,replay,self-validation,joined-validation,packing}.json`, `oesophageal-rounded-axial-sections.json`, and `lesser-omentum-oesophageal-rounded-fit-validation.json` / `lesser-omentum-rounded-complete-wall.json`. Construction: `scripts/stage-male-oesophageal-rounded-carrier.py` then `stage-male-oesophageal-rounded-native-course.py`, followed by `stage-male-lesser-omentum-oesophageal-fit.py --rounded`. Local comparison: `/oesophageal-junction-review.html?sex=male&profile=omentum`; Current and Proposed preserve the review camera. Source rights remain unverified; none of these assets are authorized for public deployment.
