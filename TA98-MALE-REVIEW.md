# TA98 male combined review

## Files

- Blender project: `outputs/atlas-ta98-male-review.blend`.
- Viewer catalogue: `.local-models/ta98-male.json`.
- Every TA98 entry and source mesh correspondence: `reports/ta98-male/audit.json`.
- Spreadsheet: `reports/ta98-male/audit.csv`.
- Missing representations and containing structures to inspect: `reports/ta98-male/README.md`.
- Actual supplement mesh topology: `reports/ta98-male/geometry.json`.

Open the local viewer at <http://localhost:3016/?model=local-ta98&tree=guest:ta98>.
The model menu calls it **Male · TA98 combined review**. It is a local development model.

## What was built

The detailed male remains the primary body. Four male catalogues were compared:
detailed, reference GLBs, standard BodyParts3D, and local reference body male. Correspondences
come from the existing reviewed TA98 concordance and unique source names; source
groups are assembled only through known component correspondences and compatible
laterality. Partial groups remain marked as subsets.

| Measurement | Count |
| --- | ---: |
| TA98 entries audited | 7,474 |
| Mesh objects in combined review | 5,656 |
| Supplement mesh objects | 448 |
| Additional terms with independent source surfaces | 225 |
| Additional complete source-component groups | 36 |
| Terms with exact nomenclature correspondences, including review candidates | 1,711 |
| Exact terms in the detailed male alone | 1,450 |

The 261 extra terms comprise independent surfaces and component groups. These
counts do not certify anatomical completeness, shape fidelity, or registration.
TA98 child terms use recovered source geometry; no arbitrary cutting plane is
assigned to a child term. The replacement stomach tract additionally includes
reconstructed junction tissue, described below. The stomach fundus, body, and pyloric part, for example, were
already separate source objects and are recovered as separate selections.

## Viewer and Blender use

Unreviewed supplements are hidden on first load and Reset. Four stomach regions
(fundus, body, pyloric part and cardia) replace the suppressed original stomach.
The source-derived exterior is repaired as a solid before constructing a nested
inner lumen. The esophagus follows cross-section centers and radii from the
primary model. An oblique inlet cut preserves the fundus. The duodenal source's
tiny unintended perforation is sealed; its original distal cap is retained.

The complete wall is partitioned into six selectable pieces. Junction positions
and packed normals are shared exactly: **397 cardiac edges** and **264 pyloric
edges**. The final mesh has **412,470 triangles**, one connected surface, Euler
characteristic 2 (genus 0), no open/nonmanifold edges, and no nonadjacent triangle
intersections across the whole surface. These checks replace the insufficient
previous edge-only check. The physiological proximal esophageal mouth remains
open through the lumen, with an annular tissue rim.

## Anatomical review limits

Cardia (A05.5.01.007) now has its own selectable region. Its boundary is explicitly
provisional. The approximately 1 mm reconstructed wall is not a reconstruction
of individual TA98 wall layers. The complete [33-entry inventory](TA98-STOMACH-PARTS.md)
records missing gastric regions, layers, landmarks and microscopic structures;
they are not represented by relabeling the outer wall. Source licenses remain
unverified and the model stays local.

The previous genus-19 surface and the rejected experiments are superseded.
Section caps now combine the shared tract partitions instead of closing each
region with an artificial fan. Saved packed views survive append-only catalogue
changes, including the new cardia; reordered catalogues remain rejected.
Reconstruction inputs are local fields in `.local-models/ta98-review/`; run
`scripts/construct-ta98-tract-wall.py` with Blender to rebuild the nested wall.
`scripts/weld-ta98-stomach.py` checks the repair and source hashes before
projecting source colors and exporting partitions. Changed input geometry must
be reconstructed and revalidated, rather than reusing the existing repair.

Select an individual candidate
to inspect it, or use its eye control. Some supplements are alternatives to
parts of a primary organ. Hide the primary representation before showing these
in context; selecting a parent and all descendants can include alternatives.
The TA98 hierarchy notes when a selection contains supplements awaiting review.

In Blender:

1. **01 Primary detailed male** contains the original detailed geometry by system.
2. **02 Supplement geometry — review before use** contains 222 term collections.
   Enable a collection and unhide its objects to examine its source surfaces.
3. **03 TA98 hierarchy — missing geometry is empty** retains all 7,474 entries
   and their source parent links. It is excluded from the view layer initially.
   Mapped entries link existing objects; they do not duplicate geometry.
   Empty entries record missing or non-surface representations.

Objects retain source IDs, licenses, provenance and anatomical review status.
Stomach objects additionally retain their source/target connection landmarks
and complete transformation in the `connection_registration` property.
Collections retain published TA98 codes, internal entity keys, coverage status,
and containing terms to inspect. Textures are packed. Coordinates are meters,
converted from the viewer's Y-up frame to Blender's Z-up frame.

## Remaining anatomy work

The audit distinguishes material structures from source-designated immaterial
entities and positional terms. Female genital structures are outside this male
model. Three source entries have contradictory material/immaterial flags and
are marked for source review.

There are 4,315 entries with containing geometry to inspect, 19 unresolved source
groups, 34 incomplete component groups, and eight entries without a located
corresponding surface or modeled ancestor in the supplied parent graph.
These categories identify review work; they do not prove geometry is absent.

The geometry report checks welded connected components and open/nonmanifold
edges. These can help locate loose fragments and defects. Neither a connected
component nor a spatial cut establishes an anatomical boundary. A named bone
landmark, cortical area, organ layer, vessel segment or muscle subdivision needs
an anatomically justified surface/volume or a reviewed face mask before receiving
its own mesh identity. An ancestor mesh is never labeled as the missing child.

The local reference body source has no established redistribution license in this
checkout. The review model stays in `.local-models` and is excluded from the
public production model menu and deployment. Original source assets are reused
by the viewer catalogue; the Blender project packs its textures.

## Rebuild and verification

```sh
npm run build:ta98-male-review
npm run test:ta98
npm run test:ta98-male-review
npm run check
/Applications/Blender.app/Contents/MacOS/Blender --background --factory-startup --python scripts/export-ta98-male-blender.py
/Applications/Blender.app/Contents/MacOS/Blender --background outputs/atlas-ta98-male-review.blend --python scripts/validate-ta98-blender.py
```

Blender's Metal initialization crashed in Codex's restricted execution sandbox.
The same installed Blender succeeded with normal host graphics access. The
export and reopen checks verify object counts, vertices/triangles, term
collections and parentage, hidden candidates, packed textures, and units.
Viewer checks verify default visibility, direct candidate inspection, and
shared URLs preserving both hidden and explicitly shown supplements.


Named hepatic/source batch: `named-hepatic-source-batch-v1/inventory.json` freezes ten additional qualified original-source identities (six hepatic branch families, tonsillar branches, diagonal branch, umbilical remnant and plantar arch), with exact source membership and authority parents. Fourteen cerebral regions and ambiguous ascending-colic/anorectal candidates remain deferred. This is source identity coverage, not anatomical completion; no geometry, source ancestry, approved heart pose, visibility, default installation or routing changes.

The named-source derivative `atlas-cardiac-source-identities-v1` has actual app/TA98/Systems validation for645 historical blood terms +4 cross-root terms, retaining32 source-heart terms and every spatial record/bank. Fifteen guarded source-family selectors now use logical Systems entries excluded from physical explosion ownership; all21 focused tests,7 plane tests,18 viewport tests in bounded batches, and exact-candidate plane/viewport checks pass. The changed TypeScript graph1044files passes; full cold-project checking exceeded1GiB and remains explicitly recorded in checkpoint.json. No routing or anatomical completion acceptance.

`next-source-coverage-v1/inventory.json` records two corrected source-parent homonyms: right pulmonary anterior/posterior segmental vessels had hepatic TA98 labels. The candidate assigns their actual pulmonary identities and retains hepatic anterior/posterior controls only as explicit right-hepatic V/VIII and VI/VII source-member unions; common stems remain unrepresented. Every installed pulmonary/hepatic binding is audited against the unchanged source forest. Geometry, pose, hierarchy, routing and anatomical acceptance are unchanged.

The source-parent correction candidate passed actual guarded TA98/Systems resolution for651 blood controls and32 source-heart terms. The647 unaffected prior memberships remain exact; the two hepatic homonyms are replaced by explicitly qualified right-hepatic source-member unions. All17 collection controls retain strict guards and selectionOnly behavior; forged pulmonary membership is rejected. Exact candidate plane1/1 and viewport4/4 pass with all banks/poses/source hierarchy literal. Not routed; anatomical completion remains unaccepted.

Correction of the previous source-union coverage claim: the named V/VIII and VI/VII daughter arteries do not prove the common anterior/posterior hepatic branch stems. A12.2.12.033/.034 are now parent-only placements, excluded from exact bindings, collection manifests and the exact numerator. Their old pulmonary homonyms remain correctly identified as A12.2.01.104/.107. The corrected base is next-source-coverage-v2; final apex profile atlas-cardiovascular-review-v4 has649 exact blood controls (645/1149 +4 cross-root),33 heart terms and15 exact source collections, with actual TA98/Systems proof. Earlier651 exact-union candidates are superseded and must not be installed. No geometry, pose, original source hierarchy or app routing changed in this correction.

`cerebral-course-regions-v1/inventory.json` inventories every ACA/PCA descendant before partitioning. Four planned pre/postcommunicating course regions require actual named communicating-mouth and source-ring witnesses, disjoint original face ownership and reciprocal section contours. All coordinates, packed attributes, fixed anatomy and original mouths are fixed; no wall layer, lumen or finite-tissue acceptance is inferred.

`cerebral-inflow-regions-v1/inventory.json` records every authoritative ICA, vertebral and MCA descendant, parent, current donor source and missing representation before any further partition. Ten possible course regions require actual carotid-canal/cavernous/dural, C6/C1/foramen-magnum or insular witnesses; clinical three-part axillary/subclavian regions are absent from TA98. All source attributes, heart/coronary poses, skeleton and adjoining branches remain fixed. Open source mouths and provisional boundaries remain explicit; no wall, lumen or microscopic representation is implied.

`ocular-course-regions-v1/inventory.json` records all authoritative central-retinal artery/vein descendants, source parts and parents before partition. Four extra/intraocular course regions require the original registered collection donor CN II ocular aperture, original retinal/scleral context and original source-edge contours; native eye coordinates are not identity evidence. Every original source face, coordinate and packed attribute remains fixed and covered once. No filled eye opening, new lumen/wall, microscopic feature, whole-body fit or anatomical acceptance is claimed.

The ocular source boundary contract now explicitly retains every source facet touching the original optic-disc aperture slab as an untyped transition; it is assigned neither extraocular nor intraocular exact identity. Actual triangle/aperture intersection loops must all lie inside the original CN II aperture. Typed regions contain only faces entirely on their own side of the complete nonplanar-rim slab. The first original edges may fan laterally inside the eye and remain coarse face-boundaries, not invented circular anatomical rims. Original whole-course selection includes all three exhaustive owners.

Ocular course candidate `ocular-course-regions-v1/candidate-v1` passes actual production guards, TA98 and Systems for 669 blood controls (665 historical male terms + 4 cross-root), preserving all 665 prior controls and 35 heart terms literally. Four source-backed extra/intraocular terms own 24,020 original triangles exactly once across 12 indexed parts; 324 ambiguous optic-disc transition facets remain untyped. All 12,190 attribute rows, world poses, normals/signs, UVs/materials, existing banks and RHA source branches replay exactly. Actual artery/vein disc-crossing loops (2/3 per side) lie inside the original donor CN II ocular aperture; actual triangle projections and side sections inspected. Original source facet granularity, 37.7 micrometre rim warp, 4–27 micrometre existing retina gap, complete region anatomy, native-eye fit and wall/lumen/patency remain unaccepted. No routing/default/app change. The other ten ICA/vertebral/MCA terms remain explicitly deferred; see their source-landmark checkpoint.

`mca-source-regions-v1/construction-obstacle.json` records the bounded original-source MCA follow-up. Faces241/243/244 meet at one exact-position edge;244/245 are a nonzero37.48-micrometre-wide sheet with distinct packed normal seams, not duplicate triangles. The diagnostic core0–243 has three6-point mouths and no nonmanifold edge, but the complete source sheet was neither deleted nor welded. Exact own-donor main/superior mouths agree within0.487micrometre; inferior mouths within1.148micrometres, both still unaccepted finite joins. The inferior terminal source extends44.51mm from original donor insular triangles, so its whole mesh cannot stand for M2. Original donor insula, actual incident-face and course projections inspected. Limen/circular-sulcus boundary and complete exterior-envelope repair remain unproved; no new term/candidate, source/bank/pose, RHA/retinal, app or route change.

`mca-source-regions-v1/inventory.json` now explicitly inventories every MCA authority descendant, parent and current source before source-envelope cleanup. Only the two original nonduplicate sheet faces244/245 per side are under repair review; the full original variant remains hidden for rollback. All three original six-point mouths, source attributes, source identities, adjoining family, heart, RHA/retinal geometry and poses remain fixed. No M1/M2 term, finite wall, cap or lumen/patency claim is added. Actual source facets and uncapped cross-sections are required for classifying the sheet.

`mca-source-regions-v1/inventory.json` now explicitly inventories every MCA authority descendant, parent and current source before source-envelope cleanup. Only the two original nonduplicate sheet faces244/245 per side are under repair review; the full original variant remains hidden for rollback. All three original six-point mouths, source attributes, source identities, adjoining family, heart, RHA/retinal geometry and poses remain fixed. No M1/M2 term, finite wall, cap or lumen/patency claim is added. Actual source facets and uncapped cross-sections are required for classifying the sheet.
