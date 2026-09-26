# Section geometry review options

The audit counts **source primitives**, not distinct anatomical defects. In the
Male reference GLBs, 2,679 of 5,094 primitives have boundary edges and 204 have
nonmanifold edges or winding conflicts. Of the open primitives, 1,592 belong to
concepts made of multiple material pieces; the pieces may close when combined.
The 1,087 single-piece open meshes also include intentional sheets, cut skin
regions, vessel ends, and hollow organ surfaces. A boundary flag alone is not
evidence that tissue should fill its interior.

At the reviewed torso axial plane (70%), 99 flagged open primitives intersect
the plane. These include the derived skin, 18 named skin-region pieces, 11
fascia pieces, 33 vascular pieces, 13 muscle pieces, 4 skeletal pieces, and 2
digestive pieces. The per-part figures are in `local-reference.json`.

## Choices

| Choice | What changes | Benefit | Anatomical risk |
| --- | --- | --- | --- |
| A. Review by tissue and concept (recommended) | Combine material pieces of one solid concept; retain narrow wall cuts for skin, fascia, vessels, and digestive lumen; repair individual solid concepts only after a screenshot and contour check. | Fixes false openings without inventing tissue. | Some source defects remain visible until reviewed. |
| B. Fill every closed-looking contour automatically | Bridge boundary loops and cap every mesh the algorithm judges solid. | Fewer immediate gaps. | Can seal intestines, skin, vessels, and spaces that must remain open; overlapping contours can create rotation artifacts. |
| C. Preserve source geometry exactly | Draw only the existing cut edge when a source primitive is open. | Avoids adding false volume. | Many muscle and bone sections remain visibly hollow. |

The current implementation follows **A** for reference bones, muscles, and lung
lobes. For example, the two rectus abdominis material pieces at y=1.05257 each
look open alone (boundary ratios 0.057 and 0.118), but the combined concept is
closed and yields a filled section. The right superior lung lobe at y=1.355
changes from three open cut edges to one filled contour when its intersecting
material pieces are combined. The duodenum remains a wall-only cut; its lumen
is not filled. `node scripts/inspect-cap-group.mjs <concept-id> <height>` prints
these contours and the cap choice for any reference concept.

## First individual review candidates

1. **L2 and L3 vertebrae:** their cartilage primitives at the axial 70% plane
   have boundary ratios 0.250 and 0.243. Verify against the neighboring bone
   primitives before adding a cartilage cut face.
2. **Serratus posterior inferior, left and right:** paired muscle/tendon
   primitives cross that plane and include a few nonmanifold edges. Review the
   combined contour from superior and turned views before any geometry repair.
3. **Longus colli, left and right:** each has over 3,400 nonmanifold edges in
   the source audit. Inspect where the problem intersects a section before
   deciding whether to repair or keep the source surface.
4. **Cavernous sinus, left and right:** roughly 4,480 nonmanifold edges each.
   Treat as a venous space; do not convert it into a solid block by default.

For every manual repair, retain source triangles and material provenance in the
import record, store the repair as a separate reproducible transform, and verify
axial, coronal, sagittal, and rotated views with click/hover selection.
