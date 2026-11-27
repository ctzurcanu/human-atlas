# Male skin, eyes and shared curved brush-hair appearance

Native ZA skin regions establish the male surface frame. The alternative FJ2810 skin remains hidden. Native coordinates, normals, UVs, anatomical IDs and camera behavior stay unchanged. Restore visible ocular material distinctions on the existing eye parts; do not add another eye pair. Scalp hair is a separate source-scoped parent representation, with 1 mm length, 0.1 mm diameter and approximately 0.6 mm neighboring-root spacing. Sample original scalp triangles against a continuous smooth boundary, rather than partitioning at coarse triangle corners. Native skin remains visible between shafts; use one shared UV-authored prototype, no hair image files. Crown, temporal, occipital and nape regions are included; auricles and face excluded.

| TA98 | Parent | Term | Current source | Representation |
|---|---|---|---|---|
| A16.0.00.014 | A16.0.00.002 | pili | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A15.2.07.037 | A16.0.00.014 | cilia | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.015 | A16.0.00.014 | lanugo | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.016 | A16.0.00.014 | capilli | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.017 | A16.0.00.014 | supercilia | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.018 | A16.0.00.014 | barba | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.019 | A16.0.00.014 | tragi | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.020 | A16.0.00.014 | vibrissae | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.021 | A16.0.00.014 | hirci | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.022 | A16.0.00.014 | pubes | REF:abdomen-pelvis:195:0 | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.024 | A16.0.00.014 | m. arrector pili | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.025 | A16.0.00.014 | flumina pilorum | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.026 | A16.0.00.014 | vortices pilorum | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |
| A16.0.00.027 | A16.0.00.014 | cruces pilorum | none | scalp shaft appearance planned only; microscopic layers/follicles remain missing |

The hair family is not a completed follicle/integumentary assembly. Its related scalp skin, auricles, brows, eyelids and native eyes remain separate anatomy. Body skin color uses shared native materials with a darker warm complexion. Ocular exteriors/iris/cornea retain distinct appearance; pigmentation does not establish anatomical completion. Check actual front/back/side renders, source-coordinate integrity, root spacing/dimensions and camera/explosion invariants. All model work stays local.


Current local implementation: 154,912 measured shafts on 108 exact native ZA head/neck skin sources; mean root spacing 0.599280 mm. Requested length 1 mm and diameter 0.1 mm are shared with female. The prior ZA head-hair mesh is preserved suppressed in the master and excluded from active runtime/review catalogues and concepts. Native eye geometry stays unchanged; off-white ocular exterior, brown iris and low-opacity cornea/optical media replace the uniform ochre appearance. Existing male skin-region color changes from [191,151,126] to [176,126,95], with one reused pore tile. Native geometry hashes and runtime preservation pass, as do the shared fiber/visibility and TypeScript checks. No public deployment.

Live close-up verification is pending: the initial browser showed the darker male body but still loaded the old hair. The latter is now retired. Browser/computer-use previews became unavailable before the final close-up, so do not treat these data/geometry checks as completed appearance acceptance.


Preview access recovered. Actual male full/close frontal and lateral views now show darker skin, visible native eyes and the replacement short brush hair. Old head hair is absent. Scalp hair selection preserves the exact DOM camera pose, including projection offset. No browser errors and zero pending geometry chunks were observed. Final images: outputs/male-curved-brush-front.png and male-curved-brush-side.png. Current male head/neck review still submits 12,867,636 triangles and holds 459,419,630 geometry-array bytes, including the existing main atlas; this is not an efficiency acceptance or FPS benchmark. The hair field shares one shaft prototype and adds approximately 4 MB of root attributes, but full model delivery/render efficiency remains unfinished.


## Brow and pubic-hair appearance correction

| TA98 | Parent | Term | Current source | Action |
|---|---|---|---|---|
| A16.0.00.017 | A16.0.00.014 | supercilia | ZA:Eyebrow.l, ZA:Eyebrow.r, ZA:Hairs of eyebrow.l, ZA:Hairs of eyebrow.r | macroscopic hair appearance; follicles/coats unrepresented |
| A16.0.00.022 | A16.0.00.014 | pubes | REF:abdomen-pelvis:195:0 | macroscopic hair appearance; follicles/coats unrepresented |

Per user instruction, exclude pubic-hair sources from the active male catalogue, preserving originals in the master/assets. Make only the two eyebrow-hair sources lighter brown and narrow their vertical extent to 60% about the existing center. Preserve X/Z, every triangle/UV and native skin/eye geometry. Transport source normals with the inverse transpose of the same positive affine map. This is a cosmetic source adjustment, not complete anatomical repair or a hair-follicle representation. Check actual frontal appearance, unchanged topology and source attribute preservation.

Actual front view confirms lighter brown, narrower brows. Active runtime/review pubic hair is absent; master original remains suppressed. Real installed brow indices and UVs match original source hashes; height replay matches the 60% target. Six appearance/visibility tests and runtime preservation pass. Screenshot: outputs/male-slim-brown-brows.png. No skin, eye or skeletal coordinates were moved.

## Paired eyes, bare ears and complexion correction, 2026-10-07

User-directed appearance correction excludes all 44 named native auricular skin sources from the male scalp sampling domain. The regenerated field has 152,995 shafts, mean nearest-root spacing 0.599067 mm, with unchanged 1 mm length and 0.1 mm diameter. Native ear and scalp geometry is retained; no head hair is sampled from an auricular source. Both lateral renders confirm bare ears. Original assets remain preserved.

Skin material changes slightly from [176,126,95] to [165,116,88]. A soft procedural vermilion mask on existing native oral skin adds muted rosy-brown lips, preserving surrounding face pigmentation and all surface coordinates. Eye exteriors use neutral white without the illustrative grain map; the posterior surface has a separate dark appearance. A mirrored native-coordinate pupil pigment field on existing eye surfaces makes both centers black, including the original iris/optical surfaces that previously showed a grey reflective center. This is optical display appearance, not a new solid pupil or acceptance of detailed ocular anatomy. No new meshes, image files or texture buffers are added for lips or pupils.

Actual frontal and both lateral views were inspected in a separate diagnostic browser tab; the user's current camera was not navigated. Screenshot evidence: outputs/male-face-corrected-front.png, male-ear-hair-corrected-side.png, male-ear-hair-corrected-opposite-side.png. TypeScript, runtime geometry/attribute preservation, six shared appearance/visibility tests, camera intent and all five explosion-plane cases pass. Whole-body anatomy and full delivery efficiency remain incomplete.

Correction after user review: the preceding eye-white claim missed the actual native `ZA:Sclera.l` and `.r` parts. Both still used the original ochre [225,188,104] material despite the corrected generic anterior surfaces. The installer now explicitly includes Sclera and assigns both white [255,255,255] materials without the illustrative grain map. An existing appearance regression checks both actual sclera IDs. The frontal render now shows white sclera on both sides with black pupil centers: outputs/male-white-sclera-corrected.png. Earlier screenshots are retained as history, not proof of this correction.
