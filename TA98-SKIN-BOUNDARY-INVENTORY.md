# Whole skin boundary: TA98 inventory and repair scope

Existing whole skin boundary only. Do not infer dermal/epidermal layers, appendages, sensory endings or microscopic details from this coarse surface.

Remove only an exact, oppositely oriented coincident triangle pair forming a zero-volume fin if local topology confirms it. Do not change any coordinates, normals, UVs, body pose or anatomical openings.

| TA98 | Parent | Latin term | Classification | Male source | Female source |
|---|---|---|---|---|---|
| A15.2.07.037 | A16.0.00.014 | cilia | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.002 | A16.0.00.001 | cutis | whole skin | FJ2810 | missing |
| A16.0.00.003 | A16.0.00.002 | sulci cutis | boundary/line (not separate tissue) | missing | missing |
| A16.0.00.004 | A16.0.00.002 | cristae cutis | boundary/line (not separate tissue) | missing | missing |
| A16.0.00.005 | A16.0.00.002 | retinacula cutis | connective attachment | missing | missing |
| A16.0.00.006 | A16.0.00.002 | retinaculum caudale | connective attachment | missing | missing |
| A16.0.00.007 | A16.0.00.002 | toruli tactiles | microscopic surface/tissue feature | missing | missing |
| A16.0.00.008 | A16.0.00.002 | lineae distractiones | boundary/line (not separate tissue) | missing | missing |
| A16.0.00.009 | A16.0.00.002 | epidermis | wall layer | missing | missing |
| A16.0.00.010 | A16.0.00.002 | dermis | wall layer | missing | missing |
| A16.0.00.011 | A16.0.00.010 | stratum papillare | wall layer | missing | missing |
| A16.0.00.012 | A16.0.00.011 | papillae | microscopic surface/tissue feature | missing | missing |
| A16.0.00.013 | A16.0.00.010 | stratum reticulare | wall layer | missing | missing |
| A16.0.00.014 | A16.0.00.002 | pili | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.015 | A16.0.00.014 | lanugo | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.016 | A16.0.00.014 | capilli | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.017 | A16.0.00.014 | supercilia | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.018 | A16.0.00.014 | barba | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.019 | A16.0.00.014 | tragi | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.020 | A16.0.00.014 | vibrissae | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.021 | A16.0.00.014 | hirci | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.022 | A16.0.00.014 | pubes | material structure; anatomical extent needs review | REF:abdomen-pelvis:195:0 | missing |
| A16.0.00.023 | A16.0.00.010 | folliculus pili | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.024 | A16.0.00.014 | m. arrector pili | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.025 | A16.0.00.014 | flumina pilorum | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.026 | A16.0.00.014 | vortices pilorum | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.027 | A16.0.00.014 | cruces pilorum | material structure; anatomical extent needs review | missing | missing |
| A16.0.00.028 | A16.0.00.002 | glandulae cutis | microscopic surface/tissue feature | missing | missing |
| A16.0.00.029 | A16.0.00.028 | glandula sudorifera | microscopic surface/tissue feature | missing | missing |
| A16.0.00.030 | A16.0.00.028 | glandula sebacea | microscopic surface/tissue feature | missing | missing |
| A16.0.00.031 | A16.0.00.002 | terminationes nervorum | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.001 | A16.0.00.002 | unguis | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.002 | A16.0.01.001 | matrix unguis | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.003 | A16.0.01.001 | vallum unguis | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.004 | A16.0.01.001 | corpus unguis | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.005 | A16.0.01.001 | lunula | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.006 | A16.0.01.001 | margo occultus | boundary/line (not separate tissue) | missing | missing |
| A16.0.01.007 | A16.0.01.001 | margo lateralis | boundary/line (not separate tissue) | missing | missing |
| A16.0.01.008 | A16.0.01.001 | margo liber | boundary/line (not separate tissue) | missing | missing |
| A16.0.01.009 | A16.0.01.001 | perionyx | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.010 | A16.0.01.001 | eponychium | material structure; anatomical extent needs review | missing | missing |
| A16.0.01.011 | A16.0.01.001 | hyponychium | material structure; anatomical extent needs review | missing | missing |

A named source mapping does not establish a distinct layer, microscopic structure or accepted representation. Classifications are planning labels pending anatomical review. No additional TA98 coverage is accepted by this topology repair.

## Additional staged boundary repair plan

After removing the coincident zero-volume scalp fin, three pre-existing skin
self-intersection pairs remain at the anterior pelvic surface. All involve one
triangle. Before changing any coordinates, test only legal edge flips within
its existing two-triangle neighborhoods. Preserve the quad perimeter, shared
world coordinates, source ownership, vertex/UV storage, face count, component
count and orientable genus. Reject flipped or degenerate faces, open/nonmanifold
edges, changed vertex links or remaining self-intersections. Record volume change
and inspect the actual render and cross-sections before any installation. An
intersection pass alone does not establish accurate genital shape or skin layers.

## Staged male boundary evidence

The current male FJ2810 skin has an opposite coincident face pair at the scalp
(original face IDs 25031/25110), one nonmanifold edge, two invalid vertex links
and three additional self-intersection pairs near the anterior pelvic surface.
Removing the zero-volume face pair preserves all coordinates and volume and
leaves a closed orientable genus-one boundary. Three self-intersection pairs
remain, all involving face 19943 in the staged index order.

One legal local edge flip of faces 19943/19944 preserves the quad perimeter,
all vertices and component/genus counts while eliminating all three crossings.
The whole 44,742-triangle boundary then has zero open/nonmanifold edges,
duplicate/degenerate faces, inconsistent winding, invalid vertex links and
self-intersections. Volume changes −6.93e−9 m³ (about −0.00693 mL). A different
legal edge flip increased crossings and was rejected; a third reversed a face
normal and was rejected. No vertex was moved and no source model was overwritten.

Actual Blender full-body, local-surface and clipped cross-section before/after
renders were inspected. Variants are never superimposed. These checks prove a
staged topology repair, not all anatomical skin details or acceptance of its
42-term TA98 subtree. Main catalogue installation and source-bound regenerated
runtime verification remain outstanding. Local evidence: skin-zero-volume-fin-
repair.json, skin-crossing-edge-flips.json and skin-topology-repair-render.json
under reports/ta98-male/. Current source skin containment was never used as an
accepted solid before this repair.

## Local male installation and female audit

The reviewed male repair is now installed in the local master, lossless runtime
and ordinary bowel-layer atlas. Every stored position, normal, UV and other
vertex attribute is preserved; only the index array changes. Skin uses 44,742
triangles, two fewer than before. All three catalogues independently replay the
reviewed coordinates and triangles. The complete male dependency forest is
rebuilt; skin-free pelvic contexts are rebound after proving their geometry is
unchanged. The runtime source identity and gzip/attribute replay pass. The
runtime builder reapplies the optional source-bound repair when local evidence
exists; stale coordinates are skipped rather than silently modified. No public
deployment or ignore-rule change occurs.

The installed main atlas was rendered and skin selection restored through its
native view API. Both immediate and settled camera values remain exact, with no
browser errors. This installs a topology repair, not additional dermal layers or
complete anatomical skin acceptance. Evidence: skin-topology-installed.json,
skin-topology-install-validation.json and skin-repair-context-refresh.json under
reports/ta98-male/.

The current female audit finds 56 open edges / 56 triangle intersection pairs in
LOCAL:female:skin:180 and 66 open edges / 160 pairs in
LOCAL:female:skin:1005, with no nonmanifold edges in either. Their openings must
be classified against actual adjoining anatomy before any closure. They remain
unchanged pending rendered regional review; no blind caps or cross-sex coordinate
transfer is authorized by this audit. Evidence: reports/ta98-female/skin-boundary-
audit.json. The full male/female TA98 goal remains incomplete.

## Female assembly classification and lossless indexing stage

The two female skin meshes are a shared body/head surface, not two whole-body
skins. Their 26 shared boundary coordinates match exactly. Welding positions
for diagnosis yields one combined surface with four open loops and 70 boundary
edges; the matching body/head seam disappears without adding a cap. Loop
neighborhoods are recorded against unchanged source organs. Nearest source
vertices locate the urethra (0.13 mm from the 8-edge loop), vagina (1.75 mm from
the 10-edge loop), mouth (2.97 mm from the 40-edge loop), and the source labeled Large intestine
anus (1.21 mm from the 12-edge loop). The rectum is farther away (22.0 mm). These measurements do not certify continuity or opening
identity. Anatomical surfaces are not reconstructed from the distances.

Six bilateral crossing clusters are on the body mesh (56 pairs), including
shoulder/axillary, foot and hand locations. The head mesh has 160 pairs near the
mouth opening. Actual regional renders are reviewed; lip apposition and genuine
fold defects need separate assessment. All openings remain uncapped and no
female source coordinates or anatomy were modified.

A separate storage-only stage indexes identical *complete* vertex records,
including position, normal, UV and every other present vertex attribute. Body
records reduce 36,294 → 7,134 and head records 12,984 → 2,407. Geometry payload
reduces 1,478,340 → 445,180 bytes (~70%). Every attribute on every triangle
replays bit-exact; texture seams and normals are preserved. Faces, openings,
source identity and anatomical crossing counts are unchanged. This is a runtime
resource improvement, not a skin-shape repair or accepted additional TA98 term.
It remains staged pending installation and ordinary-viewer verification.
Evidence: skin-assembly-openings.json, skin-defect-locations.json,
skin-boundary-render.json and skin-attribute-dedup.json under reports/ta98-female/.

## Female lossless indexing installed locally

The full-record indexing is installed in the female master, lossless runtime
and ordinary local atlas. Body/head skin retains all 16,426 triangles and every
triangle's position, normal and UV values, with 9,541 stored records instead of
49,278. Its uncompressed payload is 445,180 bytes instead of 1,478,340 (~70%
less). Unreferenced chunks are removed from catalogue accounting while old files
are retained. No textures, images, primitive identity or shape changes occur.
All three catalogue triangle-attribute hashes replay independently. Body/head
open-edge counts and self-intersection counts are unchanged. This is indexing,
not simplification, smoothing, repair of the 216 pairs or acceptance of skin
layers. The female dependency forest is rebuilt and one geometrically unchanged
context is rebound to the new catalogue hash.

The ordinary viewer renders both skin pieces together with the existing texture
finish, identifiers and hierarchy. The optional local runtime build hook preserves
this indexing when source triangle attributes match staged evidence and skips
stale data without crashing the build. The male repaired skin and all other
anatomy remain geometrically unchanged. No public deployment or ignore-rule
change occurs. Evidence: skin-indexing-installed.json, skin-indexing-validation.json
and skin-indexing-context-refresh.json under reports/ta98-female/.

## Female shoulder repair plan

Current storage-indexed skin retains all 56 body and 160 head crossing pairs.
All are noncoplanar under the existing checker; this still does not distinguish
intended lip apposition from a fold defect. Shoulder work is restricted to the
existing body surface, not the mouth or four open anatomical interfaces. The
whole body/head seam and adjoining organs remain fixed.

Test only existing shoulder quad diagonals with shared *full* vertex records
across the old edge, preserving normal/UV seams. Preserve all coordinates,
normal/UV arrays, face counts, boundary-loop coordinate sets and connected
components. Reject new crossing pairs, face-normal reversals, degenerate or
nonmanifold triangles and changes in openings. At most six improving local flips
are staged; actual textured and section renders must precede installation.
This does not assign an additional skin layer or an accepted microscopic feature.

## Female shoulders: local repair installed

All 216 reported female skin pairs are noncoplanar under the current checker.
That geometric classification alone does not determine intended lip apposition.
Six legal existing body-skin quad diagonals eliminate the 22 shoulder/axillary
crossing pairs, reducing body-skin pairs 56 → 34. No new crossing pairs are
introduced. Shared full records protect normal/UV seams; every vertex, attribute
array, triangle count, body/head seam and opening loop stays unchanged. Positive
face orientation and original connected-component topology are preserved.
Remaining body pairs are on the hands/feet; the head/lip 160 pairs are unaffected.
Thus this is a shoulder-region repair, not a completed whole-skin repair.

Actual source/staged front, back and clipped-section renders were inspected with
original normals and mutually exclusive geometry variants. The local master,
runtime and ordinary atlas now contain the repaired indices. Every installed
candidate replays, full crossing checks confirm 34 body pairs and zero shoulder
pairs, and the runtime byte/gzip checks pass. The female dependency forest and
skin-body-free pelvic reference are refreshed. Geometry remains at the reduced
storage cost; no vertex movement, new UV/image or public deployment occurs.
The optional runtime build hook retains the correction when coordinates match
and skips stale source evidence rather than failing a resized model's build.

Future local searches use changed-triangle versus full-surface checks and a final
whole-surface replay instead of repeatedly recomputing all unchanged pairs.
Evidence: shoulder-skin-retriangulation.json, shoulder-skin-render.json,
shoulder-skin-installed.json and shoulder-skin-install-validation.json under
reports/ta98-female/. Anatomical skin layers and the whole-body TA98 goal remain
incomplete.
