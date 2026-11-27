# TA98 coronary arteries: coverage and assembly plan

Complete ascending-aorta/coronary descendant inventory before any geometry reconstruction. Source terminology correspondence is distinct from certified walls, lumens and perfusion territories. Independent reference donors are recorded explicitly; they are not native branches.

| TA98 | Name | Parent | Male source records | Female source records |
|---|---|---|---|---|
| A12.2.03.001 | ascending aorta | A12.2.02.001 | ZA:Ascending aorta (exact, native) | HRA:VH_F_ascending_aorta (exact, independent reference, hidden) |
| A12.2.03.002 | aortic sinus | A12.2.03.001 | missing | missing |
| A12.2.03.003 | supravalvular ridge | A12.2.03.001 | missing | missing |
| A12.2.03.004 | aortic bulb | A12.2.03.001 | missing | missing |
| A12.2.03.101 | right coronary artery | A12.2.03.001 | ZA:Right coronary artery (exact, native) | missing |
| A12.2.03.102 | atrioventricular branches | A12.2.03.101 | missing | missing |
| A12.2.03.103 | conus branch | A12.2.03.101 | FJ2670 (exact, independent reference, hidden); FJ2676 (exact, independent reference, hidden) | missing |
| A12.2.03.104 | sinu-atrial nodal branch | A12.2.03.101 | missing | missing |
| A12.2.03.105 | atrial branches | A12.2.03.101 | missing | missing |
| A12.2.03.106 | right marginal branch | A12.2.03.101 | FJ2667 (exact, independent reference, hidden); FJ2668 (exact, independent reference, hidden); FJ2672 (exact, independent reference, hidden); FJ2674 (exact, independent reference, hidden); FJ2675 (exact, independent reference, hidden) | HRA:VH_F_right_marginal_artery (exact, independent reference, hidden) |
| A12.2.03.107 | intermediate atrial branch | A12.2.03.101 | missing | missing |
| A12.2.03.108 | posterior interventricular branch | A12.2.03.101 | FJ2692 (exact, independent reference, hidden); FJ2693 (exact, independent reference, hidden); FJ2694 (exact, independent reference, hidden); FJ2695 (exact, independent reference, hidden); FJ2696 (exact, independent reference, hidden); FJ2697 (exact, independent reference, hidden); FJ2698 (exact, independent reference, hidden); FJ2699 (exact, independent reference, hidden); FJ2700 (exact, independent reference, hidden) | missing |
| A12.2.03.109 | interventricular septal branches | A12.2.03.108 | missing | missing |
| A12.2.03.110 | atrioventricular nodal branch | A12.2.03.101 | missing | missing |
| A12.2.03.111 | right posterolateral branch | A12.2.03.101 | ZA:Right inferolateral branch of right coronary artery (exact source identity, native; wall/lumen/packing unverified) | missing |
| A12.2.03.201 | left coronary artery | A12.2.03.001 | ZA:Left coronary artery (exact, native) | missing |
| A12.2.03.202 | anterior interventricular branch | A12.2.03.201 | missing | missing |
| A12.2.03.203 | conus branch | A12.2.03.202 | FJ2643 (exact, independent reference, hidden); FJ2644 (exact, independent reference, hidden) | missing |
| A12.2.03.204 | lateral branch | A12.2.03.202 | missing | missing |
| A12.2.03.205 | interventricular septal branches | A12.2.03.202 | ZA:Septal branches of anterior interventricular artery (exact, native) | missing |
| A12.2.03.206 | circumflex branch | A12.2.03.201 | FJ2649 (exact, independent reference, hidden); FJ2650 (exact, independent reference, hidden); FJ2651 (exact, independent reference, hidden); FJ2652 (exact, independent reference, hidden); FJ2653 (exact, independent reference, hidden); FJ2654 (exact, independent reference, hidden) | HRA:VH_F_left_circumflex_artery (exact, independent reference, hidden) |
| A12.2.03.207 | atrial anastomotic branch | A12.2.03.206 | missing | missing |
| A12.2.03.208 | atrioventricular branches | A12.2.03.206 | missing | missing |
| A12.2.03.209 | left marginal artery | A12.2.03.206 | missing | HRA:VH_F_left_marginal_branch (exact, independent reference, hidden) |
| A12.2.03.210 | intermediate atrial branch | A12.2.03.206 | missing | missing |
| A12.2.03.211 | posterior left ventricular branch | A12.2.03.206 | missing | missing |
| A12.2.03.212 | sinu-atrial nodal branch | A12.2.03.206 | missing | missing |
| A12.2.03.213 | atrioventricular nodal branch | A12.2.03.206 | missing | missing |
| A12.2.03.214 | atrial branches | A12.2.03.206 | missing | missing |

## Assembly and interfaces

### Complete left-source seam consistency trial

Three explicitly named adjoining source interfaces have six mutually nearest boundary coordinates per side within 0.5 micrometres under the empirical body-frame scale: left main/anterior interventricular, left main/circumflex, and circumflex/left marginal. Measured maximum separations are 0.164690, 0.124727, and 0.100758 micrometres. Source endpoint normal discrepancies reach 7.810, 14.531, and 18.377 degrees respectively. This is source mesh correspondence evidence; no native placement, complete wall/lumen or perfusion claim follows.

Use the uncrossed triangular-gap candidate as a fixed fixture. Before accepting a seam change, require each six-point match to form a complete boundary cycle, have opposite oriented boundary edges, and be one-to-one. All original faces, UVs and material bindings remain; do not create bridge faces or overlapping caps. Snap the paired points to the same float32 mean and use a common normalized mean endpoint normal only if each side has no conflicting split-normal crease at that point. Preserve all other coordinates and normals exactly. Record every original packed vertex involved and the resulting displacement/normal angle. Reject partial boundary matches, wrong edge winding, changed isolated-part topology or increased crossing counts. Verify the joined four-object source surface and actual rendered finish; the source proximal aortic-root interface, other branch joins and native cardiac packing remain provisional. Keep earlier candidates intact.

### Local source uncrossing trial, before native fitting

The repair fixture is the entire 15-object rebuilt female arterial assembly with eight fixed source cardiac/root guides. Only two original objects are eligible for local vertex movement: atrial branches of circumflex (source node `342cdddb39bf1b1e375bfa58dc57e6dd`, faces 873/874/875) and ventricular branches of right coronary (node `7f2879963bad86ac8fd1274cc00c8115`, faces 48/98/204/205/207/502). Both are partial named arterial exterior networks, not established walls or lumens. Every other object is fixed. The original 29-term inventory above remains; no new partition or identifier is inferred from a repair patch.

All original triangle lists, UVs, material bindings, open boundary edges and cross-object exact shared coordinates are fixed. The atrial defect involves six distinct coordinates, none on its six-edge boundary. The ventricular defects involve eight coordinates, three on its fifteen-edge boundary; those three must remain fixed. Coordinate duplicates move together. Editable neighborhoods are restricted to those defect vertices and one adjoining edge ring, excluding protected coordinates. Trial displacement is bounded to 0.2 mm under the empirical body-frame scale; this is a numerical repair limit, not an anatomical tolerance or proven source calibration. Changes remain in a separate local candidate archive.

Candidate movement is conservative local averaging, with original connections and branch membership retained; no face deletion, end capping, shrink wrapping to incomplete cardiac guides or vessel merging is allowed. Transport source normals through local old/new face rotation, preserving their original differences, and leave UVs exact. Reject the trial if any protected coordinate/normal changes, component/boundary/link topology changes, degenerate or duplicate faces appear, face orientation reverses, or whole-assembly self/cross-object intersection review worsens. Render the actual candidate and the same local sections before judging improvement. Source ostia, bifurcations, tissue scopes and native myocardial/venous/pulmonary interfaces remain unresolved; this trial cannot authorize installation or anatomical acceptance.

Keep each native coronary root, left/right main arterial exterior, qualified anterior-interventricular and circumflex branches, septal group and right-coronary source branches in their original cardiac frame. A branch is adjoining vascular tissue, not a solid contained within its parent artery. Aortic sinuses/ridges are cavity/surface landmarks, not separate filler solids. Arterial walls, microscopic layers, individual septal territories and coronary anastomoses are not established by the source outer tubes.

The native female catalogue currently has no independently named coronary source meshes under these terms. Available cardiac donors are independent reference objects: do not move or substitute them without complete registration/packing evidence and licensing review. Both native Marginal artery records are colonic A12.2.12.068 sources and must stay outside coronary assemblies. The male inferolateral source branch has exact source-scoped A12.2.03.111 identity: [FIPAT TA2 Part 4 entry 4141](https://fipat.library.dal.ca/wp-content/uploads/2021/08/FIPAT-TA2-Part-4.pdf) and [TAH3784/FMA76994](https://ifaa.unifr.ch/Public/TNAEntryPage/auto/unit/LAFR/TAH3784%20Unit%20EN.htm) explicitly equate right inferolateral and right posterolateral branch names. Actual fixed-frame coronary surfaces and unchanged earlier packed-source vertices/faces were reviewed. Minimum native RCA vertex gap is 0.160296 mm with zero exact shared coordinates; joined walls, lumen continuity, dominance and perfusion remain unverified. See reports/ta98-right-coronary-source-identity-20261024.json.

Before any remodeling, define aortic coronary ostia, left-main bifurcation, septal takeoffs, epicardial sulci, vessel–myocardial/fat interfaces and venous crossings. Fix and review whole surfaces against the atrial/ventricular, aortic/pulmonary and venous neighbors before any partition. Coronary dominance and posterior-interventricular origin require qualified source evidence, not a universal right-origin assumption. Do not duplicate overlapping faces or cap vascular openings blindly.

## Source evidence

Actual packed float32 male arterial sources were decoded without edits. Exact shared-coordinate counts are recorded in local `reports/ta98-completeness/native-coronary-source-review.json`; they constrain source seams only, not biological lumen continuity. Complete topology, self-crossing, cross-section, actual-render and surrounding-tissue review remain required before acceptance. Source normals/UVs/materials and original geometry resources are preserved. Named branch paths follow [OpenStax coronary circulation](https://openstax.org/books/anatomy-and-physiology-2e/pages/19-1-heart-anatomy); no texture or diagram is copied.

## Rebuilt female source staging

The rebuilt Female_2024.glb source now has 15 coronary mesh candidates staged locally with eight fixed cardiac/root guides. The 23-object compressed archive is 567,335 bytes; coronary source geometry totals 14,708 original triangles. The GLB remains untouched. Independent direct accessor replay verifies 77 positions, normals, UV and index arrays; absent UVs remain absent. All eleven morph weights are zero in the reviewed default shape. Nonzero weights or implicit morph ingestion remain rejected. Original converted triangle indices, including zero-area source faces, are preserved for source-defect review.

| Source candidate | CID | Original node | Triangles | Open edges | Components |
|---|---|---|---|---|---|
| Right Marginal Branch of Right Coronary Artery | 17003 | ecdf5a2ac9299826043e956a42809e40 | 296 | 6 | 1 |
| Sinuatrial Nodal Branch of Right Coronary Artery | 20489 | 3df1a890cdbc842d5aa0a59e324c5ad9 | 4138 | 6 | 1 |
| Conal Branch of Right Coronary Artery | 23484 | d0e775e66e10d045fa054678df0616b7 | 458 | 6 | 1 |
| Inferior Left Ventricular Branch of Circumflex Artery of Heart | 18678 | f481dcf5612fa41cf31ba25fd7f3ffbf | 894 | 6 | 1 |
| Left Coronary Artery | 19880 | 43d530c0140a7f921b10af8401f29d75 | 240 | 16 | 1 |
| Anterior Interventricular Artery | 17349 | e890cc6cab0c386d3751375a6dd52de9 | 728 | 12 | 1 |
| Atrial Branches of Circumflex Artery of Heart | 19651 | 342cdddb39bf1b1e375bfa58dc57e6dd | 2534 | 6 | 1 |
| Diagonal Branch of Anterior Interventricular Artery | 18983 | 3380d1d3edba38229665379d952e8cfe | 356 | 6 | 1 |
| Ventricular Branches of Right Coronary Artery | 17065 | 7f2879963bad86ac8fd1274cc00c8115 | 1419 | 15 | 2 |
| Atrioventricular Nodal Branch of Right Coronary Artery | 18330 | 2d6b339cc210ac5c7f212b83ace2c019 | 186 | 6 | 1 |
| Right Inferolateral Branch of Right Coronary Artery | 22555 | fbd72a896c739b10264f78d80baebd35 | 498 | 6 | 1 |
| Left Marginal Branch of Circumflex Artery of Heart | 23934 | 93b33e9ac8521114633f1aae29c31fd3 | 418 | 6 | 1 |
| Right Coronary Artery | 16666 | 1d4df469e7c2c50c91595bca895eb60d | 1356 | 62 | 1 |
| Inferior Interventricular Artery | 24089 | f142a83496dd21c33ef47130f5a427c9 | 1003 | 9 | 1 |
| Circumflex Artery of Heart | 20694 | 05c057f070fc806e5f8b1ac5c31010e1 | 184 | 24 | 1 |

Source names and original material records are retained as local review evidence. Native female coronary coverage remains missing until registration and anatomical acceptance. Default source cardiac guides are cut/open surfaces; native atrial sources include adjoining venous stubs, so similarly named objects do not establish matching tissue scopes. The source and native heart projections use independent axes/units. No bounding-box fit, coordinate deformation, image copying, caps, resampling or native installation occurs. No new resources are loaded by the application.

Six named source-interface checks find no exact shared coordinates at aortic-root/right-left coronary, left-main/anterior/circumflex, or right-coronary/marginal joins. Circumflex/left-marginal shares four coordinates; that alone does not prove continuous tissue or lumen. Open coronary ends cannot be closed blindly; myocardial guide openings must not be treated as solid cutters. Establish coordinate units and stable matched aortic-root/annular and body landmarks, investigate source-scope differences, then review full self-crossings, actual sections and native myocardial/venous/pulmonary packing before any installation. Preserve the entire 29-term inventory and distinguish openings, regions, walls and individual branches during reconstruction.


## Female coronary body-frame trial and source crossing review

The default rebuilt source and native female body have been compared using the actual surfaces of 24 named fixed skeletal guides. Nine exterior bounds centers initialize a proper uniform similarity; fifteen are held out. These centers are scope-dependent numerical proxies, not selected homologous anatomical points. The original models and all source arrays remain untouched. The fitted scale 0.025255582001075232 is empirical, not a physical unit declaration.

Held-out proxy errors span 2.561–12.984 mm, median 5.851 mm. Sternum training error is 20.433 mm; fifth lumbar error is 16.417 mm. Actual front/side overlay reveals ventricular mismatch: source wall union lower extent lies 34.482 mm above the native whole-ventricular lower extent. Different scopes prevent interpreting that as a matched apex displacement. The body fit is rejected for installation and retained only as initialization for constrained cardiac fitting.

| Fixed guide | Source node | Native concept | Role | Proxy residual (mm) |
| --- | --- | --- | --- | --- |
| Sternum | 6534601f242f8ee467bdfd9cc8cbac78 | LOCAL:female:sternum | Training | 20.433 |
| Sacrum | 2fac851ee8a4f32930560d4f14d97645 | LOCAL:female:sacrum | Training | 3.957 |
| Third Cervical Vertebra | e59d7016f8c34b7a68d817b6c46d0340 | LOCAL:female:cervical_vertebra_c3 | Training | 8.844 |
| Fourth Cervical Vertebra | db875cc2e1dd78053bebf73d3e608e1b | LOCAL:female:cervical_vertebra_c4 | Held out | 7.099 |
| Fifth Cervical Vertebra | 2478033f8455c9b50f0248321aa2dd84 | LOCAL:female:cervical_vertebra_c5 | Held out | 5.714 |
| Sixth Cervical Vertebra | 4396ea5c767b36a69548ca488fa22a79 | LOCAL:female:cervical_vertebra_c6 | Held out | 3.791 |
| Seventh Cervical Vertebra | f79176a8183732ed3546977020b4b01b | LOCAL:female:cervical_vertebra_c7 | Training | 0.683 |
| First Thoracic Vertebra | 9bff40723c5ca47425fbe870f21b1ae2 | LOCAL:female:thoracic_vertebra_t1 | Training | 3.184 |
| Second Thoracic Vertebra | bac44e514034033b1d7348a2025f8467 | LOCAL:female:thoracic_vertebra_t2 | Held out | 2.561 |
| Third Thoracic Vertebra | 740ad1b9f688c1bd7f20b290d62aea70 | LOCAL:female:thoracic_vertebra_t3 | Held out | 2.619 |
| Fourth Thoracic Vertebra | e5201367eece9249c7a8851f6462bd88 | LOCAL:female:thoracic_vertebra_t4 | Held out | 7.269 |
| Fifth Thoracic Vertebra | 3f48c5da73131db34a7813b8a0bab765 | LOCAL:female:thoracic_vertebra_t5 | Training | 9.223 |
| Sixth Thoracic Vertebra | 0d9fe0300d564dac3d64a0eb9d05a729 | LOCAL:female:thoracic_vertebra_t6 | Held out | 9.199 |
| Seventh Thoracic Vertebra | a2781562afe3f3bb08489f54f5913414 | LOCAL:female:thoracic_vertebra_t7 | Held out | 5.851 |
| Eighth Thoracic Vertebra | 78c214e3a607d6c84eedc247c2983220 | LOCAL:female:thoracic_vertebra_t8 | Held out | 4.879 |
| Ninth Thoracic Vertebra | 55742e1d63f54f34a6427bf2958f65c6 | LOCAL:female:thoracic_vertebra_t9 | Training | 4.513 |
| Tenth Thoracic Vertebra | d650dad34d9f76922e2d911e90679273 | LOCAL:female:thoracic_vertebra_t10 | Held out | 5.892 |
| Eleventh Thoracic Vertebra | 966da3c8c51ee10305774557841a6dd5 | LOCAL:female:thoracic_vertebra_t11 | Held out | 11.578 |
| Twelfth Thoracic Vertebra | dcdbaa112504ff0ccd6c2c71fa1c4da8 | LOCAL:female:thoracic_vertebra_t12 | Held out | 12.984 |
| First Lumbar Vertebra | 8a7882a8bf68c5c138e22c3c395a388f | LOCAL:female:lumbar_vertebra_l1 | Training | 9.163 |
| Second Lumbar Vertebra | 49769a731fd6c2ce97464f635cab1742 | LOCAL:female:lumbar_vertebra_l2 | Held out | 3.555 |
| Third Lumbar Vertebra | e6aec58df72854e3645f154b33161c5c | LOCAL:female:lumbar_vertebra_l3 | Held out | 3.125 |
| Fourth Lumbar Vertebra | ef3f705b32309dad10ff5b6b71cdc24b | LOCAL:female:lumbar_vertebra_l4 | Held out | 8.649 |
| Fifth Lumbar Vertebra | 234b746630222377996bb5678535ae0f | LOCAL:female:lumbar_vertebra_l5 | Training | 16.417 |

All 15 coronary source meshes underwent whole-surface self-intersection checks and all 105 distinct arterial mesh-pair checks in the diagnostic frame. Float64 triangle interior/coplanar checks find seven self-intersection pairs: two in atrial circumflex branches (faces 873/874 and 873/875), and five in right-coronary ventricular branches (48/98, 98/205, 204/207, 204/502, 205/207). All seven actual local front/side projections and transverse sections were rendered and inspected. Zero detected intersections between distinct arteries does not prove junction continuity, complete packing, tissue-volume separation or lumen anatomy. Original faces, normals, UVs and material records are retained. No cleanup, cap, deformation or native installation is performed.

Independent original-accessor validation still passes for all 23 objects/77 arrays and confirms native atlas hash unchanged. Existing narrow-phase validator passes 864 winding/order/symmetry/atlas-coordinate cases. No runtime triangles, texture resources, draw groups or requests are added. The full 29-term inventory and missing native representations remain authoritative.

Local reproducible evidence: `scripts/review-ta98-female-coronary-registration.py`, `scripts/review-ta98-female-coronary-crossings.py`, `scripts/render-ta98-female-coronary-crossing-sections.py`; reports under `reports/ta98-completeness/female-coronary-body-registration.json`, `female-coronary-source-crossings.json`, and `female-coronary-source-crossing-sections.json`. Images: `outputs/ta98-female-coronary-body-registration.png` and `outputs/ta98-female-coronary-source-crossing-sections.png`. These local files keep their existing ignore status.

Next: establish matched root/annular and native myocardial surface constraints; preserve the whole named arterial assembly while repairing the source self-crossings; validate every junction, source-scope difference, whole native cardiac/venous/pulmonary packing and actual cuts before installation. New verified coronary geometry is still missing; anatomical acceptance remains unproven.


## Local uncrossing candidate results

### Follow-up triangular-gap trial definition

Original vertex/index review identifies the three-edge boundary implicated in the five remaining ventricular intersections as edges (78,81), (78,150), (81,150). Its longitudinal extent is approximately 1.294 mm in the empirical trial frame, with much smaller transverse extents. It is a triangular gap in the rendered local vessel surface; no anatomical opening or tissue/lumen identity is assigned. Its three coordinates share no exact positions with any of the other 22 source arterial/cardiac-guide objects. The two separate six-edge source boundaries remain distinct and fixed. Nearest sampled source guide vertices do not establish surface clearance or authorize native fitting.

For the follow-up local candidate only, permit the three explicitly listed triangular-gap coordinates to move with the same 0.2 mm bound and their local source normals to be transported. Keep the three boundary edges and original open topology; add no closing face. Freeze every other boundary and cross-object shared coordinate/normal. All original face indices, UVs, component identities and loop connectivity must remain exact; local face reversals, new degeneracy or increased whole-assembly crossings reject the candidate. Preserve the earlier fixed-boundary candidate separately. The unresolved triangular gap remains a visible source defect even if its crossings are removed, and prevents complete anatomical acceptance.

A separate local candidate removes the two atrial-circumflex self-intersection pairs at the original fork. It moves 22 distinct packed coordinates, at most 0.199902 mm in the empirical trial frame. All 2,534 atrial source faces remain. Source split normals are locally transported through old/new face rotation; indices, UVs, source materials, boundaries and cross-object seams remain exact. Sixty-five incident triangles change; their total surface area is 96.6933% of the source value and minimum individual area ratio is 56.1271%. These values quantify shape change, not validated arterial caliber or lumen anatomy.

The right-coronary ventricular source stays byte-exact: all 60 local smoothing attempts retain its five original self-crossing pairs, so no trial is adopted for that object. Three implicated boundary coordinates remain fixed. Removing or moving that boundary requires explicit source-interface review; the boundary is not assumed to be a free anatomical vessel endpoint.

Independent candidate replay verifies all 23 objects and 77 arrays; only the atrial position/normal arrays change. All original face indices and UVs, absent arrays, exact-coordinate duplicates, boundary and cross-object shared coordinates/normals pass. No face reversals or topology/component/boundary changes arise. All 15 full arterial self-checks and all 105 distinct arterial pair checks retain five self-crossing pairs and zero cross-object pairs. Numerical narrow phase still passes 864 permutation/coordinate cases plus exhaustive brute-force versus bounded sweep fixtures.

Actual original/candidate solid front/side views with transported normals and all seven original defect-patch front/side/transverse-section comparisons were inspected. The atrial crossing patch improves, but the ventricular crossing defects remain visible; this is not a complete repair, native registration, lumen certificate or anatomical acceptance. Source/native body-frame mismatch and all-neighbor cardiac packing remain unresolved. No geometry is installed in either main body, no runtime resources are added, and the full 29-term inventory remains unfulfilled.

Local candidate: `.local-models/ta98-review/female-coronary-ca2024-local-uncross-candidate.npz`, 567,331 bytes, SHA256 `e5e0ca84cf06db871e832a030a8229d77b9c254ee0a17264c93c0a86d95be564`. Original archive remains SHA256 `7d56091ee6890b4052d08a98aeef639c9de772abb356bfe229db8dedddf8e2e4`. Evidence scripts: `repair-ta98-female-coronary-source-local.py`, `validate-ta98-female-coronary-local-candidate.py`, `review-ta98-female-coronary-crossings.py --candidate`, `render-ta98-female-coronary-crossing-sections.py --candidate`, and `render-ta98-female-coronary-local-surface.py`. Local reports are `female-coronary-local-uncross-candidate.json`, `female-coronary-local-candidate-validation.json`, and `female-coronary-candidate-crossings.json`. Images are `outputs/ta98-female-coronary-candidate-crossing-sections.png` and `outputs/ta98-female-coronary-local-surface-comparison.png`. Existing ignored paths stay ignored.


## Triangular-gap follow-up candidate results

The follow-up trial removes the five ventricular source self-intersection pairs while retaining the earlier atrial correction. Sixteen distinct ventricular coordinates move, at most 0.199896 mm in the empirical trial frame. Only the three explicitly identified triangular-gap boundary coordinates are released; the gap stays open, with its same three edges and source vertex connections. The two six-edge boundaries, other object interfaces, all original faces/UVs/materials and original component identities remain exact. Source split normals are locally transported; no face is inverted, removed, added or capped. The earlier fixed-boundary candidate and original GLB/archive remain intact.

All 15 entire coronary source meshes now have zero detected proper/coplanar self-intersection pairs in this follow-up candidate. All 105 separate arterial pair checks also have zero detected pairs. Independent attribute replay verifies all 23 objects/77 arrays, with only the two repaired position/normal arrays per vessel changed. It separately verifies the explicitly permitted open triangular loop, every other boundary/cross-object coordinate and normal, original faces/UVs, exact-coordinate duplicates, 0.2 mm displacement bounds and positive old/new face orientation.

The two changed whole vessel meshes were also compared before/after against every one of the eight fixed source cardiac/root guides (16 comparisons, original and candidate). Atrial-circumflex/left-atrial surfaces retain 355 intersecting triangle pairs; every other tested guide pair remains at zero. No guide count increases. These existing intersections require tissue-scope/contact review and do not establish tissue-volume separation or source/native packing. The open/cut reference shells cannot be used as solids for clipping.

Actual same-camera solid original/candidate ventricular front/side views and all seven original defect patches with common original-position front/side/transverse cut planes were rendered and inspected. The folded surface improves; the open triangular gap remains unresolved. Numerical uncrossing does not prove a wall, lumen, ostium, correct vessel caliber, anastomosis or complete anatomical repair. Native body-frame/ventricular mismatch and whole cardiac/venous/pulmonary packing remain unresolved. None of this candidate is installed; accepted TA98 coverage is still unproven. The full 29-term inventory stays in scope.

Candidate archive: `.local-models/ta98-review/female-coronary-ca2024-triangular-loop-candidate.npz`, 567,328 bytes, SHA256 `91501b91ee038c6b455fdb1ebdb57ef9f348b13644597e25ce7b87e9c27f00a1`. Reproduce with repair, attribute validation, full crossing review and section-render scripts using `--triangular-loop`; use `check-ta98-female-coronary-candidate-guides.py` for the 16 fixed-guide comparisons and `render-ta98-female-coronary-local-surface.py --triangular-loop --ventricular` for the solid patch comparison. Local reports: `female-coronary-triangular-loop-candidate.json`, `female-coronary-triangular-loop-validation.json`, `female-coronary-triangular-loop-crossings.json`, `female-coronary-triangular-loop-guides.json`. Actual render files: `outputs/ta98-female-coronary-triangular-loop-ventricular-surface-comparison.png` and `outputs/ta98-female-coronary-triangular-loop-crossing-sections.png`. No runtime requests, geometry, textures or material groups are added; ignored files keep their existing status.

Next: characterize the remaining source gap and named source junctions, constrain registration from actual aortic root/annular and native myocardial surfaces, and review the complete native cardiac assembly in sections. Preserve variants/missing coronary terms and do not infer complete tissue coverage from a zero-intersection count.


## Complete left-source seam candidate results

The three complete six-point boundary cycles are now reconciled in a separate local candidate: left main/anterior interventricular, left main/circumflex, and circumflex/left marginal. Each is mutually one-to-one, has opposite original boundary winding and retains all source faces. Common float32 coordinates and normalized endpoint normals remove the tiny numeric cracks and independent endpoint normal jumps. Maximum trial-frame displacements are 0.082345, 0.067662 and 0.054193 micrometres respectively; each paired point and original packed vertex is recorded in the local report. Other positions/normals remain exact to the prior uncrossed candidate. No bridge, filler, wall, cap, subdivision or texture is invented.

The combined four-object left source exterior contains the same 1,570 triangles and changes from three exact-coordinate components to one. Open edges decrease from 54 to 22, leaving four intentional/unresolved source boundary groups. Combined nonmanifold/winding/duplicate/degenerate defects and invalid vertex links are zero. These results demonstrate source surface seam consistency, not complete coronary lumen/wall anatomy or biological branch coverage. Isolated source part topology remains unchanged; parts remain separately selectable future records rather than one opaque merged mesh.

Independent replay verifies all 23 objects/77 arrays and each declared complete seam coordinate/normal. All undeclared attributes, original face indices, UVs and material references remain exact, with no face reversals. Full 15-mesh self-crossing and all 105 arterial pair checks remain at zero detected pairs. All 32 before/after checks of the four changed entire arteries against eight fixed source cardiac/root guides retain their previous counts: left main/aortic root 115, anterior interventricular/left-ventricular guide 9, circumflex/left atrial guide 4 and circumflex/left-ventricular guide 74; all other tested pairs zero. These source-shell intersections still need contact/tissue-scope review; no native myocardial packing is certified.

Actual same-view before/after source surfaces were rendered with depth-tested barycentric interpolation of the stored source/candidate vertex normals, and principal-plane cuts through all three matched cycles inspected. The normals no longer jump at matched endpoints. Diagnostic principal planes are numerical source-cycle planes, not confirmed anatomical annular landmarks. The repaired joins correspond to source left-coronary/anterior-interventricular/circumflex/marginal routes (TA98 .201/.202/.206/.209); source branches remain adjoining tissue, not physical components contained within parent vessels. Native relationships and camera behavior are untouched until a correctly fitted native assembly is available.

Candidate: `.local-models/ta98-review/female-coronary-ca2024-left-seam-candidate.npz`, 567,330 bytes, SHA256 `a584183765b1fad951b7d90f2d43e30fccc8560e0fb8c93ed1cf2a84e4cbc419`. Prior candidates and original source remain intact. Reproduction: `join-ta98-female-left-coronary-source-seams.py`, `validate-ta98-female-coronary-left-seams.py`, crossing and guide scripts with `--left-seams`, and `render-ta98-female-left-coronary-seams.py`. Reports: `female-coronary-left-seam-candidate.json`, `female-coronary-left-seam-validation.json`, `female-coronary-left-seam-crossings.json`, `female-coronary-left-seam-guides.json`, `female-coronary-left-seam-render.json`. Render: `outputs/ta98-female-left-coronary-source-seams.png`. No new runtime geometry, requests, texture resources or material groups; existing ignore policy retained.

Remaining: actual native aortic-root/annular and myocardial registration; the unresolved three-edge ventricular gap; partial/other source branch/root interfaces; source tissue contact and full native cardiac/venous/pulmonary packing, walls/lumens and the entire 29-term coverage. No anatomical acceptance or native installation is claimed from the improved seams.

## Native aortic attachment review and constrained initialization

The fixed native female aorta/ventricular pair shares 24 exact packed coordinates. Its ventricular attachment is a 24-triangle annular strip with two twelve-edge closed rims. The aortic object includes 22 of those same triangles, oppositely oriented; two ventricular strip triangles have no identical aortic triangle. This is a numerical source attachment, not proof of the anatomical valve annulus, an open ostium, or a complete wall. The separate aortic and pulmonary valve objects share no exact coordinates with this attachment. Preserve all native surfaces while establishing correspondence.

The source Root of Aorta has two closed basal boundary cycles (26 and 43 edges) and two noncycle upper boundary groups. The source Aortic Fibrous Ring is a two-component open shell with four boundary cycles; it is not a single planar ring. Actual native attachment and valve cuts, and source root/ring/left- and right-coronary cuts, must be reviewed before choosing landmarks. Source ring and native valve scopes are not interchangeable.

Constrained initialization plan: retain the empirical body scale and orient each numerical basal-rim principal normal consistently with body Y. Apply only the minimal proper rotation between the source outer basal-rim normal and the native outer attachment-rim normal, followed by center translation. Lock rotation about that normal to the body-frame initialization; do not resize a root from bounds or deform any vessel. This is an uninstalled diagnostic trial, not an asserted anatomical registration. Preserve all source parts, faces, UVs, normals and seams; the transform must be common to the whole source assembly. Compare both rims and actual whole myocardial/arterial surfaces and sections. Reject installation if myocardial extent, ostial routes, neighboring tissues, missing terms or source defects remain unresolved. Keep the earlier body fit and source archives intact.

2026-10-03 result: the open Blender 5.2.1 session accepts 41 diagnostic objects in three separate review scenes. The original scene and original object pointers are preserved; no existing main file is opened or saved. Actual native and source surfaces are inspected there, with flat diagnostic shading rather than an assertion of source material finish. Four numerical source/native root-interface cuts are also inspected in `outputs/ta98-female-aortic-interface-sections.png`.

The minimal rotation trial is 18.7233 degrees with zero numerical rim-center residual. Source rim to nearest native rim **vertex** distances are 0.660–5.174 mm, median 2.352 mm; unequal polygon sampling means these are not surface-fit errors. Source ventricular lower Y extent remains 24.783 mm above native lower Y extent after the trial (different wall/whole-organ scopes). Actual common-camera arterial surfaces and horizontal/sagittal sections show source arteries still floating away from, or hidden within, the native ventricular surface. Center alignment is insufficient: the trial fails installation. Next identify homologous whole myocardial/sulcal surfaces and the root/valve scope before proposing a coupled fit; do not simply project all vessel vertices onto an arbitrary nearest heart face.

The native fixed guides retain existing defects: aorta 2 degenerate faces / 2 nonmanifold edges; ventricles 6 / 6 and a three-edge genuine open cycle plus two collapsed self-edge boundary artifacts; aortic valve 10 / 10; pulmonary valve 6 / 7 and one duplicate face. These exact-weld figures include degenerate-face effects and require localized review rather than blind hole filling. No native repair or anatomical acceptance is claimed.

Independent replay verifies all four native packed position/face arrays exactly; all 24 diagnostic source objects receive the same proper rigid transform with at most 1e-15 m numerical replay tolerance, all faces exact and face vectors rotating without inversion. It verifies 22 opposite native attachment faces and the unchanged preferred coronary candidate archive hash. Evidence: local `female-aortic-interface-review.json`, `female-aortic-blender-import.json`, `female-aortic-interface-sections.json`, `female-coronary-root-initialization.json`, `female-coronary-root-initialization-render.json`, and `female-aortic-root-review-validation.json`. Common-camera review: `outputs/ta98-female-coronary-root-initialization.png`. No runtime meshes, triangles, materials, textures or requests added; source archives and intentionally ignored files remain local. The complete 29-term goal is still unfulfilled.
