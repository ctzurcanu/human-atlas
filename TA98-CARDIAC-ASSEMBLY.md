# TA98 cardiac assembly and source inventory

## Rendered placement correction — 22 October 2026

The previous packing audits used the installed part transforms, while the main
viewer still used stale transforms embedded in materials. For female epicardium,
the material scale was 1.0352162624 and the installed part scale 0.6633665810,
with different rotations. Thus those earlier zero-overlap checks did not describe
what the viewer actually displayed. Replaying the legacy rendered finite tissue
gives 74.293549 mL overlap with the left lung and 9.378201 mL with the right.

Main-view materials now bind the same part transform used by picking, sections,
the anatomical labels and the geometry tools. The material-sharing key separates
different poses, while identical shader programs remain reusable. Installed mesh
banks, all part records, skeleton, lungs and all 113 term controls are unchanged.
The corrected female finite tissue has 0 / 0 mL lung overlap; every owned cardiac
part, including the roots, has zero detected lung or chest-wall crossing pairs.
Evidence: `female-installed-lung-audit-v678`, `female-legacy-render-overlap-v679`,
`female-correct-render-lungs-v680`, `female-correct-render-chest-v681` under the
local heart workbench; `server/organ-render-pose.test.mjs` replays the actual
indexed epicardial vertices through GPU placement uniforms against audited bounds.
These checks do not establish clinical sizing or microscopic anatomical acceptance.

## Installed thoracic packing correction — 21 October 2026

Both ordinary local TA98 catalogues now contain the reviewed corrections:

| Model revision | Heart-body overlap with left/right lung | Minimum epicardium–sternum gap |
| --- | ---: | ---: |
| `cardiac-male-clear-thorax-v525` | 0 / 0 mL | 10.729340 mm |
| `cardiac-female-clear-thorax-v493` | 0 / 0 mL | 21.864092 mm |

These occupied-tissue checks cancel only literally opposed internal facets
and reject same-facing duplicates. They test closed finite heart-body tissue,
including nested overlap, rather than relying on surface crossing counts.
The previous male v382 overlaps were 24.921961 / 14.108489 mL and female v426
overlaps were 2.395558 / 14.941997 mL. The whole physical cardiac assembly now
has zero detected surface crossings with the active native sternum, ribs,
chest-wall muscles or diaphragm in each model. Cardiac walls also clear the
bronchial source tree. Every skeletal record remains literal.

The correction is **not complete great-vessel acceptance**. Separate coarse
vascular references still cross pulmonary territories: male pulmonary trunk,
68 face pairs; female aortic arch and pulmonary trunk, 378 face pairs. These
are not overlaps of the finite heart-body tissue measured above. Exact distal
vessel joins and clinical size remain unverified. The male left apical source
also retains its original 52 crossing face pairs with left longus colli; no
new lung/rib or lung/chest-muscle crossings were introduced. Existing female
pulmonary/chest source defects were not changed by its cardiac pose correction.

### Male pulmonary space and shared-boundary correction

Uniform enlargement was rejected: a 1.35-scale trial created 52,226 thoracic
crossing face pairs. A conservative local halfspace proposal limited uniform
growth to 1.021193366, but its actual triangles still failed packing. Neither
trial was installed. The v474 global heart search reached its 90-second
deadline and produced no installable candidate; its runtime evidence remains
`570cd72122b14a05a9a7ba9433a15a78.json`.

The installed solution reshapes the lower pulmonary assembly with one smooth
field shared by 19 territories and 27 connected bronchial pieces. Below
Y 1.42 m the field uses 1.20 transverse and 0.90 depth scale, with −5 mm X/Z
translation; it transitions to identity at Y 1.475 m. A monotone vertical
field raises the middle region by at most 10 mm. All 2,286 upper source
controls remain literal. This is a packing deformation, not a physiological
lung-growth measurement. Summed model envelope volumes change from
668.626 / 792.802 mL left/right to 663.900 / 807.550 mL.

The continuous field initially caused 15 discrete coarse/fine chart crossings.
A common-pool repair moves five internal controls by at most 30.001 µm,
pinning every outer lung control and preserving every source face, UV and
normal chart. Independent audits then pass all 19 closed cells, all 171
occupied-volume/disjointness pairs, and 42 reciprocal chart interfaces.
The bronchial source remains one connected component, its 859 source open
edges are retained, and all 20 named segment destinations pass. This does not
establish finite airway walls or histological territory boundaries.

The male heart uses a proper common similarity of scale 0.8061126277 relative
to the original shared cardiac bank; its internal tissue interfaces are literal.
The female uses scale 0.7370739788 relative to that same bank. Her duplicate
long coarse inferior-caval reference is retained only as an open cardiac cuff
above Y 1.2605 m; the original hepatic/caval source below it stays in place.
No cap or invented finite vascular wall closes that reference mouth.

All evidence below is under `.local-models/geometry-workbench/heart/`:

- Male: `male-pulmonary-shared-chart-repair-v517`,
  `male-pulmonary-cell-quality-v518`, `male-pulmonary-owned-check-v519`,
  `male-final-heart-lung-volume-v520`, `male-final-cardiac-chest-v521`,
  `male-final-cardiac-bronchi-v522`, `male-final-cardiac-roots-v523`,
  `male-final-pulmonary-chest-v524`, `male-cardiopulmonary-final-install-v525`.
- Female: `female-thoracic-caval-cuff-v479`,
  `female-thoracic-clear-chest-v480`, `female-thoracic-clear-bronchi-v481`,
  `female-thoracic-clear-lungs-v482`, `female-thoracic-clear-volume-v492`,
  `female-thoracic-final-install-v493`.

The male guarded job passed 84 changed-part audits and 162 interface/clearance
checks in 39.94 seconds; local installation completed in 59.30 seconds. Actual
atlas exterior and coronal-section images are retained with the reviewed job.
Rollback manifests are in `.local-models/ta98-review/cardiac-male-clear-thorax-v525/`
and `cardiac-female-clear-thorax-v493/`. No external deployment was performed.

## Local installation — 21 October 2026

The male `cardiac-male-clear-thorax-v525` and female `cardiac-female-clear-thorax-v493`
revisions provide **41 parts, 49 controls and 42 directly mapped A12.1 terms**
in each model. Two additional controls use parent/partial placement; together
the controls refer to 44 distinct authority terms. **71 of the 113 A12.1 terms
still lack a distinct exact representation**. Mapping is not anatomical
acceptance, and collective controls are not extra meshes.

Installed scope: four chamber territories; interatrial and atrioventricular
septa; muscular and membranous interventricular septum; finite endocardium,
myocardium and epicardium; four valve groups with 11 source leaflets; five
papillary parts; five separately selectable native chordal fans; two finite
annular support regions at measured AV leaflet hinges; and eight
retained or placed great-vessel/root parts. The
Systems tree groups these under Heart. Canonical TA98 parents are preserved.
The whole-heart control includes the same real parts as its peeling controls.
Retired duplicate chamber/valve sources are suppressed and retained locally.
Skeleton and diaphragm geometry remain unchanged; the male pulmonary packing
correction is described above.

The fourteen wall/pulmonary-root cells, fifteen atrioventricular leaflet/chordal/
papillary cells, two annular support regions and the finite SVC use one
**24,210,874-byte geometry bank**, shared by both models. All 32 exported
Float32 cells are closed, without detected self-crossings or occupied-volume
overlaps. There are 117 reciprocal part-pair interfaces, including source
leaflet-to-host, annular-to-leaflet/myocardial and all five papillary-base
contacts. Each part-pair certificate includes all its subdivided charts.
Explicit ownership preserves the source union. Leaflet coaptation/sliding
contacts remain adjacent rather than fused. Papillary and annular regional
extents remain provisional, including one small native RV papillary patch at
the outer coat. Source controls remain literal; only newly constructed cut
controls were rounded, at most 19.813791 nm against their immutable construction
reference, with a global crossing guard. Cut-corner normal/UV fields retain
source-face interpolation. No new image assets are needed.

Six named exterior controls (base, apex and four surfaces) select actual
source-bound outer facets. These are surface footprints with provisional
boundaries, not duplicate tissue solids. The two new collective controls
group the existing atria and ventricles. Source-qualified LV papillary
homology and complete fine anatomy remain unresolved.

The right pulmonary cusp's inherited seven-edge breach is closed by five
triangles using its existing controls. Original faces, coordinates and
attributes remain unchanged. This patch closes finite leaflet tissue, not the
pulmonary valve orifice. Its 70,331-byte bank is also shared. The mixed-material
normal-sign fallback now preserves normals on ordinary meshes without a signed
attribute, fixing their previously flat grey appearance.

The ten finite wall cells and four finite pulmonary-vein roots previously used
the **11,979,849-byte bank**, now superseded by the coupled bank above.
They have closed oriented boundaries and zero detected self-crossings. All
occupied-volume pairs with overlapping bounds pass (63 male, 61 female);
52 declared tissue interfaces have reciprocal coordinates and normals.
A separate shared-facet/winding certificate resolves
the Boolean engine's small shared-boundary overlap classification. Native
residuals and the certificate are retained in the evidence; tolerances were
not relaxed. Endocardial/epicardial 40 µm target thickness is illustrative;
the initial coat construction left 1.678% of difficult source surface as
uncoated myocardium. The original
anatomical carrier controls are preserved. Calculation tetrahedra are local
construction data and are not viewer assets.

Female placement reuses the literal shared bank with one positive uniform
similarity, currently scale 0.7370739788 relative to the original bank. The
original landmark-only scale 1.0352162624 was too large and anterior for the
recipient chest. The v426 scale 0.8281730099 still intersected lung tissue;
the installed v493 pose includes an additional 0.89 relative scale and proper
rotation/translation. It does not re-encode finite cardiac vertices or normals. The
native recipient aortic arch is retained; the four finite pulmonary-vein roots
share the same cardiac placement. This fit is not a clinical size
measurement. Distal recipient vessel joins remain provisional.

Remaining anatomical limits: regional septal extents; fibrous annular
histological extent/discontinuities and
fine attachment anatomy beyond the checked macroscopic contacts; complete fine chordal distribution;
the unresolved inferior/lateral LV papillary homology; full left
anterior papillary extent; detailed valve patency and distal recipient joins; native
vascular-source self-crossings; conduction tissues, finer valve
features and pericardial anatomy. No claim is made that all 113 terms or the
entire circulatory system are finished. Actual exterior, peeling and coronal
section views were inspected before installation; numerical passes alone did
not establish anatomical acceptance.

Reproducible evidence and rollback:

- `.local-models/geometry-workbench/heart/annular-finished-bank-v378/`
- `.local-models/geometry-workbench/heart/annular-male-install-v382/`
- `.local-models/geometry-workbench/heart/annular-female-install-v384/`
- `.local-models/ta98-review/cardiac-male-annuli-v382/manifest.json`
- `.local-models/ta98-review/cardiac-female-annuli-v384/manifest.json`
- `.local-models/geometry-workbench/heart/apparatus-owned-bank-v242/`
- `.local-models/geometry-workbench/heart/coupled-attachment-bank-v256/`
- `.local-models/geometry-workbench/heart/coupled-male-install-v260/`
- `.local-models/geometry-workbench/heart/coupled-female-install-v262/`
- `.local-models/ta98-review/cardiac-male-attachments-v260/manifest.json`
- `.local-models/ta98-review/cardiac-female-attachments-v262/manifest.json`
- `.local-models/geometry-workbench/heart/apparatus-owned-male-install-v244/`
- `.local-models/geometry-workbench/heart/apparatus-owned-female-install-v246/`
- `.local-models/ta98-review/cardiac-male-apparatus-v244/manifest.json`
- `.local-models/ta98-review/cardiac-female-apparatus-v246/manifest.json`
- `.local-models/geometry-workbench/heart/pulmonary-cusp-patch-v223/`
- `.local-models/geometry-workbench/heart/pv-root-male-install-v199/`
- `.local-models/geometry-workbench/heart/chordal-male-install-v218/`
- `.local-models/geometry-workbench/heart/chordal-female-install-v220/`
- `.local-models/geometry-workbench/heart/chordal-finite-candidate-v214/`
- `.local-models/ta98-review/cardiac-male-chordal-controls-v218/manifest.json`
- `.local-models/ta98-review/cardiac-female-chordae-v220/manifest.json`
- `.local-models/geometry-workbench/heart/pv-root-female-install-v201/`
- `.local-models/geometry-workbench/heart/pv-root-bank-v196/`
- `.local-models/ta98-review/cardiac-{male-roots-v199,female-roots-v201}/manifest.json`
- `.local-models/geometry-workbench/heart/male-install-v106/`
- `.local-models/geometry-workbench/heart/male-controls-install-v115/`
- `.local-models/geometry-workbench/heart/female-install-v113/`
- `.local-models/geometry-workbench/heart/regional-bank-v102/`
- `.local-models/geometry-workbench/heart/female-instance-v112/`
- `.local-models/ta98-review/cardiac-{male-local-v106,female-local-v113}/manifest.json`

The source inventory below was recorded before reconstruction against the
complete authoritative A12.1 section, including pericardial terms whose
original parents sit outside the whole-heart root. Source columns retain that
baseline. The final column records the installed representation status in both
models. Local donor redistribution rights remain unverified; these revisions
were installed locally and were not published.

| TA98 | Term | Original parent | Representation class | Male sources | Female sources | Missing active exact M/F | Installed male/female control |
|---|---|---|---|---|---|---|---|
| A12.1.00.001 | heart | A12.0.00.000 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.00.002 | base of heart | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Selectable source-facet exterior footprint; provisional boundary, no extra solid |
| A12.1.00.003 | anterior surface | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Selectable source-facet exterior footprint; provisional boundary, no extra solid |
| A12.1.00.004 | diaphragmatic surface | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Selectable source-facet exterior footprint; provisional boundary, no extra solid |
| A12.1.00.005 | right pulmonary surface | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Selectable source-facet exterior footprint; provisional boundary, no extra solid |
| A12.1.00.006 | right border | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.007 | apex of heart | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Selectable source-facet exterior footprint; provisional boundary, no extra solid |
| A12.1.00.008 | notch of cardiac apex | A12.1.00.007 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.009 | anterior interventricular sulcus | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.010 | posterior interventricular sulcus | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.011 | coronary sulcus | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.012 | ventricle | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | LOCAL:female:ventricles:1506 (parent, native) | True/True | Collective control of the stored named chamber tissues; no extra solid |
| A12.1.00.013 | interventricular septum | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | HRA:VH_F_interventricular_septum (exact, reference, hidden) | True/True | Stored geometry / collective control; provisional extents |
| A12.1.00.014 | muscular part | A12.1.00.013 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.00.015 | membranous part | A12.1.00.013 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.00.016 | atrioventricular septum | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.00.017 | atrium | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Collective control of the stored named chamber tissues; no extra solid |
| A12.1.00.018 | auricle | A12.1.00.017 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.019 | interatrial septum | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.00.020 | trabeculae carneae | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.021 | vortex of heart | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.022 | papillary muscles | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Inferior papillary muscle of left ventricle (parent, native); ZA:Inferior papillary muscle of right ventricle (parent, native) | missing | True/True | Parent-only / partial source tissue |
| A12.1.00.023 | chordae tendineae | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Five native finite chordal fans and collective control; reciprocal neck/insertion interfaces; finer distribution provisional |
| A12.1.00.024 | false chordae tendineae | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.025 | right fibrous trigone | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.026 | left fibrous trigone | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.027 | right fibrous ring | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Finite owned annular support at measured leaflet hinges; reciprocal interfaces, histological extent provisional |
| A12.1.00.028 | tendon of infundibulum | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.029 | tendon of valve of inferior vena cava | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.030 | triangle of sinu-atrial node | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.00.031 | facies pulmonalis sinistra | A12.1.00.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Selectable source-facet exterior footprint; provisional boundary, no extra solid |
| A12.1.00.032 | left fibrous ring | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Finite owned annular support at measured leaflet hinges; reciprocal interfaces, histological extent provisional |
| A12.1.01.001 | right atrium | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Right atrium (exact, native) | LOCAL:female:r_atrium:1512 (exact, native) | False/False | Stored geometry / collective control; provisional extents |
| A12.1.01.002 | right auricle | A12.1.01.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.003 | crista terminalis | A12.1.01.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.004 | openings of smallest cardiac veins | A12.1.01.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.005 | fossa ovalis | A12.1.01.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.006 | limbus fossae ovalis | A12.1.01.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.007 | foramen ovale | A12.1.01.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.008 | musculi pectinati | A12.1.01.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.009 | opening of coronary sinus | A12.1.01.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.010 | opening of inferior vena cava | A12.1.01.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.011 | opening of superior vena cava | A12.1.01.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.012 | sinus of venae cavae | A12.1.01.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.013 | sulcus terminalis cordis | A12.1.01.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.014 | intervenous tubercle | A12.1.01.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.015 | valve of inferior vena cava | A12.1.01.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.01.016 | valve of coronary sinus | A12.1.01.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.001 | right ventricle | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Right ventricle (exact, native) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.002 | right atrioventricular orifice | A12.1.02.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.003 | tricuspid valve | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Inferior leaflet of right atrioventricular valve (parent, native) | LOCAL:female:tricuspid_atrioventricular_valve:1503 (exact, native); HRA:VH_F_tricuspid_valve (exact, reference, hidden) | True/False | Stored geometry / collective control; provisional extents |
| A12.1.02.004 | anterior cusp | A12.1.02.003 | Named whole or regional tissue; extent, interfaces and anatomy unverified | FJ2421 (exact, reference, hidden) | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.02.005 | posterior cusp | A12.1.02.003 | Named whole or regional tissue; extent, interfaces and anatomy unverified | FJ2433 (exact, reference, hidden) | missing | True/True | Source-qualified posterior identity; attachments provisional |
| A12.1.02.006 | septal cusp | A12.1.02.003 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Septal leaflet of right atrioventricular valve (exact, native); FJ2436 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.007 | supraventricular crest | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.008 | conus arteriosus | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.009 | opening of pulmonary trunk | A12.1.02.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.010 | pulmonary valve | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | LOCAL:female:pulmonary_semilunar_valve:1505 (exact, native); HRA:VH_F_pulmonary_valve (exact, reference, hidden) | True/False | Stored geometry / collective control; provisional extents |
| A12.1.02.011 | right semilunar cusp | A12.1.02.010 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Right semilunar leaflet of pulmonary valve (exact, native); FJ2434 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.012 | left semilunar cusp | A12.1.02.010 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Left semilunar leaflet of pulmonary valve (exact, native); FJ2427 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.013 | anterior semilunar cusp | A12.1.02.010 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Anterior semilunar leaflet of pulmonary valve (exact, native); FJ2417 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.014 | nodules of semilunar cusps | A12.1.02.010 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.015 | lunules of semilunar cusps | A12.1.02.010 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.016 | commissures of semilunar cusps | A12.1.02.010 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.017 | anterior papillary muscle | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Anterior papillary muscle of right ventricle (exact, native) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.018 | posterior papillary muscle | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | FJ2430 (exact, reference, hidden) | HRA:VH_F_papillary_muscle_of_heart_pos (exact, reference, hidden) | True/True | Source-qualified posterior identity; attachments provisional |
| A12.1.02.019 | septal papillary muscle | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Septal papillary muscle of right ventricle (exact, native) | HRA:VH_F_papillary_muscle_of_heart_med (exact, reference, hidden) | False/True | Stored geometry / collective control; provisional extents |
| A12.1.02.020 | septomarginal trabecula | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.02.021 | trabeculae carneae | A12.1.02.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.03.001 | left atrium | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Left atrium (exact, native) | LOCAL:female:l_atrium:1507 (exact, native) | False/False | Stored geometry / collective control; provisional extents |
| A12.1.03.002 | left auricle | A12.1.03.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.03.003 | musculi pectinati | A12.1.03.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.03.004 | openings of pulmonary veins | A12.1.03.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.03.005 | valve of foramen ovale | A12.1.03.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.001 | left ventricle | A12.1.00.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Left ventricle (exact, native) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.04.002 | left atrioventricular orifice | A12.1.04.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.003 | mitral valve | A12.1.04.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | LOCAL:female:l_bicuspid_atrioventricular_valve:1502 (exact, native); HRA:VH_F_mitral_valve (exact, reference, hidden) | True/False | Stored geometry / collective control; provisional extents |
| A12.1.04.004 | anterior cusp | A12.1.04.003 | Named whole or regional tissue; extent, interfaces and anatomy unverified | FJ2420 (exact, reference, hidden) | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.04.005 | posterior cusp | A12.1.04.003 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Posterior leaflet of left atrioventricular valve (exact, native); FJ2432 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.04.006 | commissural cusps | A12.1.04.003 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.007 | anterior papillary muscle | A12.1.04.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | HRA:VH_F_papillary_muscle_of_heart_ant (exact, reference, hidden) | True/True | Parent-only / partial source tissue |
| A12.1.04.008 | posterior papillary muscle | A12.1.04.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.009 | aortic vestibule | A12.1.04.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.010 | aortic orifice | A12.1.04.001 | Opening; never a solid plug | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.011 | trabeculae carneae | A12.1.04.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.012 | aortic valve | A12.1.04.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | LOCAL:female:aortic_semilunar_valve:1504 (exact, native); HRA:VH_F_aortic_valve (exact, reference, hidden) | True/False | Stored geometry / collective control; provisional extents |
| A12.1.04.013 | right semilunar cusp | A12.1.04.012 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Right coronary leaflet (exact, native); FJ2435 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.04.014 | left semilunar cusp | A12.1.04.012 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Left coronary leaflet (exact, native); FJ2426 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.04.015 | posterior semilunar cusp | A12.1.04.012 | Named whole or regional tissue; extent, interfaces and anatomy unverified | ZA:Non-coronary leaflet (exact, native); FJ2431 (exact, reference, hidden) | missing | False/True | Stored geometry / collective control; provisional extents |
| A12.1.04.016 | nodules of semilunar cusps | A12.1.04.012 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.017 | lunules of semilunar cusps | A12.1.04.012 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.04.018 | commissures of semilunar cusps | A12.1.04.012 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.05.001 | endocardium | A12.1.00.001 | Tissue layer; independent volume and shared interfaces required | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.06.001 | myocardium | A12.1.00.001 | Tissue layer; independent volume and shared interfaces required | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.06.002 | conducting system of heart | A12.1.06.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.06.003 | sinu-atrial node | A12.1.06.002 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.06.004 | atrioventricular node | A12.1.06.002 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.06.005 | atrioventricular bundle | A12.1.06.002 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.06.006 | right bundle | A12.1.06.005 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.06.007 | left bundle | A12.1.06.005 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.06.008 | subendocardial branches | A12.1.06.004 | Fine structure or tissue; whole organ/valve mesh does not establish it | missing | missing | True/True | Not implemented distinctly |
| A12.1.07.001 | transverse pericardial sinus | A07.1.03.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.07.002 | oblique pericardial sinus | A07.1.03.001 | Boundary, region or cavity landmark; not established by a generic outer surface | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.001 | pericardium | A07.1.02.105 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.002 | fibrous pericardium | A12.1.08.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.003 | sternopericardial ligaments | A12.1.08.002 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.004 | bronchopericardial membrane | A12.1.08.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.005 | serous pericardium | A12.1.08.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.006 | parietal layer | A12.1.08.005 | Tissue layer; independent volume and shared interfaces required | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.007 | visceral layer | A12.1.08.005 | Tissue layer; independent volume and shared interfaces required | missing | missing | True/True | Stored geometry / collective control; provisional extents |
| A12.1.08.008 | serosa | A12.1.08.005 | Tissue layer; independent volume and shared interfaces required | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.009 | subserosa | A12.1.08.005 | Tissue layer; independent volume and shared interfaces required | missing | missing | True/True | Not implemented distinctly |
| A12.1.08.010 | fold of left vena cava | A12.1.08.001 | Named whole or regional tissue; extent, interfaces and anatomy unverified | missing | missing | True/True | Not implemented distinctly |

## Assembly interfaces

- Fix native chambers, valves and great vessels while establishing source correspondence.
- Keep left/right chambers distinct. A combined ventricular envelope does not independently represent either chamber.
- Leaflets are tissue; valve openings are empty boundaries, never plugs.
- Semilunar valves have no chordal/papillary apparatus.
- Review actual atrial/venous and chamber/valve/outflow surfaces, sections and lumen routes before accepting contacts.
- Do not project coronary tubes indiscriminately onto a nearest ventricular face; preserve sulcal routes, ostia and venous crossings.
- No source reference installation until geometry, scope, attributes, whole packing and local redistribution requirements pass.

## Native source identity review

Unmapped native cardiac records are retained in the JSON evidence. Bind any corrections to their source concept IDs. The female ventricular envelope remains a broader group; do not assign separate left/right chamber identities. Inferior/posterior papillary or tricuspid terminology needs source-qualified review rather than an unqualified synonym guess.

Primary terminology: [TA98 right ventricle](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.1.02%20TA98%20EN.htm), [TA98 left ventricle](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.1.04%20TA98%20EN.htm). Flow/valve context: [NIH NHLBI](https://www.nhlbi.nih.gov/health/heart/blood-flow).

See TA98-CORONARY-ARTERIES.md for all 29 ascending-aorta/coronary terms, original/default source archives, source repairs, failed registration trials and remaining acceptance requirements. Runtime geometry and materials remain fixed.

## Native cardiac identity and inflow corrections — 2026-10-03

Eighteen source-qualified identity rules resolve all previously unmapped native cardiac records: fourteen exact terms (eight male leaflet identities; six female valve/atrial identities) and four broader parent placements. The combined female ventricular envelope does not represent separate left/right chambers. The male inferior tricuspid leaflet and two inferior papillary sources remain unresolved against posterior terminology; their parent assignments do not certify that homology.

Named same-frame caval and pulmonary venous returns now have reciprocal before/after atrial links. The female catalogue adds superior/inferior caval inflow pairs (four directed links). Existing male routes already supplied these targets and remain deduplicated. Final reciprocal relation totals: male 16,062; female 4,116. Suppressed, donor, cross-sex and inappropriate-system targets are excluded. Named valve tissue and unresolved groups display their representation limits in the normal inspector.

Validation: ten targeted cardiac/coronary tests, TypeScript, component checks across all four catalogues, all 7,474 authority terms across eight hierarchy catalogues, reciprocal links and production build pass. Live female aortic valve shows TA98 A12.1.04.012; male left coronary leaflet shows A12.1.04.014. Showing related anatomy preserves byte-identical actual position, target, up, quaternion and projection offset in both live checks. Share-frame normalized coordinates may change with selection bounds without moving the camera; actual pose is recorded in reports/ta98-completeness/native-cardiac-live-camera.json.

Full audit: male active direct material mappings 1,222 / 5,910 (20.7%); female 608 / 5,903 (10.3%). This pass adds eight male and six female exact identities, not new meshes. All four master/runtime source hashes remain unchanged. No geometry, textures, accepted anatomical terms, commits or deployment added. Coronary fitting, myocardial/valvular interfaces, lumen continuity and main-body geometry acceptance remain unfinished.

## Pulmonary venous root extraction inventory — 21 October 2026

Only existing connected proximal source components are selected for the female heart. The remaining tributaries stay retained locally; a root is not counted as its complete vein. Complete source vein terminology and descendants:

| TA98 | Term | Representation |
|---|---|---|---|---|---|---|---|
| A12.3.02.001 | pulmonary veins | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.101 | right superior pulmonary vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.102 | apical vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.103 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.104 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.105 | anterior vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.106 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.107 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.108 | posterior vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.109 | infralobar part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.110 | intralobar part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.111 | middle lobe vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.112 | lateral part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.113 | medial part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.201 | right inferior pulmonary vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.202 | superior vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.203 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.204 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.205 | common basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.206 | superior basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.207 | anterior basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.208 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.209 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.210 | inferior basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.301 | left superior pulmonary vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.302 | apicoposterior vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.303 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.304 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.305 | anterior vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.306 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.307 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.308 | lingular vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.309 | superior part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.310 | inferior part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.401 | left inferior pulmonary vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.402 | superior vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.403 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.404 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.405 | common basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.406 | superior basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.407 | anterior basal vein | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.408 | intrasegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.409 | intersegmental part | Proximal connected source root only; distal subdivision not newly built |
| A12.3.02.410 | inferior basal vein | Proximal connected source root only; distal subdivision not newly built |

## Apparatus correction prepared — 21 October 2026

A current-frame replay exposed a 4.35943 mm vertical registration error in the
two added anterior AV cusps. Their earlier translation had not been composed
with the recorded 183 cm whole-body similarity. The candidate composes those
transforms algebraically and carries the identical source correction through
the existing female placement. It changes no vertex/normal words, adds no
geometry bank, and keeps all chamber walls, other leaflets, vessel roots, lungs,
diaphragm and skeletal parts fixed. Both anterior cusps remain finite source
tissues; their annular/chordal interfaces remain adjacent/unjoined until
independently established. This is a placement correction, not a fabricated
attachment or anatomical acceptance.

A source-qualified comparison identifies the native inferior tricuspid leaflet
with BodyParts3D FJ2433 (posterior leaflet; maximum all-vertex surface distance
0.3074 mm), and the inferior RV papillary source with FJ2430 (posterior papillary
muscle; 0.7339 mm). These are sampled source-shape distances, not whole-surface
Hausdorff or attachment proofs. The named source identities, complete source
frame composition and attitudinal terminology together support the two TA98
placements. The unresolved inferior/lateral LV papillary source is retained
as partial/parent anatomy. The native source shapes are not modified.

References inspected: [Anatomy and Physiology of the Tricuspid Valve, JACC](https://www.jacc.org/doi/10.1016/j.jcmg.2018.07.032)
(search-indexed author text states the conventional posterior/attitudinal
inferior correspondence; full publisher page returned HTTP 403), and the
[100-specimen tricuspid anatomy study](https://pubmed.ncbi.nlm.nih.gov/26133960/)
(abstract describes inferior leaflet, inferior papillary variation and septal
cord attachments). No source illustrations or proprietary mesh pixels copied.

Prepared candidates and numerical source replay:
`.local-models/geometry-workbench/heart/apparatus-correction-v123/` and
`apparatus-correspondence-v121/`. The earlier v120 comparison omitted the
documented height transform and is invalid; v122 had a female frame-label typo
and is superseded by v123. Neither superseded report authorizes installation.

The diagnostic root-gating attempt did not prove separated blood circuits.
Its native-plane gates fail to meet tissue at two pulmonary-vein roots; that
may reflect an anatomical gap or an unsuitable diagnostic plane and needs
actual section/contact review. No successful patency claim follows from the
six-handle topology alone. The complete installed wall also contains a tiny
12-face disconnected remnant (about 1.16e-7 mL); it remains an explicit repair
item rather than evidence of an extra anatomical structure.

The male apparatus correction is installed as `cardiac-male-apparatus-v124`;
the corresponding female revision is `cardiac-female-apparatus-v125`. The
independent replay checks every emitted cusp triangle against the composed
current-frame source transform: coordinate error stays below the 1.863 nm
analytical Float32 bank bound. Female placement replays the corrected male
triangles through the proper similarity within 5.31e-16 m. No added geometry
banks, textures or parts. Detailed valve attachment and complete circuit
patency remain open work; this revision does not establish them.

## Pulmonary-vein root wall and interface plan — 21 October 2026

The complete installed ten-cell wall is fixed, as are the skeleton, diaphragm,
lungs, chamber/septal boundaries and the other vessel roots. The four pulmonary
veins retain their actual source outer courses. Their branching source pieces
are united before deriving finite walls; overlapping source pieces must not
become overlapping tissue cells. The inner courses use the same explicitly
illustrative 0.8 radial allocation that established the existing atrial bores.
This is an effective wall allocation, not measured histological thickness.

The root tissue must stop at the atrial blood boundary, and must not continue
inside the left atrium. Existing myocardial occupancy takes precedence. Actual
cut fragments at the root/myocardium contact must be shared once with reciprocal
orientation and normals. End rims close wall thickness only; no disk may span
a blood channel. End disks used by Boolean envelopes or circuit diagnostics
are computational objects and are never installed or counted as anatomy.

Before partition or installation, inspect each complete root wall for connected
components, open edges, self-crossings and lumen continuity; check all wall/root
and root/root occupied volumes. After packing, repeat those checks on the actual
emitted Float32 bank and review the rendered outer surface and sections. Female
transfer uses the established positive similarity only after the male candidate
passes. Native distal recipient joins remain explicitly provisional until
measured. No circulatory work beyond the communicating cardiac roots is implied.

The earlier native-plane circuit diagnostic is superseded by measured lumen
contours in `heart/blood-spaces-v130/blood-spaces.json`. With computation-only
gates, the actual installed boundary separates right atrium/ventricle from left
atrium/ventricle and from the exterior. This proves that diagnostic separation,
not complete valve patency. Complete curve-to-wall distance in
`heart/root-mouths-v127/` and section contours in `heart/root-necks-v129/` locate
the source mouth discrepancies; distance alone does not establish whether a
mouth is outside tissue or intrudes into atrial blood space.

The first exterior-only ownership trials (v131–v135) exposed slivers where an
existing atrial bore and the native inner vein course nearly coincide. They
are rejected, even where a numerical cell pass alone occurred: the high genus
and fragments do not describe a coherent root. The v136 complete-root prototype
instead retains the continuous finite root wall, trims it at the actual left
blood-space boundary, and assigns the small shared atrial/root tissue to the
root. All four complete roots have closed boundaries and zero self-crossings.
Before installation, the corresponding volume must be removed reciprocally
from the original atrial coats/myocardium, with common cut charts. This changes
local tissue ownership only; it does not move the heart, chamber boundaries,
vessels or skeleton, or add occupied volume inside the blood space. Complete
assembly tests and rendered sections are still required.

The first reciprocal re-import trials fail at the nearly coincident atrial
lining. They are not installed. The next ownership candidate keeps a 5 µm band
of the existing atrial lining behind the blood-space boundary: only root
controls derived from the diagnostic fluid cut move inward into the root
tissue. The original atrial cells and blood surface stay fixed; the band remains
atrial tissue, so it is not an anatomical gap. Five micrometres is an explicit
illustrative material-partition inset within the existing 40 µm illustrative
lining, not an audit tolerance or a measured histological thickness. Full root
and emitted shared-wall checks remain mandatory.

Root/host source faces also coincide at the exterior. A further candidate uses
an explicit positive root-only scale of 0.99995 about each root's own centroid
after the lining allocation. The maximum control displacement and actual
blood-space overlap must be measured before use. This is a bounded geometric
source repair, not permission to increase audit tolerances; all original heart
wall controls remain fixed. A failed candidate must not replace the installation.

After exact source-feature replay and shared storage repair, eight crossings
remain in two contact neighborhoods. The next local repair may adjust retained
heart-source controls in those neighborhoods by at most 10 nm; all other
retained source controls stay literal. New intersection controls retain the
20 nm total bound from the initial encoded candidate. This explicitly replaces
the earlier blanket carrier-control freeze only at these diagnosed contacts.
Skeleton, lungs, diaphragm, body placement and camera remain fixed. Each move
still needs the complete global no-new-crossing guard and final independent
cell, interface and occupied-volume checks.

The final three reversed fragments share original atrial source facet 4002.
Grid, grouped and link-preserving contractions did not produce an acceptable
export. A new source candidate translates only the right superior pulmonary
vein root 0.25 µm inward along that facet normal, clearing the near-corner cut.
This is an explicit source-placement repair; the other roots and complete
original heart wall remain fixed. Actual blood-space overlap, source topology,
shared construction and emitted geometry must all be checked again. The earlier
three-defect candidates are rejected and must not be installed.

For the final source-parent orientation repair, the explicit geometric bounds
may be widened locally to 250 nm on derived cut controls and 50 nm on retained
source controls in diagnosed contact neighborhoods. These are actual shared
coordinate edits, not enlarged intersection or acceptance tolerances. Every
move remains subject to the complete no-new-crossing guard, source-parent
orientation, closed-cell topology and subsequent occupied-volume checks.
All controls outside the diagnosed contact neighborhoods, other organs,
skeleton and camera stay fixed. Neither the failed coordinate-reuse trial
with doubled offsets nor any incompletely checked candidate may be installed.

## Checked pulmonary-vein root assembly — 21 October 2026

The emitted `pv-root-bank-v196` uses 11,979,849 bytes for the ten original wall
cells and four finite pulmonary-vein roots. All fourteen cells are closed,
have outward source-parent orientations and zero detected self-crossings.
All 63 occupied-volume pairs whose bounds overlap pass; the remaining pairs
have disjoint bounds. Fifty-two reciprocal tissue interfaces use common
coordinates and opposite normal signs. Eight complete source-centerline paths,
including every straight segment, reach the left atrial blood space without
crossing the actual emitted tissue. This checks those paths, not detailed
valve motion, every fluid streamline or distal joins to recipient lung vessels.

The accepted storage repair keeps all original carrier/root controls literal.
Derived intersection controls move at most 59.663 nm from their first encoded
positions. The separately declared 0.25 micrometre right-superior-root source
translation remains recorded in `pv-root-source-clearance-v177`.
Orientation solving and the global crossing guard both pass in
`pv-root-final-contact-v195`; failed intermediate candidates are not installed.
Actual emitted checks and path evidence are in `pv-root-bank-v196`.
No new texture, separate duplicate organ or additional TA98 coverage claim is
introduced by this root repair. The exact A12.1 coverage remains 31/113 terms.

Both revisions are installed locally, with exterior and coronal-section review
bound to their accepted jobs. Male catalogue SHA-256:
`8c4b38eadb35fbbedf17b1bd4880a0ff010f95d409f379e42c79b8b6562a6360`.
Female catalogue SHA-256:
`c8d59a50722414e95f41e0afed52418a8661b967070fc9a794ad9c20b3140cd8`.
The female emitted geometry independently passes all fourteen cell audits,
61 relevant occupied-volume pairs and the same 52 tissue interfaces. Its
positive uniform instance reuses the source bank and the checked fluid paths;
it adds no second geometry bank. Skeleton, lungs and diaphragm are unchanged.

## Native chordal partition plan — 21 October 2026

The current five atrioventricular cusp records contain visible, branching
chordal strands as well as the leaflet sails. This is existing source tissue,
not a missing shape that needs a synthetic replacement. Independent complete
audits in `apparatus-source-v203` find all five closed and without detected
self-crossings. Their genera (7–12) reflect the visible source chordal webs;
they are not evidence of extra chambers. The separate right pulmonary
semilunar cusp has an inherited closure defect and is a separate repair item.

The complete 113-term inventory above still applies. The planned additional
collective control is A12.1.00.023, chordae tendineae, below Heart; its source
subsets retain their left/right valve ownership. False chordae, annuli,
trigones, nodules and commissural features are not established by this source
finding. Papillary muscle names and existing partial qualifications stay intact.

First measure intrinsic neighbourhood shape and inspect the whole leaflet and
its strands. Then establish each actual sail-to-strand neck before splitting.
Every existing exterior source triangle region, source normal and UV must be
preserved exactly once. The shared neck cross-sections close actual solid
tissue with reciprocal coordinates and normal signs; they must not span a
valve orifice. No invented chordal routes, moved skeleton or overlapping
duplicate valve solids are permitted. Measured papillary contacts must be
distinguished from merely adjacent, unjoined source ends. Classification alone
does not establish an anatomical boundary or accepted representation.

References consulted: Gunnal et al., [116-heart cadaveric chordal morphology
study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4379635/) (indexed text;
the full PMC page returned a browser challenge); Mangine et al.,
[chordal construction study](https://arxiv.org/html/2411.09599v1), Methods 2.1.
The latter distinguishes synthetic branching from specimen-derived geometry.
Its figure files did not load through the research browser, so they have not
been visually accepted as reference evidence. The actual installed source
apparatus was inspected together in the local 3D atlas.

The first partition candidates use measured sections toward each leaflet sail,
with 8/10/8/9/8 separate small solid strand contours at the five proposed necks
(`chordal-sections-v207`). The planes are computational ownership boundaries
within the continuous source attachments; their histological locations remain
provisional. No source strand or leaflet shape is replaced. Thirteen inherited
nonoutward normal corners in the three native cusp meshes may be corrected to
their original face normals; all other existing source words and UVs remain
literal. New cut-corner fields are interpolated from their source faces.

## Installed native chordal controls — 21 October 2026

All ten emitted leaflet/chordal cells are closed, with zero detected
self-crossings. Each source pair has disjoint tissue volume and one literal
reciprocal neck interface. Independent exported-bank replay checks every
source control, original exterior normal corner, source-face area coverage
and UV interpolation. Only the thirteen documented nonoutward source normal
corners change. The new Float32 cut controls round by at most 59.602 nm in
the male source frame. The complete five-pair bank is 951,401 bytes and is
shared by both models; no new image assets are used.

The ordinary TA98 and Systems controls expose the collective chordae term
A12.1.00.023 and its five left/right source subsets. Peeling includes those
parts exactly once. The recovered partial LV papillary head is grouped under
the left ventricle; its explicit four-source registration qualifies adjacency
in the native male frame without certifying an anatomical attachment.

At this v218/v220 stage the papillary ends remained adjacent, unjoined source tissue. Annuli, finer
commissural anatomy and the other listed missing features are not supplied
by this split. A measured strand-neck ownership plane is not a histological
boundary. The then-unrepaired right pulmonary cusp and papillary ends are
superseded by the installed repairs described below.

Male catalogue SHA-256:
`0fd2945366c9dc9a4fcadaf3b623b498a416f2a11c294116a9608e07a1484244`.
Female catalogue SHA-256:
`93288408841b8d7b91af4df4ab41d6ddc659dbefde39aa9c3776e95a2cb1cd0d`.
Actual assembled apparatus and coronal sections were inspected in the local
atlas before installation. Emitted geometry, reciprocal interfaces, hierarchy,
peeling, camera intent and explosion checks pass. A first replay run exceeded
its 25-second deadline because of a repeated lookup; the corrected replay
passes in under one second at the same limit. Failed runs are not acceptance
evidence and were not installed.

## Right pulmonary cusp repair plan — 21 October 2026

A12.1.02.011 is present but its native source has a single seven-edge surface
defect, diameter 0.8474 mm. This is a breach in the finite leaflet tissue;
it is not the pulmonary valve orifice. Repair is limited to a five-triangle
patch on those existing boundary controls. Preserve every original face,
coordinate, UV and valid normal word. Candidate triangulations must pass the
complete source/candidate self-intersection test, closed genus-zero topology
and emitted-bank replay. Inspect the leaflet exterior, section and three-cusp
assembly before local installation. Valve/root placement, chamber walls,
skeleton, lungs and diaphragm remain fixed. Transfer the same literal bank
through the existing female leaflet placement only after male checks pass.

## Native chordal–papillary interface plan — 21 October 2026

Complete triangle contacts in `papillary-contacts-v228` find actual intersecting
surfaces at all five declared source insertions (185/331/313/184/169 crossing
face pairs). Prior point-to-vertex distances did not establish physical gaps.
These source tissues overlap without a conforming interface. Preserve their
union and overall course: assign the overlap to papillary tissue, trim the
chordal end reciprocally, and subdivide the original papillary surface into
one shared contact chart and its remaining free surface. Do not bridge an
invented gap or move papillary heads. The existing leaflet–chordal necks stay
shared and unchanged. Before installation, check complete cell topology,
self-crossings, occupied-volume ownership, source ancestry, Float32 storage,
actual exteriors and sections. Existing registration limits remain explicit.

The first complete fifteen-part candidate additionally exposed five inherited
tricuspid overlaps: posterior/anterior chordae, posterior/septal leaflet,
posterior leaflet/anterior chordae, posterior/anterior leaflet and the two
papillary muscle bases. That failed intermediate candidate was not installed. Resolve ownership
without moving the source union: anterior apparatus takes precedence over
posterior apparatus, then the septal apparatus over their remaining overlap.
Papillary basal boundaries remain illustrative regions in continuous muscle.
Leaflet-to-leaflet contacts are coaptation/sliding interfaces, not fused tissue.
They must not acquire a continuous-attachment relationship merely because
their geometric charts are reciprocal. Recheck the complete apparatus after
both ownership operations, rather than accepting the five isolated triplets.

## Installed cusp and complete apparatus repair — 21 October 2026

The right pulmonary cusp repair passed whole-source self-intersection,
closed genus-zero topology and actual three-cusp surface/section inspection,
then was installed in male v225 and female v227. Those revisions are retained
as provenance and are superseded by the current apparatus installation.

`apparatus-owned-bank-v242` passes the complete fifteen-cell Float32 audit,
all occupied-volume pairs with overlapping bounds and fifteen reciprocal
interfaces. The two explicit tricuspid ownership operations preserve the
source union and all previously declared interfaces. Ordinary atlas exterior
and coronal section views of both candidates were inspected before installing
male v244 and female v246. Female placement reuses the same literal bank with
the existing positive similarity; all fifteen recipient origins agree within
2e-15 m. Skeleton and unrelated anatomy remain fixed.

Current male catalogue SHA-256:
`c8b7bc5ce9d43a6214d722cb987c4bdb65376bfd91b942c6b0f6c176c2b4a64d`.
Current female catalogue SHA-256:
`b2ec6b89463cde5bc82304214eec9e57535f499226f8e5a806af4527c1f3399f`.
Apparatus bank SHA-256:
`94f6f5e087455fa970e08850a10031c8fab74f305a0904f3398a2154efd00a95`.

Independent local installation jobs passed in 20.32 s (male) and 20.41 s
(female). Installation and hierarchy/forest regeneration finished in 53.97 s
and 17.24 s respectively, within their recorded 90-second deadlines. Thirty-eight
installed-control, relation, shader, camera-intent and explosion checks pass;
TypeScript and the production build also pass. Numerical acceptance applies to
the checked assemblies, not to all cardiac anatomy or physiological valve
motion. Coverage remains 32 directly mapped A12.1 terms out of 113.

## Annular, basal and vessel continuation plan — 21 October 2026

The complete 113-term inventory above remains the construction inventory.
Current-frame full-triangle measurements (`wall-attachments-v248`) establish
crossing source surfaces at each of the five papillary bases against its
ventricular myocardium and endocardium. These are unresolved tissue ownership,
not missing bridges. Define the whole coupled assembly as the fourteen current
wall/pulmonary-root cells plus the fifteen leaflet/chordal/papillary cells.
Resolve their source overlap through reciprocal interfaces while preserving
their union, external course and every existing internal chart. The skeleton,
lungs, diaphragm and native vessel courses are fixed. Papillary-to-myocardial
contacts are continuous tissue; leaflet coaptation remains sliding. Contacts
with other chambers are not automatically accepted as attachments.

Measure the full leaflet-to-host contact charts to identify the real annular
hinges. A right/left fibrous ring is tissue (A12.1.00.027/.032), whereas an
atrioventricular or outflow opening is an orifice, never a solid disk. Actual
boundary curves and exterior surface patches may supply lightweight controls
for the corresponding TA98 boundary/opening terms. Do not count a generic
outer face as histology or a circular primitive as an anatomical annulus.

Anatomical references: [mitral apparatus imaging review, figures 2–4 and 9](https://pmc.ncbi.nlm.nih.gov/articles/PMC7712333/)
distinguishes the hinge line from the discontinuous posterior fibrous annulus
and shows trabecular papillary bases.
[Tricuspid annular histology study](https://pmc.ncbi.nlm.nih.gov/articles/PMC9318885/)
was available as indexed abstract text; its full PMC page returned a challenge.
These support preserving actual source attachment anatomy rather than making
thick circular rings. All new finite cells, full assembly contacts, stored
coordinates/normals, patent mouths and rendered sections must pass before male
installation. Only then transfer the same bank through the existing proper
female placement, review the recipient joins and install that revision.

## Native caval-root construction plan — 21 October 2026

The next adjoining source is `ZA:Superior vena cava` (A12.3.03.001), a native
open outer tube with thirteen twelve-control source sections. Its individual
outer sheet passes the crossing check, but its terminal section folds back
against the next section when a finite rim is made. The failed v263, v264 and
v266 finite trials were not installed. The v267 candidate reorients each source
section by a proper rigid rotation to its centreline tangent, preserving every
section centre, radius and area. A 0.6 mm corresponding inner sheet supplies
illustrative finite wall thickness; end closures are annular tissue only.
This is a finite root, not a new heart chamber or a plug in the blood space.
It passes its standalone closed-wall/crossing check in 0.21 s. Host attachment,
open communication, emitted storage and actual section inspection remain
required before acceptance. Resolve its tissue overlap with the complete
29-cell cardiac assembly using reciprocal source charts. Do not move the
heart, skeleton, diaphragm, lungs or hepatic dependency tree. The existing
thoracic IVC is owned by the hepatic communicating network and must not be
capped or replaced as an isolated closed object.

The ascending aorta (A12.2.03.001; 25 twelve-control source stations) and
pulmonary trunk (A12.2.01.001; six 23-control stations) are the next two
explicit adjoining roots. Retain their existing centre courses and sectional
areas; a source-section rotation may correct a folded station without reducing
the lumen. Construct corresponding 1 mm illustrative inner walls with annular
rims. Check each root, its correct ventricular outflow, both root-to-root
relationships, and its source leaflet contacts before coupling. Separate
semilunar leaflet contacts must not become fused coaptation. The great-artery
orifices are boundary annotations, never cap disks. This work is limited to
the heart roots and does not rebuild the downstream circulatory system.

The initial arterial trials v274/v275 failed complete self-crossing checks;
v277 separately failed the planar-section assumption. None was installed.
The aortic v276 sweep preserves axial areas to 1.9e-15 relative error; pulmonary
v278 preserves its nonplanar three-dimensional source contours by proper
rigid transforms rather than flattening the native collar. Both standalone
finite walls pass. Their full fluid envelopes nevertheless overlap, so the
pair is not accepted. Define two dependency groups, each root plus its three
semilunar cusps. Test a bounded common rigid separation of those groups before
coupling to the fixed myocardium. Do not repair an arterial crossing by
declaring a fused arterial interface. Semilunar source coaptation also has
occupied-volume overlaps; resolve explicit ownership without fusing the cusps.

The two-group v283 fit passes standalone geometry and mutual arterial
clearance, but the full coupled v285 diagnostic reveals a false aorta–SVC
contact. Both remain uninstalled. The revised fit must protect the fixed SVC
envelope explicitly, test all three vessel pairs, preserve the arterial groups
rigidly, and keep their individual translations within 8 mm in the transverse
plane. A successful arithmetic Boolean is not acceptance of a false vessel join.

## Installed coupled and caval stage — 21 October 2026

Further annular work uses the user's local `Male_2024.glb` only as a shape and
placement reference. Its named structure metadata identifies right/left
fibrous rings (cid 22809/18598), right/left trigones (20644/23777), auricles
(20454/23809), septomarginal trabecula (19988), fossa/limbus (22621/21663) and
pectinate muscles (18054). The original complete 113-row authority inventory
above remains the scope; these references do not count as implemented items.
Read the explicitly zero-weight default morph only. Do not copy its mesh or
texture into the atlas. Fit its landmarks to the actual native valve apparatus
before using them to constrain newly constructed regions.

Planned annular regions are finite connective tissue owned within the current
myocardial host, bounded by measured leaflet hinges and informed by the local
reference. They replace that local host volume, with reciprocal interfaces;
they must not be overlapping torus solids. The right/left orifices and arterial
mouths are separate boundary controls, with no filled disk. Source pectinate,
trabecular and conducting structures remain missing unless their actual
relief/tissue is constructed and checked, rather than represented by an
unmodified smooth chamber surface. Skeleton and the existing heart pose are
fixed. A passing male assembly precedes literal positive-similarity transfer
and independent recipient-interface checks in the female model.

Male v271 and female v273 retain the same 40 directly mapped A12.1 terms;
73 remain without a distinct exact control. The 30-cell bank has 108 reciprocal
charts, zero detected stored-cell self-crossings and passing complete occupied-
volume checks. All five native papillary bases and the five AV source leaflets
have shared host interfaces. This establishes macroscopic source attachments,
not complete fibrous annuli or clinical fine anatomy. The SVC has finite wall
attachment to right atrial tissue. Its centre path and eight 3 mm offset tracks
reach the right atrial blood space without intersecting the actual exported
tissue; this is a bounded passage check, not valve-flow simulation.

New SVC source sections were rotated at most 0.869373 mm in displacement;
centres and section areas stay fixed. Both recipients use the same literal
bank with their existing proper placements. Actual exterior and section views
were inspected. Numerical review/installation records and rollback are in
`svc-male-install-v271`, `svc-female-install-v273` and the matching named
`.local-models/ta98-review/cardiac-*-caval-*` manifests.

The first coupled installation job v258 reached its 100 s deadline and was not
installed. Exact geometry/attributes/native-version-bound audit replay then
allowed bounded replacement jobs v260/v262 to pass in about 17 s. Partial
checks were never treated as a passing candidate. The later full SVC bank
checks completed in saved 40.69 s, 40.49 s and 3.55 s batches; the complete
30-cell combination then passed in 3.27 s.

## Arterial continuation under validation

The uninstalled v388–v402 candidate constructs communicating aortic/LV and
pulmonary/RV entries while protecting both arterial roots and all six native
semilunar cusps. A 10 µm computational cutter clearance enters the protected
root wall only; the physical root wall and lumen do not change. Eighteen centre
and 3 mm offset tracks pass against the exported tissues with static closed
semilunar cusps excluded. This does not establish valve motion or distal joins.
The complete stored-cell checks still found small face crossings; the candidate
is not installed. The female native aortic arch contains five mouths and three
branches and must not be replaced by the male short ascending-root part. Two
zero-area source faces can be omitted to recover its closed, genus-four wall,
but a branch-preserving distal join remains required.

## Female thoracic packing correction

The user reports anterior cardiac tissue beyond the sternum. Check the actual
installed epicardial surface against the native female sternum, costal cage,
lungs, diaphragm and oesophagus. Keep the skeleton fixed. A common positive
similarity must carry every cardiac wall, leaflet, chordal, papillary and owned
root member; use the complete ownership closure. Preserve all shared bank data
and exact internal interfaces. The native complete aortic arch and independently
owned recipient vascular trees remain fixed and require separate distal joins.
Measure sternum-projected forward excess and complete 3D contacts before and
after any candidate; scaling from landmark centroids alone is insufficient.

### Installed female thoracic correction v426

The corrected common placement carries 40 owned heart parts and all of their
shared internal interfaces. The native full aortic arch and every nonmember,
including the skeleton, lungs, diaphragm and oesophagus, remain literal. No
geometry bank, UV or normal field was regenerated. The dependency forest now
owns the coordinated cardiac assembly explicitly instead of treating its
separate source objects as independent fit targets.

Before correction, the epicardial surface projected up to 24.214275 mm forward
of the native sternum posterior surface. After correction it is behind that
surface over every checked overlapping triangle projection; the independently
measured minimum 3D sternum gap is 0.947087 mm. Native chest-bone crossing pairs
fall from 2,538 to zero, diaphragmatic pairs from 1,514 to zero, and oesophageal
pairs from 697 to 583. Pulmonary pairs change from 2,632 to 2,629; that packing
remains unfinished and the counts do not establish anatomical acceptance.
The 0.9-scale candidate was rejected because it increased pulmonary crossings.

All 40 changed parts and 118 declared interfaces/clearances passed the guarded
22.01-second installation job. Actual exterior and sagittal section views were
inspected. Local installation completed in 62.30 seconds, within its recorded
90-second deadline. Rollback is `.local-models/ta98-review/cardiac-female-thoracic-v426/manifest.json`;
evidence is `.local-models/geometry-workbench/heart/female-packing-install-v426/` plus the four v421–v424 context
reports. Exact term-control coverage remains 42 of 113, with 71 still missing.


## Cardiac label expansion plan — 21 October 2026

The user's minimum is 100 distinct selectable TA98 labels in each main model.
The existing 113-term inventory remains authoritative. This expansion keeps
three separate counts: stored tissue controls, source-bound surface/opening
controls, and provisional point/curve reference representations. A point or
surface annotation does not establish microscopic tissue, a new wall layer,
or a patent opening. Do not count reference representations as finished meshes.

Use local Male_2024 named anatomy only to measure centroids and landmark
relationships, fitted to the installed valve/papillary apparatus. No reference
faces, UVs or texture pixels become runtime assets. Keep both complete cardiac
placements, every vessel, both lungs and the skeleton fixed. New controls reuse
existing carriers and the standard selectable labels. Openings label an empty
space or rim and add no solid disk. Fine structures use explicit reference
points or curves, with their representation and limits retained in the details.
Missing/variable structures remain excluded rather than assigned arbitrary
coincident labels. All 113 rows below retain their tissue implementation status;
label coverage will be reported separately from tissue completion.

## Representation correction and thirteen remaining terms — 21 October 2026

The v610 sidecars provide 100 selectable names, not 100 completed anatomical
representations. Their 58 added controls are provisional label locations on
parent tissues. They must not enter the completed geometry count. A surface
requires an actual surface domain, a boundary an actual curve, an opening an
unobstructed rim/channel, and a material structure its own finite tissue or
explicitly owned region. The full 113-row inventory above remains the baseline.
No new reconstruction in this section has yet been accepted or installed.

The installed geometric sources are male v594 and female v603, with the
obsolete female native pulmonary source hidden by v615. The complete cardiac
similarities remain fixed; all skeletal, muscular and pulmonary context is
fixed. Existing shared coordinates/normals, lumens and tissue union must survive
any partition. Male construction precedes the same proper similarity transfer
to female, followed by separate female packing checks. Reference-model faces,
UVs and texture pixels are not copied.

| TA98 | Required representation and ownership | Adjoining interface / condition |
|---|---|---|
| A12.1.00.008 | Apical notch: actual cleft/boundary in the ventricular exterior | Between ventricular apical regions; a point on an unnotched surface is insufficient |
| A12.1.00.024 | False chordae: finite fibromuscular band | LV septum to papillary muscle or free wall; no valve-leaflet insertion; variable anatomy |
| A12.1.00.028 | Infundibular tendon: owned finite fibrous tissue | Between arterial outflow support structures; preserve both outflow lumens |
| A12.1.00.029 | Tendon of Todaro: owned finite fibrous band | Inferior caval valve region to central fibrous support; not an exposed free string |
| A12.1.00.030 | Koch triangle: bounded endocardial surface domain | Tendon of Todaro, septal tricuspid hinge and coronary-sinus opening; TA98's sinu-atrial wording is a documented terminology error |
| A12.1.01.004 | Smallest cardiac-vein openings: real small channels/rims | Endocardial mouths with an actual empty lumen; no disks or dots counted as holes |
| A12.1.01.007 | Foramen ovale: conditional opening/closed adult site | No patent fetal communication added to the default adult assembly; an open variant needs a separate compatible septal state |
| A12.1.01.014 | Intervenous tubercle: superior interatrial fold/rim domain | Superior to fossa ovalis; do not invent an independent nodule unsupported by human dissection |
| A12.1.04.006 | Commissural cusps: owned finite mitral commissural leaflet regions | Anterolateral and posteromedial commissures; preserve leaflet union and chordal insertions |
| A12.1.08.004 | Bronchopericardial membrane: finite collagenous sheet | Carina/main bronchi, dorsal pericardium, pulmonary venous adventitia and diaphragm |
| A12.1.08.008 | Serosa: finite serosal tissue coat | Shared interface with subserosa; physiological cavity remains empty |
| A12.1.08.009 | Subserosa: finite connective tissue layer | Shared interfaces with serosa and the underlying existing tissue |
| A12.1.08.010 | Fold of left vena cava: finite serosal fold | Left atrial appendage/left pulmonary venous ridge toward left pulmonary artery |

Reference basis: [FIPAT TA98 heart](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.1.00%20TA98%20EN.htm),
[FIPAT pericardial definitions](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/12.1.08%20TA98%20EN.htm),
[FIPAT corrected AV-node triangle identity](https://ifaa.unifr.ch/Public/TNAEntryPage/auto/unit/LAEN/TAH3643%20Unit%20EN.htm),
[Loukas et al., 100-heart intervenous-fold dissection](https://pubmed.ncbi.nlm.nih.gov/22109383/),
[false-tendon autopsy study](https://pubmed.ncbi.nlm.nih.gov/3711514/),
[intact-pericardium dissections](https://pmc.ncbi.nlm.nih.gov/articles/PMC8655319/).
These references establish anatomical type and relationships, not exact
individual coordinates, thicknesses or source-model acceptance.

## Installed term controls — 21 October 2026, v640

Both ordinary local TA98 models now load 113 distinct cardiac/pericardial term controls. Heart details and Show components expose 112 descendants (the whole heart is the 113th term). Systems shows 113; canonical TA98 Heart contains 101 entries, with 12 further terms under pericardium and related spaces. Counts describe terminology controls, not new tissue meshes.

The original 42 geometry-backed term controls remain unchanged. All 71 additional references are explicitly provisional; object, surface, opening and microscopic geometry remain unfinished where stated. This metadata install does not include the failed v639 serosal prototype. A closed adult foramen-ovale reference creates no patent hole.

Evidence and both screenshots: `.local-models/geometry-workbench/heart/cardiac-all-terms-v640/`. Recoverable sidecar backups: `.local-models/ta98-review/cardiac-all-terms-v640/manifest.json`. Geometry, materials, skeleton and cardiac poses are unchanged.

## Typed reference correction — 21 October 2026, typed-reference-plan-v641

The 71 provisional controls below are audited before construction. The complete 113-term authority inventory above remains binding. No material structure is accepted as a point. Collections select their actual children; cavities and orifices remain empty. A default adult foramen-ovale site is closed. Source-space registration alone is not anatomical acceptance.

| TA98 | Structure | Authority parent | Required representation | Current carrier(s) |
|---|---|---|---|---|
| A12.1.00.006 | Right border | A12.1.00.001 | boundary-curve | TA98HEART:male:myocardium:right-atrium |
| A12.1.00.008 | Notch of cardiac apex | A12.1.00.007 | boundary-curve | TA98HEART:male:epicardium |
| A12.1.00.009 | Anterior interventricular sulcus | A12.1.00.001 | boundary-curve | TA98HEART:male:epicardium |
| A12.1.00.010 | Posterior interventricular sulcus | A12.1.00.001 | boundary-curve | TA98HEART:male:epicardium |
| A12.1.00.011 | Coronary sulcus | A12.1.00.001 | boundary-curve | TA98HEART:male:epicardium |
| A12.1.00.018 | Auricle | A12.1.00.017 | collection | TA98HEART:male:myocardium:right-atrium, TA98HEART:male:myocardium:left-atrium |
| A12.1.00.020 | Trabeculae carneae | A12.1.00.001 | collection | TA98HEART:male:endocardium, TA98HEART:male:myocardium:right-ventricle, TA98HEART:male:myocardium:left-ventricle |
| A12.1.00.021 | Vortex of heart | A12.1.00.001 | finite-mesh | TA98HEART:male:myocardium:left-ventricle |
| A12.1.00.022 | Papillary muscles | A12.1.00.001 | collection | ZA:Anterior papillary muscle of right ventricle, ZA:Inferior papillary muscle of left ventricle, ZA:Inferior papillary muscle of right ventricle, ZA:Septal papillary muscle of right ventricle, TA98CARDIAC:FJ2418 |
| A12.1.00.024 | False chordae tendineae | A12.1.00.001 | finite-mesh | TA98HEART:male:myocardium:left-ventricle |
| A12.1.00.025 | Right fibrous trigone | A12.1.00.001 | finite-mesh | TA98HEART:male:fibrous:right-annulus |
| A12.1.00.026 | Left fibrous trigone | A12.1.00.001 | finite-mesh | TA98HEART:male:fibrous:left-annulus |
| A12.1.00.028 | Tendon of infundibulum | A12.1.00.001 | finite-mesh | TA98HEART:male:myocardium:right-ventricle |
| A12.1.00.029 | Tendon of valve of inferior vena cava | A12.1.00.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.00.030 | Triangle of atrioventricular node (Koch) | A12.1.00.001 | surface-region | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.002 | Right auricle | A12.1.01.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.003 | Crista terminalis | A12.1.01.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.004 | Openings of smallest cardiac veins | A12.1.01.001 | opening-rim | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.005 | Fossa ovalis | A12.1.01.001 | surface-region | TA98HEART:male:myocardium:interatrial-septum |
| A12.1.01.006 | Limbus fossae ovalis | A12.1.01.001 | boundary-curve | TA98HEART:male:myocardium:interatrial-septum |
| A12.1.01.007 | Foramen ovale | A12.1.01.001 | surface-region | TA98HEART:male:myocardium:interatrial-septum |
| A12.1.01.008 | Pectinate muscles of right atrium | A12.1.01.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.009 | Opening of coronary sinus | A12.1.01.001 | opening-rim | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.010 | Opening of inferior vena cava | A12.1.01.001 | opening-rim | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.011 | Opening of superior vena cava | A12.1.01.001 | opening-rim | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.012 | Sinus of venae cavae | A12.1.01.001 | surface-region | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.013 | Sulcus terminalis cordis | A12.1.01.001 | boundary-curve | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.014 | Intervenous tubercle | A12.1.01.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.015 | Valve of inferior vena cava | A12.1.01.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.01.016 | Valve of coronary sinus | A12.1.01.001 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.02.002 | Right atrioventricular orifice | A12.1.02.001 | opening-rim | TA98HEART:male:fibrous:right-annulus |
| A12.1.02.007 | Supraventricular crest | A12.1.02.001 | finite-mesh | TA98HEART:male:myocardium:right-ventricle |
| A12.1.02.008 | Conus arteriosus | A12.1.02.001 | finite-mesh | TA98HEART:male:myocardium:right-ventricle |
| A12.1.02.009 | Opening of pulmonary trunk | A12.1.02.001 | opening-rim | TA98HEART:male:myocardium:right-ventricle |
| A12.1.02.014 | Nodules of semilunar cusps of pulmonary valve | A12.1.02.010 | finite-mesh | ZA:Right semilunar leaflet of pulmonary valve |
| A12.1.02.015 | Lunules of semilunar cusps of pulmonary valve | A12.1.02.010 | surface-region | ZA:Right semilunar leaflet of pulmonary valve |
| A12.1.02.016 | Commissures of semilunar cusps of pulmonary valve | A12.1.02.010 | point-set | ZA:Right semilunar leaflet of pulmonary valve |
| A12.1.02.020 | Septomarginal trabecula | A12.1.02.001 | finite-mesh | TA98HEART:male:myocardium:right-ventricle |
| A12.1.02.021 | Trabeculae carneae of right ventricle | A12.1.02.001 | finite-mesh | TA98HEART:male:myocardium:right-ventricle |
| A12.1.03.002 | Left auricle | A12.1.03.001 | finite-mesh | TA98HEART:male:myocardium:left-atrium |
| A12.1.03.003 | Pectinate muscles of left atrium | A12.1.03.001 | finite-mesh | TA98HEART:male:myocardium:left-atrium |
| A12.1.03.004 | Openings of pulmonary veins | A12.1.03.001 | opening-rim | TA98HEART:male:myocardium:left-atrium |
| A12.1.03.005 | Valve of foramen ovale | A12.1.03.001 | finite-mesh | TA98HEART:male:myocardium:interatrial-septum |
| A12.1.04.002 | Left atrioventricular orifice | A12.1.04.001 | opening-rim | TA98HEART:male:fibrous:left-annulus |
| A12.1.04.006 | Commissural cusps | A12.1.04.003 | finite-mesh | TA98HEART:male:mitral-anterior-cusp |
| A12.1.04.007 | Anterior papillary muscle of left ventricle | A12.1.04.001 | finite-mesh | TA98CARDIAC:FJ2418 |
| A12.1.04.008 | Posterior papillary muscle of left ventricle | A12.1.04.001 | finite-mesh | ZA:Inferior papillary muscle of left ventricle |
| A12.1.04.009 | Aortic vestibule | A12.1.04.001 | space-region | TA98HEART:male:myocardium:left-ventricle |
| A12.1.04.010 | Aortic orifice | A12.1.04.001 | opening-rim | TA98HEART:male:myocardium:left-ventricle |
| A12.1.04.011 | Trabeculae carneae of left ventricle | A12.1.04.001 | finite-mesh | TA98HEART:male:myocardium:left-ventricle |
| A12.1.04.016 | Nodules of semilunar cusps of aortic valve | A12.1.04.012 | finite-mesh | ZA:Right coronary leaflet |
| A12.1.04.017 | Lunules of semilunar cusps of aortic valve | A12.1.04.012 | surface-region | ZA:Right coronary leaflet |
| A12.1.04.018 | Commissures of semilunar cusps of aortic valve | A12.1.04.012 | point-set | ZA:Right coronary leaflet |
| A12.1.06.002 | Conducting system of heart | A12.1.06.001 | collection | TA98HEART:male:myocardium:interventricular-muscular-septum, TA98HEART:male:myocardium:interatrial-septum, TA98HEART:male:myocardium:right-atrium, TA98HEART:male:endocardium |
| A12.1.06.003 | Sinu-atrial node | A12.1.06.002 | finite-mesh | TA98HEART:male:myocardium:right-atrium |
| A12.1.06.004 | Atrioventricular node | A12.1.06.002 | finite-mesh | TA98HEART:male:myocardium:interatrial-septum |
| A12.1.06.005 | Atrioventricular bundle | A12.1.06.002 | finite-mesh | TA98HEART:male:myocardium:interventricular-muscular-septum |
| A12.1.06.006 | Right bundle | A12.1.06.005 | finite-mesh | TA98HEART:male:myocardium:interventricular-muscular-septum |
| A12.1.06.007 | Left bundle | A12.1.06.005 | finite-mesh | TA98HEART:male:myocardium:interventricular-muscular-septum |
| A12.1.06.008 | Subendocardial branches | A12.1.06.004 | finite-mesh | TA98HEART:male:endocardium |
| A12.1.07.001 | Transverse pericardial sinus | A07.1.03.001 | space-region | TA98HEART:male:epicardium |
| A12.1.07.002 | Oblique pericardial sinus | A07.1.03.001 | space-region | TA98HEART:male:epicardium |
| A12.1.08.001 | Pericardium | A07.1.02.105 | collection | TA98HEART:male:epicardium |
| A12.1.08.002 | Fibrous pericardium | A12.1.08.001 | finite-mesh | TA98HEART:male:epicardium |
| A12.1.08.003 | Sternopericardial ligaments | A12.1.08.002 | finite-mesh | TA98HEART:male:epicardium |
| A12.1.08.004 | Bronchopericardial membrane | A12.1.08.001 | finite-mesh | TA98HEART:male:epicardium |
| A12.1.08.005 | Serous pericardium | A12.1.08.001 | collection | TA98HEART:male:epicardium |
| A12.1.08.006 | Parietal layer | A12.1.08.005 | finite-mesh | TA98HEART:male:epicardium |
| A12.1.08.008 | Serosa | A12.1.08.005 | finite-mesh | TA98HEART:male:epicardium |
| A12.1.08.009 | Subserosa | A12.1.08.005 | finite-mesh | TA98HEART:male:epicardium |
| A12.1.08.010 | Fold of left vena cava | A12.1.08.001 | finite-mesh | TA98HEART:male:epicardium |

Current skeleton, heart, valves, vessel roots and surrounding organs are fixed for this reference pass. New reference geometry requires rendered and section review; the rejected v639 serosal bank must not be installed.

## Installed typed references — cardiac-typed-install-v677

Both ordinary local TA98 models now load 113 selectable cardiac terms: the
unchanged 42 established geometry controls and 71 typed reference controls.
The latter comprise 38 finite-mesh references (36 new buffers and two existing
LV papillary meshes), six indexed endocardial/leaflet surface regions, seven
boundary curves, nine unfilled opening-rim sets, two three-point commissural
sets, three unfilled space envelopes and six collections of actual children.
The canonical TA98 Heart branch has 101 terms; the Systems Heart also includes
the 12 related pericardial terms. Whole-Heart Contains and component label
plans have 112 distinct descendants. On-screen callout placement remains
limited by viewport space; the inspector list contains all descendants.

The surface references use actual free endocardial or leaflet facets, not
one small disconnected myocardial calculation face. The fossa and adult
foramen-ovale site remain surfaces; no patent adult hole is inserted. Mouth
rims never create plugs. Space-envelope crop faces are displayed as wireframes,
not tissue. Finite reference volumes have their own depth buffer pass so rear
facets cannot paint over their front surfaces.

All 39 new finite buffers (36 material references and three empty-space
boundaries) passed complete Float32 topology and self-intersection checks in
`typed-final-buffers-v665/float32-audit.json`. Rejected/incomplete candidates
were excluded. Surface construction is recorded in
`typed-endocardial-surfaces-v674/typed-geometries.json`. The installed source
and rollback manifests are in `.local-models/ta98-review/cardiac-typed-install-v677/`;
the installed verification and complete typed plan are in
`.local-models/geometry-workbench/heart/cardiac-typed-install-v677/`.

Original male/female catalogue bytes, the shared cardiac bank, skeleton and
all existing body-part poses remain unchanged. These are typed reference
representations, not a claim that 71 additional histologically established
tissue partitions have been completed. Individual smallest-vein channels,
conduction histological extents/terminal density, variable false tendons and
pericardial reflection boundaries remain qualified. Female external
pericardial insertions are transferred references and remain provisional.
