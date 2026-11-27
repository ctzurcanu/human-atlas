# Urinary bladder: full TA98 construction inventory

Male is the construction master; transfer follows review of the complete male assembly.
The source tables below record the inventory before construction. The current local installation is recorded at the end; source correspondence alone does not establish a finite layer.

| TA98 | Parent | Term | Representation | Male source mappings | Female source mappings |
| --- | --- | --- | --- | --- | --- |
| A08.3.01.001 | A08.0.00.000 | urinary bladder | whole-organ | Missing | LOCAL:female:bladder:187 |
| A08.3.01.002 | A08.3.01.001 | apex of bladder | point-landmark | Missing | Missing |
| A08.3.01.003 | A08.3.01.001 | body of bladder | tissue-region | TA98PART:ZA:Urinary bladder:A08.3.01.003 | Missing |
| A08.3.01.004 | A08.3.01.001 | fundus of bladder | tissue-region | TA98PART:ZA:Urinary bladder:A08.3.01.004 | HRA:VH_F_fundus_of_urinary_bladder_base; HRA:VH_F_fundus_of_urinary_bladder_dome |
| A08.3.01.005 | A08.3.01.001 | neck of bladder | tissue-region | TA98PART:ZA:Urinary bladder:A08.3.01.005 | Missing |
| A08.3.01.006 | A08.3.01.001 | uvula of bladder | tissue-region | Missing | Missing |
| A08.3.01.007 | A08.3.01.001 | median umbilical ligament | tissue-region | Missing | Missing |
| A08.3.01.008 | A08.3.01.001 | serosa | wall-layer | Missing | Missing |
| A08.3.01.009 | A08.3.01.001 | subserosa | wall-layer | Missing | Missing |
| A08.3.01.010 | A08.3.01.001 | muscular layer | wall-layer | Missing | Missing |
| A08.3.01.011 | A08.3.01.010 | trigonal muscles | tissue-region | Missing | Missing |
| A08.3.01.012 | A08.3.01.011 | superficial trigone | wall-layer | Missing | Missing |
| A08.3.01.013 | A08.3.01.011 | deep trigone | wall-layer | Missing | Missing |
| A08.3.01.014 | A08.3.01.010 | detrusor | tissue-region | Missing | Missing |
| A08.3.01.015 | A08.3.01.014 | unstratified part | tissue-region | Missing | Missing |
| A08.3.01.016 | A08.3.01.014 | bladder neck part | tissue-region | Missing | Missing |
| A08.3.01.017 | A08.3.01.016 | external longitudinal layer | wall-layer | Missing | Missing |
| A08.3.01.018 | A08.3.01.016 | circular layer | wall-layer | Missing | Missing |
| A08.3.01.019 | A08.3.01.016 | internal longitudinal layer | wall-layer | Missing | Missing |
| A08.3.01.020 | A08.3.01.010 | pubovesicalis | tissue-region | Missing | Missing |
| A08.3.01.021F | A08.3.01.010 | vesicovaginalis | tissue-region | Not applicable | Missing |
| A08.3.01.021M | A08.3.01.010 | vesicoprostaticus | tissue-region | Missing | Not applicable |
| A08.3.01.022 | A08.3.01.001 | submucosa | wall-layer | Missing | Missing |
| A08.3.01.023 | A08.3.01.001 | mucosa | wall-layer | Missing | Missing |
| A08.3.01.024 | A08.3.01.004 | trigone of bladder | tissue-region | Missing | HRA:VH_F_trigone_of_urinary_bladder |
| A08.3.01.025 | A08.3.01.024 | interureteric crest | boundary-curve | Missing | Missing |
| A08.3.01.026 | A08.3.01.024 | ureteric orifice | opening | Missing | HRA:VH_F_ureteral_orifice_L; HRA:VH_F_ureteral_orifice_R |
| A08.3.01.027M | A08.3.01.024 | internal urethral orifice | opening | Missing | Not applicable |
| A08.3.01.028F | A08.3.01.024 | internal urethral orifice | opening | Not applicable | Missing |
| A09.2.03.001 | A08.3.01.028F | filling internal urethral orifice (female) | opening state | Not applicable | Missing |
| A09.2.03.002 | A08.3.01.028F | voiding internal urethral orifice (female) | opening state | Not applicable | Missing |
| A09.4.02.001 | A08.3.01.027M | filling internal urethral orifice (male) | opening state | Missing | Not applicable |
| A09.4.02.002 | A08.3.01.027M | emptying internal urethral orifice (male) | opening state | Missing | Not applicable |
| A04.5.03.017 | A04.5.03.016 | rectovesicalis (bladder context cross-reference) | tissue-region | Missing | Missing |

The apex is a point landmark and the interureteric crest is a ridge; neither is an independent filler solid. The trigone is a selected inner-wall region bounded by two ureteric ostia and the internal urethral ostium. Filling/voiding entries describe states of the same opening, not additional holes. The female state parent above is corrected against the primary female-urethra section; the original cached male parent is retained in correction provenance. Rectovesicalis is explicitly cross-referenced by the bladder section while retaining its canonical pelvic visceral-muscle parent.

## Assembly and adjoining interfaces

One shared tissue bank must hold the whole wall, detrusor/trigonal and neck muscle subdivisions, connective coats and mucosa. Peritoneal serosa/subserosa covers the reviewed superior surface only; uncovered external tissue is adventitia. Mucosal folds taper out at the trigone. Named features must use distinct, anatomically located selections and ordinary labels.

Required interfaces: bilateral oblique intramural ureters with patent ostia; shared urethral neck; median umbilical ligament attachment; pubovesical and sex-specific muscle attachments; prostate or anterior vagina, rectum, pelvic floor, peritoneum, native vessels/nerves and fixed pelvic skeleton. Vessels and nerves may occupy the wall but cannot enter the urinary lumen. No skeleton motion, independent donor overlay, blind lumen disk or overlapping layer faces is allowed.

All new anatomical representations and interfaces remain unfinished until actual surfaces, sections, finite cells, source/normal replay, urinary patency and whole pelvic packing are verified. Source review meshes remain local.

Sources: [1](https://ifaa.unifr.ch/Public/EntryPage/TA98%20Tree/TA98%20EN/08.3.01%20TA98%20EN.htm), [2](https://www.histologyguide.org/slideview/MHS-214-bladder/16-slide-1.html), [3](https://pmc.ncbi.nlm.nih.gov/articles/PMC4566164/), [4](https://journals.biologists.com/dev/article/134/20/3763/64555/The-development-of-the-bladder-trigone-the-center)

Bound source inventories: `.local-models/geometry-workbench/bladder-completion/source-v1/male-inventory.json` and `female-inventory.json`.

## Male outlet interfaces before neck reconstruction

The adjoining male urethra and prostate are inventoried in full before constructing the neck interface. This step only builds the urinary wall/bore needed for a patent bladder outlet; it does not certify all urethral or prostatic anatomy.

| TA98 | Parent | Term | Representation | Current male source meshes |
| --- | --- | --- | --- | --- |
| A08.5.01.001M | A08.0.00.000 | male urethra | whole-organ | TA98PART:ZA:Urethra:A09.4.02.004; TA98PART:ZA:Urethra:A09.4.02.017; TA98PART:ZA:Urethra:A09.4.02.021 |
| A09.4.02.003 | A08.5.01.001M | intramural part | tissue-region | Missing |
| A09.4.02.004 | A08.5.01.001M | prostatic urethra | tissue-region | TA98PART:ZA:Urethra:A09.4.02.004 |
| A09.4.02.017 | A08.5.01.001M | intermediate part of urethra | tissue-region | TA98PART:ZA:Urethra:A09.4.02.017 |
| A09.4.02.021 | A08.5.01.001M | spongy urethra | tissue-region | TA98PART:ZA:Urethra:A09.4.02.021 |
| A09.4.02.030 | A08.5.01.001M | external urethral orifice | opening | Missing |
| A09.4.02.005 | A09.4.02.004 | proximal part | tissue-region | Missing |
| A09.4.02.006 | A09.4.02.004 | distal part | tissue-region | Missing |
| A09.4.02.007 | A09.4.02.004 | urethral crest | tissue-region | Missing |
| A09.4.02.008 | A09.4.02.004 | seminal colliculus | tissue-region | Missing |
| A09.4.02.009 | A09.4.02.004 | prostatic utricle | tissue-region | Missing |
| A09.4.02.010 | A09.4.02.004 | prostatic sinus | opening | Missing |
| A09.4.02.011 | A09.4.02.004 | muscular layer | wall-layer | Missing |
| A09.4.02.015 | A09.4.02.004 | mucosa | wall-layer | Missing |
| A09.4.02.018 | A09.4.02.017 | muscular layer | wall-layer | Missing |
| A09.4.02.020 | A09.4.02.017 | mucosa | wall-layer | Missing |
| A09.4.02.022 | A09.4.02.021 | navicular fossa | opening | Missing |
| A09.4.02.023 | A09.4.02.021 | valve of navicular fossa | tissue-region | Missing |
| A09.4.02.024 | A09.4.02.021 | urethral lacunae | opening | Missing |
| A09.4.02.025 | A09.4.02.021 | urethral glands | microscopic-or-duct-network | Missing |
| A09.4.02.026 | A09.4.02.021 | para-urethral ducts | microscopic-or-duct-network | Missing |
| A09.4.02.027 | A09.4.02.021 | muscular layer | wall-layer | Missing |
| A09.4.02.029 | A09.4.02.021 | mucosa | wall-layer | Missing |
| A09.4.02.012 | A09.4.02.011 | circular layer | wall-layer | Missing |
| A09.4.02.014 | A09.4.02.011 | longitudinal layer | wall-layer | Missing |
| A09.4.02.019 | A09.4.02.018 | longitudinal layer | wall-layer | Missing |
| A09.4.02.028 | A09.4.02.027 | longitudinal layer | wall-layer | Missing |
| A09.4.02.013 | A09.4.02.012 | internal urethral sphincter | tissue-region | Missing |
| A09.3.08.001 | A09.3.00.001 | prostate | whole-organ | Missing |
| A09.3.08.002 | A09.3.08.001 | base of prostate | surface-landmark | Missing |
| A09.3.08.003 | A09.3.08.001 | proximal part | tissue-region | Missing |
| A09.3.08.004 | A09.3.08.001 | peri-urethral gland zone | tissue-region | Missing |
| A09.3.08.005 | A09.3.08.001 | distal part | tissue-region | Missing |
| A09.3.08.006 | A09.3.08.001 | apex of prostate | surface-landmark | Missing |
| A09.3.08.007 | A09.3.08.001 | anterior surface | surface-landmark | Missing |
| A09.3.08.008 | A09.3.08.001 | posterior surface | surface-landmark | Missing |
| A09.3.08.009 | A09.3.08.001 | inferolateral surface | surface-landmark | Missing |
| A09.3.08.010 | A09.3.08.001 | lobes of prostate | tissue-region | Missing |
| A09.3.08.016 | A09.3.08.001 | isthmus of prostate | tissue-region | Missing |
| A09.3.08.017 | A09.3.08.001 | capsule of prostate | wall-layer | Missing |
| A09.3.08.018 | A09.3.08.001 | parenchyma | tissue-region | Missing |
| A09.3.08.019 | A09.3.08.001 | prostatic ducts | microscopic-or-duct-network | Missing |
| A09.3.08.020 | A09.3.08.001 | muscular tissue | tissue-region | Missing |
| A09.3.08.021 | A09.3.08.001 | trapezoid area | surface-landmark | Missing |
| A09.3.08.011 | A09.3.08.010 | inferoposterior lobule | tissue-region | Missing |
| A09.3.08.012 | A09.3.08.010 | inferolateral lobule | tissue-region | Missing |
| A09.3.08.013 | A09.3.08.010 | superomedial lobule | tissue-region | Missing |
| A09.3.08.014 | A09.3.08.010 | anteromedial lobule | tissue-region | Missing |
| A09.3.08.015 | A09.3.08.010 | middle lobe | tissue-region | Missing |

The open native urethral source supplies course and a cross-sectional guide; it does not establish finite mucosal/muscular/adventitial coats.

Retain the distal native course. Bridge both actual bladder-neck lips to an explicitly finite urethral wall below the neck, using paired annular strips only. Prostate bore construction follows a passing shared urinary volume, with source exterior UVs retained outside the bore.

Skeleton, kidneys and the reviewed renal-root charts stay fixed. The prostate outer shape and source UVs away from the shared urinary contact must remain unchanged. No independent solid may fill an opening. Ejaculatory ducts are retained and require their own verified outlets before complete urethral acceptance.

Source and full interface bindings: `.local-models/geometry-workbench/bladder-completion/male-outlet-sources-v1/inventory-and-assembly.json`.

## Adjoining perivesical nodes before packing correction

| TA98 | Parent | Term | Representation | Current source meshes |
| --- | --- | --- | --- | --- |
| A13.3.04.026 | A13.3.04.025 | postvesical nodes | whole composite of distinct node solids | ZA:Postvesical nodes |
| A13.3.04.027 | A13.3.04.024 | lateral vesical nodes | whole composite of distinct node solids | ZA:Lateral vesical nodes.l; ZA:Lateral vesical nodes.r |

Both entries are terminal in the cached authoritative descendant inventory. Preserve every node surface, normal and UV; permit only one translation per complete source dependency group. The source-bound plan and anatomical reference are in `male-perivesical-node-fit-v1/plan.json` in the local bladder workbench. The urinary wall, prostate and skeleton stay fixed.

## Initial male local installation

Revision `male-urinary-coats-v1` installs 22 bladder cells with 19 TA98 tissue/region controls, six ureter coats, nine adjoining urethral coats and the source prostate with a shared urinary bore. All 60 other kidney tissue surfaces and their normal fields remain literal; 620 renal-root facets are unchanged. The decoded wall has one connected fluid domain and no prostate occupancy in that domain. Relative coat depths, source contours and regional boundaries remain provisional.

Three retained perivesical node groups are translated as complete forest-owned objects. Their original triangles, normals and UVs remain unchanged. Whole swept-neighborhood surface tests, occupied-volume comparisons and actual Blender context views were checked; the skeleton stays fixed. Early outward-only and posterior-only trials were rejected. The installed poses use `male-perivesical-node-fit-v4`, not the earlier trials.

The ordinary local male model exposes selection, components, named tissue controls and seven wall-peel stages. Old coarse bladder and urethral sources remain suppressed. Different coordinate-bank origins now preserve shared-interface visibility correctly. Camera-intent, explosion-plane, explosion-viewport, TypeScript and targeted urinary/renal controls checks pass. The local server delivers the exact installed catalogue and all eight new banks. Blender surfaces/sections were inspected; direct browser-canvas inspection was unavailable because local browser control is disabled.

Still unfinished: apex/opening landmark controls, interureteric crest, uvula, median umbilical ligament, pubovesical and vesicoprostatic attachments, mucosal folds and refined surface/region review. Detailed urethral/prostatic anatomy and full-body urinary packing are incomplete. Female construction/transfer has not started. This installation is a local review revision, not anatomical acceptance or completion of the bladder task.

Bound installation, replay, server and rollback evidence: `.local-models/geometry-workbench/bladder-completion/male-proposal-v1/`. Initial revision catalogue SHA256: `4ddaff7497dd228c3e50651e26c382ab19ca56eeeb7ac056ec511757f9896cbb`.

## Trigonal feature construction

The construction embeds the interureteric crest and uvula in the existing shared tissue wall. Actual ostial curves come from the boundary of the inherited reservoir lining against the constructed patent collars; a computational coat ownership seam is not substituted for an anatomical orifice. Source coordinates are read from the fitted body, with its original connectivity, rather than the earlier unfitted carrier coordinates. All three ostial curves, the outer wall, renal roots and skeleton are fixed.

The crest follows the posterior lining between the two ureteric openings; the male uvula is a small median elevation above the internal urethral opening. Shape review uses [Carter's original bladder illustration](https://commons.wikimedia.org/wiki/File:Gray1140.png). A [114-patient MRI study](https://pubmed.ncbi.nlm.nih.gov/32449124/) reports a mean midpoint AP diameter of 2.4 mm for the crest. This diameter is not treated as a measured ridge height. Initial review inputs are 1 mm crest height, 2.4 mm transverse support width and 0.8 mm uvula height; these are illustrative construction parameters, not donor histometry.

Local conforming refinement is limited to the feature footprints. A single continuous, monotone height map repairs the complete source envelope before its finite tissue partitions are rebuilt. Applying that map independently to the old, very thin partition triangles caused folded faces and was rejected. The repaired source envelope passes the complete wall/outer/fluid checks; the anatomical openings remain open. It adds 2,335 source controls. The current construction has 34,526 positively oriented tetrahedra and preserves every repaired boundary facet. These tetrahedra are not runtime geometry.

Source-envelope surface and diagnostic section renders in the existing Blender window show the actual crest and uvula relief. Source review request: `dc5c334e9e425e28d695f96a80e50877c03014a9a8e0d9370f4df297312abf16`. The gross source contour and ostial surrounds remain coarse. The final finite-cell and normal review is recorded below. Source stage: `.local-models/geometry-workbench/bladder-completion/male-trigonal-envelope-v2/`; construction coordinates: `male-trigonal-coordinates-v5/`.

The mesher now uses a reversible binary coordinate conversion, with source controls replayed literally. Its actual oriented boundary and shared facets are checked before layer construction; uncertain determinant signs use exact binary rational arithmetic. A previous decimal conversion generated nearly flat/inverted construction elements after source restoration and was rejected. The current volume complex has no flat transitions and no inverted tetrahedra. A synthetic interface-flip regression verifies the reusable orientation checker without copying private anatomical source geometry into test fixtures.

## Trigonal male local revision (superseded by the ostial revision below)

Revision `male-urinary-trigonal-v2` is installed locally. It has 25 bladder tissue cells and 21 bladder TA98 controls, including the interureteric crest and uvula. They remain finite regions of the same mucosal wall. All seven peel stages include each cell once. The adjoining six ureter coats, nine urethral coats, prostate and 60 retained kidney cells form 101 stored cells across five banks. The renal-root facets, retained kidney surfaces/normal fields, 1,044 skeleton descriptors and three perivesical node descriptors are unchanged.

Actual encoded tissues pass individual geometry, occupied-volume and whole-wall/fluid-domain checks. Thirteen cross-chart interfaces have exact reciprocal triangles and opposed normal words; the 26 directional forest declarations now have those measured face parameters. Rendering uses barycentric normal fields from the unsplit source surface to avoid tissue-partition shading seams. Final surface, coronal and sagittal views were inspected in Blender request `00a2a2031474881df52da1caa63c95a82f216deca5c3a6672c198ded90513fca`; its diagnostic cuts are explicitly uncapped. The browser canvas was not inspected because browser control is disabled. The local server serves the exact catalogue and five raw/compressed banks. All 44 targeted controls, interface, camera-intent and explosion tests pass.

The bladder remains incomplete: apex/opening controls, median umbilical ligament, pubovesical/vesicoprostatic attachments, mucosal folds, refined contours, detailed urethral/prostatic anatomy and full-body urinary packing remain. Female transfer has not begun. This is local-preview acceptance only.

Installation, literal bank replay, interface measurements and server evidence are in `.local-models/geometry-workbench/bladder-completion/male-trigonal-proposal-v2/`. Rollback: `.local-models/ta98-review/male-urinary-trigonal-v2/manifest.json`. Current male catalogue SHA256: `0a776e88bf3f7149062300d010512b5452c2947abb618652b1dfa1450078979a`; the female catalogue is unchanged.

## Ostial controls before construction

The refreshed 33-term authoritative inventory is `source-v3/male-inventory.json` under the bladder workbench. The next partitions are finite portions of the existing mucosa surrounding the two inherited ureteric rim curves and the inherited internal urethral rim. Each opening remains empty; its named control selects the actual surrounding wall and is bound to its clinical curve, rather than the computational coat seam inside the collar. A 1.25 mm support width is an explicit selection footprint, not a claimed histological thickness. It does not alter the repaired source envelope, renal roots, neck lumen or skeleton. Crest/uvula tissue ownership retains precedence where their supports meet.

The source curves have 19/23 ureteric controls and 19 neck controls; their nonplanarity is retained. Principal source spans are approximately 4.31 × 3.72 mm left, 5.47 × 4.10 mm right and 5.95 × 5.31 mm at the neck. These coarse contours still require shape refinement; creating their controls will not establish final anatomical acceptance. Their three-point trigonal relationship and oblique intramural course are reviewed against [OpenStax's urinary transport anatomy](https://openstax.org/books/anatomy-and-physiology-2e/pages/25-2-gross-anatomy-of-urine-transport). The complete whole wall remains the construction source before partitioning its tissue.

The first unrestricted support field unnecessarily subdivided deeper coats and produced duplicate/degenerate faces in three raw cells; it was rejected before packing or installation. The revised support can be positive only where all existing bladder-mucosal ownership prerequisites are positive. A strict binary margin separates its plane from existing coat interfaces. This confines the additional construction to mucosa and preserves the deeper tissue surfaces.

Both radial supports (1.25 mm and 1 mm) were subsequently rejected by the encoded intersection gate. Their nonlinear fields produced nearly coplanar internal/exterior fragments near the right ostium. The current construction instead uses affine halfspaces about each unchanged clinical rim, with an eight-sided convex selection footprint, 1 mm lateral margin and 4 mm axial margin. These are selection extents over the existing tissue domain, not a solid spanning the opening. All incident tissues receive the same conforming subdivisions; only mucosa receives the opening ownership labels. The stored curve must replay on actual exposed lining faces before the named control is installed.

The regionalized affine trial also failed raw topology. The next construction uses Cartesian selection windows with outward-rounded binary extents (2^-20 m grid). Its field values, intersection coordinates and normals use the exact axis constraints, including on nearly flat construction elements. None of the 18 planes cuts any of the 620 fixed renal-root facets. This changes the selection footprint, not the clinical rim or anatomical opening. The physical support cells retain the mucosa identifier; the opening identifier belongs to the curve-witnessed concept, with no opening solid added. A synthetic near-flat-element regression checks exact plane incidence and closed cells.

## Current male local ostial revision

Revision `male-urinary-ostia-v3` is installed locally. Its 28 bladder tissue cells expose 23 bladder TA98 controls, including both ureteric orifices and the internal urethral orifice. Each opening is a boundary over actual exposed mucosa, not a solid or lumen disk. All source curve edges are sampled with actual face and barycentric witnesses; maximum source offsets are below 0.84 nanometres. The surrounding tissue retains the mucosa identity. Ordinary selection, labels, components and the existing seven peel stages include the new cells.

All 104 encoded urinary cells pass the individual geometry checks. The assembled wall is one closed genus-two tissue boundary and its fluid domain is one connected closed genus-zero boundary. All actual lining faces remain; four reference-only fragments have exactly repeated stored coordinates and no area. One folded internal axis-plane patch was retriangulated without changing its directed boundary, holes or coordinates. No clinical rim, kidney root or outer bladder surface was moved by that repair. Fourteen external interfaces have literal reciprocal triangles and opposed normal words.

The five shared banks occupy 22,561,412 raw bytes. No texture image was added. Independent replay preserves all 60 kidney surfaces, normal fields and UVs, 1,044 skeleton descriptors and the three installed perivesical node descriptors. Actual whole, coronal and sagittal Blender views were inspected in request `80830c68f8ef0706dd6a0dd705ed468deb019039a363fd30ba7076e2f02bafab`; sections are diagnostic uncapped cuts. Original Blender objects/transforms were preserved and its main file was not saved. Direct browser-canvas inspection remains unavailable. The local server supplies the exact catalogue and all five raw/compressed banks. All 45 targeted urinary/renal, shared-interface, camera-intent and explosion tests pass.

This remains an incremental local revision. Coarse ostial contours and visible lining/relief faceting still need refinement. Apex identity, median umbilical ligament, muscular attachments, folds and complete pelvic packing remain unfinished; detailed urethral/prostatic anatomy and female transfer are also incomplete. The revision does not establish anatomical completion.

Bound job, Blender review, actual geometry/attribute replay, interface parameters and serving evidence: `.local-models/geometry-workbench/bladder-completion/male-ostial-proposal-v1/`. Rollback: `.local-models/ta98-review/male-urinary-ostia-v3/manifest.json`. Current male catalogue SHA256: `fe692fb230657fa9670b97a1eed7a39e9a399b43e2c460b8babbc2a07ce4bccc`. The female catalogue remains unchanged at `ad8ebd03b3ba815cf3105e8ed11649a812c08dade8da58e5d405573d433abe77`.

## Apical attachment assembly before construction

The refreshed 33-term source inventory is `source-v4/male-inventory.json` under the bladder workbench. It distinguishes 23 existing controls from primary material identities; a composite control is not an additional tissue cell. Apex and median umbilical ligament remain absent. `male-apical-sources-v1/sources-and-assembly.json` binds the current catalogue, full ownership forest, actual bladder exterior, transversalis fascia, linea alba, both umbilical skin patches, pubic symphysis and prostate. These sources and the skeleton are fixed.

The planned median umbilical ligament is one finite fibromuscular cord. Its bladder contact uses an explicit connected patch of the existing anterior adventitial surface, with literal opposed cap triangles and unchanged source coordinates. The apex is a landmark on this attachment, not an invented apical filler solid. Its course follows the posterior side of the actual transversalis fascia towards the height of the native umbilicus. The upper contact references the fascia beneath the umbilical skin patches; this open source sheet does not establish a finite fascial tissue interface or an independently dissected deep umbilical ring. That limitation must survive construction and review. No bladder lining or urine passage is cut.

The course is checked against the published anatomical description and MRI caption in [Imaging of the Urachus](https://pubs.rsna.org/doi/10.1148/rg.240122): the remnant ascends between transversalis fascia and parietal peritoneum to the umbilicus and has a round or mildly oval section. The article text and caption were accessible; the requested image/PDF screenshots were unavailable, so no image-based acceptance is asserted. Section sizes are explicit illustrative inputs, not measurements of this donor. Whole-cell geometry, literal contacts, lumen exclusion, surrounding meshes and actual Blender views must pass before installation.

## Adjoining greater omentum before a local corridor refit

The apical cord has a closed, connected genus-zero surface and joins the unchanged urinary wall by 32 literal cap facets. The assembled boundary remains one closed genus-two surface; the original fluid geometry is unchanged and has no contact or containment overlap with the cord. Both source contact patches and all source normals replay literally. Blender whole/context/section request `e4089f97821a47fabc31e8580796299b0bef3d9d899e0b637cf262c97da31020` was inspected. This candidate is not installed: the complete 70-neighbour scan found greater-omental and short ileal wall crossings.

At 62 sampled sections the unchanged ellipse cannot translate into the measured fascia–omentum corridor with 0.1 mm clearance on each side; the worst deficit is 1.738 mm. This is a witness against that translation-only proposal, not proof that every possible course is impossible. Retain the prescribed section instead of shrinking it solely to pass. A denser fascia-following course is being tested, followed by a bounded posterior-only refit of the local greater omentum. The skeleton and bowel remain fixed. The explicit dependency closure contains the single whole omental source and no fixed blocker; its upper gastric/colic attachment region and lateral limits are pinned. The native omentum is an open sheet with 134 boundary edges, not a completed finite peritoneal wall. A refit must preserve that source topology and all UVs; it cannot establish the missing wall thickness or its descendants.

The complete 11-term authoritative inventory, source identities and dependency closure are bound in `.local-models/geometry-workbench/bladder-completion/male-omentum-corridor-v1/`. No descendant is silently represented by the coarse whole source.

| TA98 | Parent | Term | Class | Existing exact source mappings |
| --- | --- | --- | --- | --- |
| A10.1.02.201 | A10.1.02.002 | greater omentum | whole-organ | ZA:Greater omentum (explicit named open source; no separate child geometry inferred) |
| A10.1.02.202 | A10.1.02.201 | gastrophrenic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.203 | A10.1.02.201 | gastrosplenic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.204 | A10.1.02.201 | presplenic fold | tissue-region | Missing exact source mapping |
| A10.1.02.205 | A10.1.02.201 | gastrocolic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.206 | A10.1.02.201 | phrenicosplenic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.207 | A10.1.02.201 | splenorenal ligament | tissue-region | Missing exact source mapping |
| A10.1.02.208 | A10.1.02.201 | pancreaticosplenic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.209 | A10.1.02.201 | pancreaticocolic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.210 | A10.1.02.201 | splenocolic ligament | tissue-region | Missing exact source mapping |
| A10.1.02.211 | A10.1.02.201 | phrenicocolic ligament | tissue-region | Missing exact source mapping |

## Male apical attachment local revision v4

Revision `male-urinary-apical-v4` is installed locally, with 29 bladder/attachment cells and 24 bladder TA98 controls. The median umbilical ligament is a separately selectable finite cord. Its bladder attachment uses 32 literal adventitial facets with opposite winding and normal words. Shared caps disappear only when their declared neighbour is assembled in the same pose; separating a tissue restores its boundary. The upper contact is one literal native fascial facet with the same normal field. It is explicitly an open-reference-sheet contact, not a finite fascial wall or a reconstructed deep umbilical ring. The apex landmark itself remains unfinished.

The corrected fascia-following course clears all 69 unchanged neighbouring sources. A conforming local greater-omental refinement and bounded Z refit clear the cord without moving the bowel or skeleton. The maximum omental shift is 2.586 mm; the emitted Float32 surface introduces zero new crossing/coplanar-overlap pairs across the complete 111-source shortlist. Original crossed reference regions were pinned and remain anatomically unaccepted. The omentum retains its open source topology and UV field; normal words are transported through the actual refit. No missing omental child is claimed by this parent-sheet correction.

Actual whole, oblique and sagittal Blender views were inspected in request `a452eb93a6acb0cb1b81b45ecce81feef688ee4576ac399b98da1b5a78f3fbdf`; the section is diagnostic and uncapped. Original Blender objects/transforms were preserved and the main file was not saved. Independent replay of the installed banks preserves all 28 previous bladder cells and their attributes, matches the reviewed ligament/omentum arrays and preserves all 1,044 skeleton records. The ligament adds 4,322 triangles. Active raw-bank payload increases by 245,472 bytes, with no new image; this is file accounting, not a measured GPU-memory or performance claim.

The ordinary atlas now includes attachment selection, its standard label, components and an additional peel stage. The complete forest owns all 6,781 male parts; the ligament moves with the bladder, while fascia and omentum retain explicit separate refit obligations. All 46 targeted urinary/renal, interface, camera-intent and explosion tests pass on the installed revision, and TypeScript passes. The server returns the literal catalogue and all three new raw/gzip banks. Direct browser-canvas inspection remains unavailable. The female catalogue is unchanged.

This does not complete the bladder: apex identity, pubovesical/vesicoprostatic attachments, mucosal folds, contour/faceting refinement, complete urinary packing and female transfer remain. Bound numerical, render-replay, installation, registry, forest and serving evidence: `.local-models/geometry-workbench/bladder-completion/male-apical-proposal-v2/`. Rollback: `.local-models/ta98-review/male-urinary-apical-v4/manifest.json`. Current male catalogue SHA256: `a59a31582f5e34cfa888c9153a3e01a56205ca8443d96deb470a16411df44660`.

## Muscular extension assembly before construction

The refreshed complete 33-term inventory is `source-v5/male-inventory.json`. TA98 A08.3.01.020 is the paired **m. pubovesicalis**, and A08.3.01.021M is **m. vesicoprostaticus**. Both belong to the muscular layer; a support ligament cannot substitute for either term. A08.3.01.021F is the female vesicovaginal muscle and is excluded from male construction. Filling/voiding descendants remain states of an opening, not extra solid parts.

The posterior free-space loft is rejected: the supplied irregular contact patches produce self-crossing strips, and a ring approximation is unsuitable for this short attachment. The replacement construction domain is the existing finite posterior neck compartment, between its literal external-longitudinal-muscle and prostatic-base interfaces. Diagnostic partition `male-vesicoprostatic-partition-v1` finds one connected 1.894 mm³ volume with 17.514 mm² muscle contact and 5.439 mm² prostatic-base contact. It touches no other tissue peer. Its explicit construction window is X ±3 mm, Y 889.5–893.5 mm, Z posterior to −21 mm in the unchanged source frame. These are illustrative regional bounds, not measured histological dimensions or proof of anatomical completion. Before installation, partition that entire finite compartment with reciprocal interfaces, carry the induced subdivisions into the adjoining source cells without changing their physical surface, and inspect the actual assembled/sectioned result. The source lining, openings, prostate lumen and skeleton remain fixed; classification of this finite support region as muscle is a new construction requiring anatomical review, not a claim that native adventitia establishes muscle identity.

The proposed pubovesical extensions run from the anterior external longitudinal neck muscle towards the posterior inferior pubic bodies, on each side of the symphysis. The posterior vesicoprostatic extension joins the posterior neck to the prostatic base in the region between the seminal-gland/vas insertions. These relationships are checked against the accessible text of [Sklinda et al., Normal 3T MR Anatomy](https://onlinelibrary.wiley.com/doi/10.1155/2019/3040859) and the surgical anatomy in [Rassweiler et al., 2026](https://www.frontiersin.org/journals/surgery/articles/10.3389/fsurg.2026.1804051/full). Figure links were located but image retrieval failed; no inspection of those figures is claimed.

`male-muscular-attachment-sources-v1/sources.json` binds the current catalogue and literal muscle, adventitia, prostate, seminal glands, deferent ducts, pubis and pelvic skeleton. The source skeleton and ducts remain fixed. Construction uses explicit source contact facets and transverse elliptical section inputs; these are illustrative dimensions, not donor histometry. Each extension must share a true muscle attachment. Because the source muscle is covered by adventitia, the overlying adventitia must be repartitioned at the emerging fibres. An overlapping added muscle cannot be installed as a completed join. The entire existing urinary lining, ostial curves and renal roots remain fixed. Whole-cell, tissue-occupancy, source-interface, body-packing and actual surface/section review remain required.

## Current male muscle attachment local revision v5

`male-urinary-muscles-v5` is installed locally. The bladder now has 32 finite wall/attachment cells and 26 bladder TA98 controls; the shared urinary registry has 108 cells. The added parts are paired pubovesicalis [A08.3.01.020] and vesicoprostaticus [A08.3.01.021M]. Each is selectable in TA98 and included exactly once in muscular peeling and whole-bladder selection. Left/right names use the ordinary label protocol. No additional image texture is allocated. Active replaced geometry banks grow by 157,754 raw bytes (about 154 KiB); this is not a measured GPU-memory or frame-time result.

The paired anterior muscles retain their literal source-muscle and pubic contact facets, while the covering adventitia is partitioned around the emerging fibres. Their terminal bone contacts are explicitly typed `fixed-native-surface-contact`; they do not move the skeleton or claim a separately dissected bone insertion. The rejected posterior ring loft is replaced by a thin finite regional partition along the existing neck–prostatic-base compartment. Its construction window extends to Y 894.0 mm to avoid detached remnants at the earlier diagnostic boundary. Its dimensions and histological classification remain illustrative and require further anatomical review.

The emitted arrangement has 32 independently checked cells, all passing topology and self-intersection gates. Its 105 decoded reciprocal interfaces replay literal opposed triangles and normal words. The three new muscles pass all 112 comparisons in the complete active-body bounds shortlist, with no proper crossings or undeclared coplanar overlaps. Twenty-six unaffected source tissue surfaces remain literal, including the lining and clinical-opening support cells. All skeletal descriptors, renal roots and the female catalogue remain unchanged. One explicitly diagnosed 0.550 nm derived intersection edge was identified within its measured 2.039 nm combined encoding-error bound; no original control moved. Two coplanar internal diagonals were changed without moving their boundaries. Full emitted-cell and contact checks were rerun afterwards.

Actual Blender surfaces and uncapped diagnostic sections `9ada600884c2706e082d70fd8580e8ea07081d8f61418ff07765431ec45c180c` were inspected and replayed against the proposed banks, including the actual indexed attributes. The section views retain the urine passage. Forty-three installed urinary/renal, camera-intent and explosion tests pass; the type check passes. The served catalogue and all four raw/compressed banks match the installed files. Browser canvas, procedural runtime-material appearance and device performance were not inspected.

This is a reviewed local first version, not completed bladder anatomy. Apex identity, state-specific opening representations, mucosal folds, coarse contours/ostia, wider anatomical acceptance and female transfer remain unfinished. Source-bound proposal and installation evidence: `.local-models/geometry-workbench/bladder-completion/male-muscular-proposal-v1/`; numerical delivery: `male-muscular-delivery-v3/`; interface/body checks: `male-muscular-verification-v1/`. Current male catalogue SHA256: `f816b8c63d70f1118616a21be3147c53008a564694bca369efeb743207847e10`. Rollback: `.local-models/ta98-review/male-urinary-muscles-v5/manifest.json`.

## Planned lining refinement after v5

The authoritative inventory and owned regions above remain the construction scope. The current body mucosal sheet has actual edges up to 22.85 mm; a synthetic normal field on that mesh alone cannot supply a convincing finite fold profile. Refine only the free body/fundus lining where folds are supported, propagating every marked edge through all incident canonical tissue surfaces. Keep outer coverings, muscular attachments, prostatic interfaces, renal roots, clinical ostial curves and the skeleton fixed. No new hole, solid fold overlay, independently overlapping mucosa or microscopic-structure claim is permitted.

The proposed rugae are finite relief of the existing mucosal tissue, with a continuous shared-control displacement towards the actual fluid compartment. An explicit 0.35 mm maximum height and 4 mm nominal spacing are illustrative review inputs, not donor measurements. The trigonal, crest, uvular, ostial and neck supports receive zero displacement; a 4 mm smooth falloff avoids a hard boundary. The supporting source mucosa and lamina propria are folded together as one tissue region. The relative wall dimensions remain provisional. Inspect the resulting whole lining, junctions and longitudinal/transverse cuts, then check each encoded finite cell, reciprocal normal field and the complete fluid boundary before installation. Partial visual relief does not establish a completely empty or completely distended bladder state.

Reference: [University of Leeds bladder histology](https://histology.leeds.ac.uk/home/urinary/bladder/) describes folded mucosa accommodating volume change. [Hickling et al., Anatomy and Physiology of the Urinary Tract](https://pmc.ncbi.nlm.nih.gov/articles/PMC4566164/) supplies the smooth trigone and opening relationships. The existing local male/female bladder color and normal atlases were visually consulted for fold character only; their pixels and UV layouts are not copied.

The first 0.8 mm lining refinement exceeded the fixed 100,000-triangle ceiling and is rejected. A 1.2 mm refinement retains that ceiling and produces 0.35 mm actual relief. Its complete fluid cell remains connected, but deforming only exposed controls leaves 64 crossings against internal mucosal partitions. That version is also rejected and uninstalled. The subsequent radial field transports the complete finite mucosal/submucosal band, including regional interfaces. It decays to zero at the verified outer boundary of that combined band; outer coats and attachments remain fixed. The guide is the literal union of 14 source tissue cells and the existing fluid cavity, with zero volume discrepancy. Open3D supplies ray proposals, not a certificate of first-hit correctness or anatomical acceptance.

`male-rugae-v7` restores unchanged-face normal fields from the installed source rather than an earlier rejected deformation. All source-corner witnesses and packed outward-normal gates pass. Independent complete-fluid validation passes in 24.88 seconds: one closed genus-zero component, no detected self-crossings, all 14,877 source lining facets replaced and the original 9,255 nonlining facets retained. Computational terminal disks remain proof-only geometry. Actual emitted-normal coronal and sagittal CPU cuts were inspected; these are diffuse diagnostics, not Blender or atlas-material acceptance.

All 29 emitted wall-cell checks are retained individually: 26 pass; body mucosa has 79 crossing pairs, body submucosa 14 and fundus mucosa 10. The candidate is rejected and uninstalled. A measured adaptive continuation splits 191 source edges across every incident canonical surface (212 source triangles affected), with no global resolution increase. New controls require reevaluation of the same wall field; splitting the already deformed triangles alone would not repair a crossing. Current installed male SHA256 remains `f816b8c63d70f1118616a21be3147c53008a564694bca369efeb743207847e10`; female remains unchanged. Construction, checks and renders use separate hard time limits and retain elapsed-time evidence. The 120-second monolithic and 90-second batched checks stopped at their limits; the final remaining cell completed separately in 20.69 seconds without repeating the 28 saved checks.

## Installed male lining revision v6

The subsequent source-ray refinements did not consistently remove crossings. They remain rejected. A continuous physical-coordinate relief field now uses separate C2 gates for the protected clinical lining (3.5 mm transition) and the outer boundary of the finite mucosal/submucosal band (1.2 mm illustrative transition). The first combined 3.5 mm distance gate produced only 0.078 mm relief and was rejected after visual inspection as too shallow. The installed candidate reaches 0.346456 mm actual displacement while retaining the outer wall, trigone, ostial supports, neck and all external urinary contacts. These are explicit construction parameters, not donor-specific fold measurements or a verified filling state.

Measured adaptive source refinements reduced the stronger field's remaining crossing pairs from 28 to 8 to zero. The continuous field reference is immutable; previously verified point distances are reused only after source geometry, exact point coordinates and saved query hashes replay. For the 19-point refinement, clinical/covering distance stages completed in 0.81/1.12 seconds instead of recomputing more than 82,000 points. A nearly tangent Float64 averaged normal rounded to a tangent Int16 word at one corner; only that corner now uses the actual outward facet normal. The strict emitted-normal gate remains unchanged and a thin-triangle regression covers this rounding case.

The installed candidate is `male-rugae-v14`. All 29 emitted cells pass topology/self-crossing checks, all 406 occupied-volume comparisons pass, and the complete fluid boundary passes as one closed genus-zero component. All 18 external canonical interfaces preserve literal coordinates and normal words, including 51,275 original fixed controls. All 101 reciprocal interfaces in the proposed chart and its adjoining tissues replay with exact opposing stored normal words. The three existing clinical opening curves retain their points and bind to actual refined free-lining facets. Actual emitted-normal front/right surfaces and coronal/sagittal cuts were inspected using a CPU z-buffer renderer; diffuse diagnostic images do not establish browser-material or frame-time acceptance.

The guarded installation audit passed in 134.52 seconds, then the reversible local installation completed in 52.77 seconds. Its initial request included constraints between geometrically unchanged parts and was rejected; the corrected request binds only the 14 affected existing verified seams. No numerical gate was weakened. New bank cost is 17,729,537 raw bytes / 5,751,581 gzip bytes, with no new image assets. These are file sizes, not measured GPU allocation or rendering performance.

Installed male SHA256: `9c955abf768c4f6a4d31b725523c2e7e7b78e18c45922996c4e98b6d5fbbaa23`. Revision: `male-urinary-lining-v6`; installation ID: `7e4897cceb1b5fae9794d24789eaa16688bb12fc4eb51987fdd29d6bd02e56eb`. Rollback: `.local-models/ta98-review/male-urinary-lining-v6/manifest.json`. Evidence: `.local-models/geometry-workbench/bladder-completion/male-rugae-v14/local-proposal-v2/`. The current registry imports all 6,784 parts and the complete ownership forest retains 888 fixed skeletal items. Thirty-two current urinary/pubic shared-interface declarations are independently checked and bound; unrelated provisional source declarations remain unaccepted. The bladder attachment closure contains 36 items and no fixed skeleton member. All 43 urinary/renal, camera-intent and explosion regressions pass. The installed catalogue and raw/gzip bank replay with exact hashes over HTTP through a temporary bounded local server, which was then closed. No browser canvas inspection is claimed.

This updates existing represented tissue; it adds no distinct TA98 term or anatomical-completion count. Apex identity, state-specific opening configurations, coarse native contours, tissue-depth/filling-state review and female transfer remain unfinished. Female SHA256 remains `ad8ebd03b3ba815cf3105e8ed11649a812c08dade8da58e5d405573d433abe77`. Nothing was deployed externally.


## Installed male apex reference revision v7

A08.3.01.002 now has an explicit apex landmark on the actual bladder-to-median-umbilical-ligament attachment cap. It uses the ordinary anatomical label and selection controls. The stored point, host triangle and barycentric coordinates replay against the literal installed geometry bank with zero Float64 roundtrip error. Selecting its concept URL retains the apex identity and exact anchor. The existing host wall is selection context; this adds no separate tissue solid, triangles, geometry bank, texture or material. Clipping the landmark hides its callout rather than moving it to another surface. The actual CPU front/right surface views were inspected with the diagnostic callout; browser canvas inspection is not claimed.

The authoritative terminology correction places the female filling/voiding opening configurations under A08.3.01.028F; the male configurations remain under A08.3.01.027M. The raw terminology cache remains untouched. The primary bladder table also cross-references A04.5.03.017, rectovesicalis. Its canonical muscular parent is retained and a contextual bladder occurrence is provided. No distinct rectovesicalis geometry is present yet. The sex-aware inventory records 34 terms/context references, 30 applicable to each sex; this is an inventory count, not a completion count. Filling/voiding terms remain configurations of their existing opening, not additional holes or solids.

Male catalogue SHA256: `892f84f1084885e0b57aa917a628b3389cfd46b66a6afc319e583d94436caafe`. Installation ID: `c1a9ca5488b5f228e09ed713b436b7efd8c3916c9bc0b08a449b52f777f2f9ef`; rollback: `.local-models/ta98-review/male-urinary-apex-v7/manifest.json`. Evidence: `.local-models/geometry-workbench/bladder-completion/male-apex-reference-v1/`. All parts, chunks and materials equal lining revision v6 literally. All 6,784 registry parts and the complete ownership forest are current, with 888 fixed skeletal items. All 32 current urinary/pubic interface bindings pass literal coordinate and opposed-normal checks. The attachment closure retains 36 members and zero fixed skeletal members. Unrelated provisional declarations remain unaccepted.

The installation completed in 52.10 seconds under its 120-second cap. Registry refresh took 42.78/120 seconds, forest import 6.13/60 seconds and current-interface verification 5.03/60 seconds. The 38 installed urinary/hierarchy/camera/explosion tests pass in 4.69 seconds under a 45-second cap. All seven runtime-supervisor regressions pass in 1.96 seconds under a 20-second cap, covering deadline termination, inherited deadlines, cancellation and child cleanup. The initial run with system Python lacked the PostgreSQL driver; rerunning in the existing geometry environment corrected the environment without extending the cap. All commands retain separate elapsed-time records; no automatic retries or external deployment occurred.

State-specific opening configurations, rectovesicalis representation, coarse native contour and tissue-depth/filling-state review, and female transfer remain unfinished. Female catalogue SHA256 remains `ad8ebd03b3ba815cf3105e8ed11649a812c08dade8da58e5d405573d433abe77`.


## Opening configurations: actual junction source

Before constructing a filling/voiding variant, retain the complete v7 wall, renal junctions, prostate, urethra and fixed pelvic context. The male-urethra authority section defines A09.4.02.001 and .002 as states of the same opening: a closed intramural neck with a flat baseplate during filling, and an open funnel incorporating the neck during voiding. The existing installed geometry has an open funnel; a complete closed filling configuration has not been built.

The true bladder-to-prostatic-urethral lining junction is a single literal 45-edge loop on the reciprocal 108-triangle mucosal annulus. This is distinct from the previously bound generic neck landmark approximately 4.53 mm superior to it. The new source reference uses the actual common lumen edge, not a disk, proximity join, separately invented opening or a translated copy of that earlier landmark. The two complete tissue meshes and their opposed stored normal fields must replay before this loop can establish a configuration reference. Evidence is staged under `male-opening-states-v1`; no configuration is installed or anatomically accepted by this source extraction alone.


## Installed male voiding-reference revision v8

A09.4.02.002, the emptying internal urethral orifice, now selects the existing open-funnel configuration through the ordinary TA98 hierarchy, direct concept URL and rim-anchored anatomical callout. It reuses the parent clinical ostial mucosa and its 76 actual free-lining witnesses. Maximum encoded-facet replay error is 1.18e-16 m. The callout points to an actual rim sample; it is not dispersed to a generic layer point. No new tissue solid, coordinate bank, triangles, texture, material or independently superimposed hole is added. Selecting the configuration does not change anatomy or move the camera. This identifies the existing open-neck configuration; it is not a whole-bladder pressure/filling simulation.

The deeper histological bladder/urethral junction is independently extracted as one 45-edge loop on a reciprocal 108-triangle tissue annulus, with literal opposed render-normal words. Its vector-projected area is 12.579999 mm² and its maximum residual from that projection plane is 0.148705 mm; these are geometric diagnostics, not clinical axial measurements. The actual emitted-normal sagittal render shows this tissue join lies below the prostatic base, while the existing clinical ostium remains at the base. The proposal to use the deeper tissue join as the clinical opening is rejected and explicitly retained in `histological-plan-review.json`. The clinical rim is not relocated. The new reusable `shared_lumen_rim` helper rejects a micron-scale gap, reversed annular winding, mismatched render-normal words and ambiguous inner/outer free rims.

Installed male SHA256: `c9abc4bfeee0596c799d6d05b082d2273ad39e53157a70d7974db98d3abf7932`; revision `male-urinary-voiding-v8`; installation ID `1956b85598c2911a48ad6671e450efe3f82766d58be14f42f5ad4e47eb54a920`. Rollback: `.local-models/ta98-review/male-urinary-voiding-v8/manifest.json`. Preconstruction inventory: `source-v7/male-inventory.json`. Evidence: `.local-models/geometry-workbench/bladder-completion/male-opening-states-v1/`. All parts, materials and chunks remain literal relative to v7; female SHA256 remains `ad8ebd03b3ba815cf3105e8ed11649a812c08dade8da58e5d405573d433abe77`.

The guarded catalogue check took 0.83/30 seconds, installation 51.81/120 seconds, registry refresh 42.86/120 seconds, forest import 6.19/60 seconds and current interface replay 4.80/60 seconds. All 32 current urinary/pubic interface bindings pass; the attachment closure remains 36 members with zero fixed skeleton members. The forest retains all 6,784 parts and 888 fixed skeletal items. Unrelated provisional declarations remain unaccepted. All 42 installed viewer/hierarchy/camera/explosion tests and 14 geometry/source tests pass; TypeScript passes. The clinical witness guard rejects altered bank identities, a forged filling-state claim, moved witnesses and changed ownership. Local HTTP verification passes for the exact installed catalogue, the real Vite-resolved control module and the direct-URL label; its temporary server is closed afterward. This is source/module verification, not browser-canvas inspection. The initial HTTP module check used the wrong Vite root, and a standalone replay initially lacked its workspace import path; both checks were corrected under their original caps.

The closed filling configuration A09.4.02.001, rectovesicalis representation, remaining coarse contour/tissue-depth review and female transfer are still unfinished. The existing open funnel is explicitly rejected as evidence of a closed filling state. No external deployment or full-organ anatomical acceptance is claimed.


## Filling-state construction before deformation

The v8 installed open-neck chart is the construction source; source-v8-male-inventory.json binds the complete authority inventory to it. The filling candidate must close the neck through apposition of its own mucosal walls and form the flatter baseplate described by TA98. No independent lumen disk or extra plug is permitted. The complete 29-cell bladder arrangement is transported coherently. External ureter/urethral roots, prostate/support contacts, median ligament, pubic attachments and the complete fixed skeleton retain their literal geometry. Opening curves must be rebound to the resulting actual lining.

First certify the source lining as two nonoverlapping half-graphs in a frame fixed to the clinical opening. A bounded partial apposition is an intermediate numerical diagnostic only, not a filling-state representation or an installable candidate. Before complete closure, the ostial mucosa must be partitioned into independently finite apposing cells; any median-plane caps close tissue thickness only. The resulting contact must be declared explicitly and both upstream and downstream fluid boundaries checked in the closed configuration. Full source replay, normals, tissue volume, crossings, occupied-volume contacts and actual surface/section review remain mandatory before a state can be installed. Transition distances used during construction are illustrative parameters, not patient-specific physiological measurements.


The planar half-graph candidates are rejected: real clinical lining facets overlap in both principal-axis projections. Exact facet witnesses are retained in the v1/v2/v3 filling workbench folders; no tolerance is increased and none is installed. All 1,441 selected clinical facets replay literally on the already verified fluid boundary, so this is not permission to delete the posterior lining. A radial chart instead passes its strictly inward radial-normal test on all 1,681 clipped facets. The next numerical diagnostic uses a radius-squared apposition map driven by literal source-face ray intersections. Its untapered core preserves each tissue sector's cross-sectional area; an outward monotone C2 taper freezes the actual adjoining contacts. This is not yet a filling-state representation. Complete closure still needs topology changes at coaptation, two closed fluid domains, a reviewed baseplate and the actual finite-cell/normal/volume/contact gates. The unsuccessful upper-neck region-edge extraction does not establish an anatomical rim and is not used to position a completed state.

## Partial neck-apposition diagnostic

The initial half-area radial construction is rejected: its encoded ostial cell loses about 23.74% of its volume and has 132 self-crossing pairs. Its trigonal cell has three crossing pairs; four other changed cells were not checked after that failure. A positive continuous Jacobian did not establish valid emitted geometry. The original source is refined before reevaluating the same immutable lining-ray field, rather than subdividing already deformed triangles. Shared source-edge subdivisions propagate to every incident tissue surface, with barycentric transport of the original corner-normal field. All 18 external attachment groups, comprising 3,950 canonical facets and their source normal words, remain literal.

The first local refinement splits 4,675 source edges and adds 9,490 canonical triangles. Complete checking finds 33 ostial, seven trigonal, eight internal-longitudinal and 15 submucosal crossing pairs; this version remains rejected. Refining the actual crossing witnesses splits another 89 edges: seven of the eight changed cells pass, leaving four ostial pairs. A final ten-edge refinement removes those four pairs. `male-filling-radial-half-refined-v3/partial-verification.json` binds all eight passing actual Float32 cell checks, outward Int16 normal words, the unchanged 21 cells, the current source catalogue and the inspected uncapped sagittal render. Five byte-identical cell checks are reused only after their mesh hashes, implementation signature and source-volume bindings replay. No numerical crossing gate is weakened.

The inspected result is a localized, still-patent hourglass narrowing. It is explicitly uninstalled and is not accepted as the closed filling configuration or flat baseplate. Cross-cell occupied-volume checks and the new whole-fluid proof remain incomplete. Establish the actual upper-neck/baseplate and intramural contact extent before constructing complete coaptation; retain explicit finite tissue ownership and verify the two resulting fluid domains. Male remains revision v8, SHA256 `c9abc4bfeee0596c799d6d05b082d2273ad39e53157a70d7974db98d3abf7932`; female remains SHA256 `ad8ebd03b3ba815cf3105e8ed11649a812c08dade8da58e5d405573d433abe77`. Refinement stages completed in 6.79, 13.31 and 7.71 seconds under separate 60-second caps; construction, normal transport, cell checks and rendering retain their own finite runtime records. No model was installed or deployed by this diagnostic.
