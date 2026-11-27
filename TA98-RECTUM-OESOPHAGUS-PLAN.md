# Rectum and oesophagus: dual-model TA98 construction inventory

Complete authority descendants plus shared rectal wall terms. Sources are organ-specific correspondences, not anatomical acceptance. Shared wall terminology mapped on colon segments does not represent rectal tissue; those mappings are retained separately as sameTermOnOtherOrgans in the JSON report.

| Organ | TA98 | Parent | Term | Representation | Male sources | Female sources |
|---|---|---|---|---|---|---|
| oesophagus | A05.4.01.001 | A05.0.00.000 | oesophagus | whole organ or tissue region | `OESOPHWALL:male:native:external-connective-coat`<br>`OESOPHWALL:male:native:longitudinal-muscle`<br>`OESOPHWALL:male:native:circular-muscle`<br>`OESOPHWALL:male:native:submucosa`<br>`OESOPHWALL:male:native:muscularis-mucosae`<br>`OESOPHWALL:male:native:lamina-propria`<br>`OESOPHWALL:male:native:epithelium`<br>`OESOPHWALL:male:connector:external-connective-coat`<br>`OESOPHWALL:male:connector:longitudinal-muscle`<br>`OESOPHWALL:male:connector:circular-muscle`<br>`OESOPHWALL:male:connector:submucosa`<br>`OESOPHWALL:male:connector:muscularis-mucosae`<br>`OESOPHWALL:male:connector:lamina-propria`<br>`OESOPHWALL:male:connector:epithelium` | `OESOPHWALL:female:native:external-connective-coat`<br>`OESOPHWALL:female:native:longitudinal-muscle`<br>`OESOPHWALL:female:native:circular-muscle`<br>`OESOPHWALL:female:native:submucosa`<br>`OESOPHWALL:female:native:muscularis-mucosae`<br>`OESOPHWALL:female:native:lamina-propria`<br>`OESOPHWALL:female:native:epithelium` |
| oesophagus | A05.4.01.002 | A05.4.01.001 | cervical part | whole organ or tissue region | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.003 | A05.4.01.001 | thoracic part | whole organ or tissue region | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.004 | A05.4.01.003 | thoracic constriction | regional contour | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.005 | A05.4.01.003 | diaphragmatic constriction | regional contour | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.006 | A05.4.01.001 | abdominal part | whole organ or tissue region | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.007 | A05.4.01.001 | serosa | wall tissue | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.008 | A05.4.01.001 | subserosa | wall tissue | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.009 | A05.4.01.001 | adventitia | wall tissue | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.010 | A05.4.01.001 | muscular layer | wall tissue | `OESOPHWALL:male:native:longitudinal-muscle`<br>`OESOPHWALL:male:native:circular-muscle`<br>`OESOPHWALL:male:connector:longitudinal-muscle`<br>`OESOPHWALL:male:connector:circular-muscle` | `OESOPHWALL:female:native:longitudinal-muscle`<br>`OESOPHWALL:female:native:circular-muscle` |
| oesophagus | A05.4.01.011 | A05.4.01.010 | crico-oesophageal tendon | attachment tissue | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.012 | A05.4.01.010 | broncho-oesophageus | attachment tissue | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.013 | A05.4.01.010 | pleuro-oesophageus | attachment tissue | Missing distinct representation | Missing distinct representation |
| oesophagus | A05.4.01.014 | A05.4.01.001 | submucosa | wall tissue | `OESOPHWALL:male:native:submucosa`<br>`OESOPHWALL:male:connector:submucosa` | `OESOPHWALL:female:native:submucosa` |
| oesophagus | A05.4.01.015 | A05.4.01.001 | mucosa | wall tissue | `OESOPHWALL:male:native:muscularis-mucosae`<br>`OESOPHWALL:male:native:lamina-propria`<br>`OESOPHWALL:male:native:epithelium`<br>`OESOPHWALL:male:connector:muscularis-mucosae`<br>`OESOPHWALL:male:connector:lamina-propria`<br>`OESOPHWALL:male:connector:epithelium` | `OESOPHWALL:female:native:muscularis-mucosae`<br>`OESOPHWALL:female:native:lamina-propria`<br>`OESOPHWALL:female:native:epithelium` |
| oesophagus | A05.4.01.016 | A05.4.01.015 | muscularis mucosae | wall tissue | `OESOPHWALL:male:native:muscularis-mucosae`<br>`OESOPHWALL:male:connector:muscularis-mucosae` | `OESOPHWALL:female:native:muscularis-mucosae` |
| oesophagus | A05.4.01.017 | A05.4.01.001 | oesophageal glands | microscopic structure | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.01.002 | A05.7.01.001 | serosa | wall tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.01.003 | A05.7.01.001 | subserosa | wall tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.01.004 | A05.7.01.001 | muscular layer | wall tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.01.005 | A05.7.01.001 | submucosa | wall tissue | `LOCAL:male:CAECOWALL:rectum:submucosa` | `LOCAL:female:CAECOWALL:rectum:submucosa` |
| rectum | A05.7.01.006 | A05.7.01.001 | mucosa | wall tissue | `LOCAL:male:CAECOWALL:rectum:muscularis-mucosae`<br>`LOCAL:male:CAECOWALL:rectum:lamina-propria`<br>`LOCAL:male:CAECOWALL:rectum:epithelium` | `LOCAL:female:CAECOWALL:rectum:muscularis-mucosae`<br>`LOCAL:female:CAECOWALL:rectum:lamina-propria`<br>`LOCAL:female:CAECOWALL:rectum:epithelium` |
| rectum | A05.7.01.007 | A05.7.01.006 | muscularis mucosae | wall tissue | `LOCAL:male:CAECOWALL:rectum:muscularis-mucosae` | `LOCAL:female:CAECOWALL:rectum:muscularis-mucosae` |
| rectum | A05.7.01.008 | A05.7.01.006 | intestinal glands | microscopic structure | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.001 | A05.7.01.001 | rectum | whole organ or tissue region | `LOCAL:male:CAECOWALL:rectum:serosa`<br>`LOCAL:male:CAECOWALL:rectum:subserosa`<br>`LOCAL:male:CAECOWALL:rectum:longitudinal-muscle`<br>`LOCAL:male:CAECOWALL:rectum:circular-muscle`<br>`LOCAL:male:CAECOWALL:rectum:submucosa`<br>`LOCAL:male:CAECOWALL:rectum:muscularis-mucosae`<br>`LOCAL:male:CAECOWALL:rectum:lamina-propria`<br>`LOCAL:male:CAECOWALL:rectum:epithelium` | `LOCAL:female:CAECOWALL:rectum:serosa`<br>`LOCAL:female:CAECOWALL:rectum:subserosa`<br>`LOCAL:female:CAECOWALL:rectum:longitudinal-muscle`<br>`LOCAL:female:CAECOWALL:rectum:circular-muscle`<br>`LOCAL:female:CAECOWALL:rectum:submucosa`<br>`LOCAL:female:CAECOWALL:rectum:muscularis-mucosae`<br>`LOCAL:female:CAECOWALL:rectum:lamina-propria`<br>`LOCAL:female:CAECOWALL:rectum:epithelium` |
| rectum | A05.7.04.002 | A05.7.04.001 | sacral flexure | regional contour | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.003 | A05.7.04.001 | lateral flexures | regional contour | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.004 | A05.7.04.003 | superodextral lateral flexure | regional contour | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.005 | A05.7.04.003 | intermediosinistral lateral flexure | regional contour | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.006 | A05.7.04.003 | inferodextral lateral flexure | regional contour | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.007 | A05.7.04.001 | transverse folds of rectum | inner-wall relief | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.008 | A05.7.04.001 | rectal ampulla | regional contour | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.009 | A05.7.04.001 | muscular layer | wall tissue | `LOCAL:male:CAECOWALL:rectum:longitudinal-muscle`<br>`LOCAL:male:CAECOWALL:rectum:circular-muscle` | `LOCAL:female:CAECOWALL:rectum:longitudinal-muscle`<br>`LOCAL:female:CAECOWALL:rectum:circular-muscle` |
| rectum | A05.7.04.010 | A05.7.04.009 | longitudinal layer | wall tissue | `LOCAL:male:CAECOWALL:rectum:longitudinal-muscle` | `LOCAL:female:CAECOWALL:rectum:longitudinal-muscle` |
| rectum | A05.7.04.011 | A05.7.04.010 | rectococcygeus | attachment tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.012 | A05.7.04.010 | anorectoperineal muscles | attachment tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.013 | A05.7.04.012 | rectoperinealis | attachment tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.014 | A05.7.04.012 | anoperinealis | attachment tissue | Missing distinct representation | Missing distinct representation |
| rectum | A05.7.04.015 | A05.7.04.009 | circular layer | wall tissue | `LOCAL:male:CAECOWALL:rectum:circular-muscle` | `LOCAL:female:CAECOWALL:rectum:circular-muscle` |
| rectum | A05.7.04.016 | A05.7.04.001 | lateral ligament of rectum | attachment tissue | Missing distinct representation | Missing distinct representation |

## Assembly/interface plan

- oesophagus: hypopharynx/cricoid → cervical part
- oesophagus: thoracic inlet → thoracic part
- oesophagus: aortic/bronchial constriction; not an artificial solid sphincter
- oesophagus: diaphragmatic hiatus → abdominal part
- oesophagus: distal lumen and tissue interfaces → cardia
- rectum: sigmoid lumen/wall → rectosigmoid junction
- rectum: sacral/lateral flexures and ampulla: shared regions of the same tube
- rectum: transverse folds: relief carried by mucosa and submucosa
- rectum: puborectalis/pelvic floor and mesorectal neurovascular interfaces
- rectum: rectal outlet → anal canal; lumen stays open

## Construction and efficiency

- Use one shared carrier with ordered inner/outer tissue boundary fields.
- Generate tissue appearance procedurally; retain coherent UVs only where needed.
- Peel cells independently and cap tissue intervals without capping the lumen.
- Oesophageal serosa/subserosa are regional peritoneal coverage, not full-length adventitia replacements.
- Rectum has a continuous longitudinal coat; do not give it colonic taeniae.
- Supplementary epithelium/lamina propria/oesophageal circular and longitudinal coats are distinct from exact TA98 entries.
- Do not invent glands, attachments or contour landmarks from an arbitrary partition.

## Placement and dependency constraints

- Preserve accepted stomach geometry and native sigmoid/anal openings.
- Rectum, its actual attachments and shared wall cells must use the correct source-frame forest group.
- Retain vessels and nerves within wall interfaces; never pass the luminal boundary.
- Known pelvic and oesophageal neighbor defects remain unresolved and require repair.

## Current construction evidence

Female oesophageal and both rectal candidates each contain seven finite cells sharing eight interfaces. Their independent field replay and local geometry gates pass; these are staged constructions, not installed/accepted organs. Male oesophageal layer correspondence retains every exterior source control/face and both original distal wall rings, with about 2.3 MB of shared fields and no image assets. Its full global cell checks and inner-chart approximation review remain pending. The ordinary male atlas already includes a distal gastric connector; both sections must form one oesophageal assembly and its shared join must suppress an isolated end rim.

Blender render/section evidence and hashes are in `reports/rectum-oesophagus-staged-render-review.json`. Distinct tissue colors and patent rectal ends are visible. The male annular rim independently replays as source-exact tissue closure with an area equal to outer minus inner ring area. No lumen disk is added. Regional peritoneal allocation, rectal folds/flexures/ampulla, glands, named attachments, body packing, runtime installation and complete hierarchy controls remain unfinished. Historical inventory snapshots retain the original construction-report hash bindings. Shared coat mappings on other colon organs are not rectal representations.

Reference constraint: oesophageal peritoneal coverage is limited to the abdominal segment (https://pmc.ncbi.nlm.nih.gov/articles/PMC7263794/); identify the actual hiatus and cardial source interfaces before allocating serosa/subserosa. Do not use the entire diaphragm bounding box as a hiatus landmark.

The seam-preserving male inner chart independently measures all 106 native planes: maximum sampled bidirectional vertex-to-segment profile error is 0.0129143 mm. This is a source-plane approximation measurement, not clinical wall accuracy or a continuous 3D error bound. Full global cell checks remain pending; the superseded ray-only chart passed its first complete cell check before that job was stopped to preserve the original distal join.

## Continued whole-organ repair

All seven seam-preserving male oesophageal cells now pass the complete topology, local prism and self-crossing checks (`oesophageal-layer-correspondence-seam-preserved-replay.json`). The broad phase uses a conservative two-ended interval search; actual finite-wall and penetrating-duplicate candidate/crossing pairs independently match exhaustive comparison (`wall-intersection-sweep-replay.json`). This changes evidence computation only.

Female gastric conformity now uses the current accepted gastric bank rather than the older wall probe. Four boundary-adjacent source refinements make a 17-control common mouth. A common 20 mm source-surface support harmonic solve adjusts the eight boundary fields, maximum displacement 0.106826 mm. The final outer/inner mouth coordinates are copied literally from the gastric bank; no stomach coordinate changes. Independent source interpolation, fixed support and harmonic replay are exact. Seven closed tissue cells and the whole joined macro wall have zero detected self-crossings, one closed genus-one wall and no topology faults. True plane-clipped sagittal junction and complete exterior inspected in existing Blender; no lumen-blocking internal rim. Images/hashes in the conforming render-review report.

The current native atlas screen still finds inner-boundary crossings with oesophageal arteries (24 face pairs), right intercostal arteries (8), left vagus (24), and liver (43). Wall-entry contacts remain legitimate candidates; these innermost contacts do not pass lumen exclusion. None of these structures is hidden or deleted. Clinical wall allocation, regional parts, distinct glands/attachments, shared normals/runtime caps and body packing remain unfinished.

Current interface diagnostics also show that the retained translated male native rectum does not match the current sigmoid: maximum nearest-control distance is 23.3204 mm. This is a sampling diagnostic, not a continuous surface-gap measurement. Its whole adjoining source family must be fitted coherently before catalogue replacement; do not silently append a second rectum or insert a coarse bridge. Both ordinary atlas catalogues remain unchanged.

The completed packing screen also reports right vagus innermost crossings (5 face pairs) and cricoid-cartilage crossings (32). Thus the proximal oesophageal/pharyngeal boundary needs landmark correction as well as the distal functional packing. These are additional unresolved contacts, not accepted wall entries. The full contact set is retained in the packing JSON.

## Existing male lower connector layer assembly

Continue the existing installed `GASTRIC:male:cardiac-connection` as the distal portion of the same oesophagus. Use its literal outer and luminal coordinates and indexed ten-station carrier; no replacement course, overlapping tube or added solid opening. Eight shared depth fields meet all eight native oesophageal fields at the upper end. At the accepted gastric mouth their allocation meets gastric interfaces 0, 2, 3, 5, 6, 7, a supplementary mucosal interpolation, and 8. Thus the oesophageal circular coat meets the gastric circular/oblique interval without a barrier. The combined external connective interval remains unallocated pending actual peritoneal/hiatus evidence. Internal allocation blends smoothly along source station arc length while the exterior/lumen remain fixed. Conditional section/disassembly rims cover only visible tissue intervals; no assembled internal rim is allowed. Validate every finite connector cell and the complete welded macro wall before any runtime replacement. Clinical thickness, source registration/packing and microscopy are separate unresolved obligations.

## 2026-10-08 male connector continuation evidence

The existing ten-station connector now has seven staged finite cells sharing eight boundary fields (780,480 field bytes; zero new images/UVs). The first absolute-offset allocation was rejected: its muscularis-mucosae interval reversed locally and produced 5,862 self-crossing face pairs. Proportional depth interpolation within each actual paired wall profile fixes that allocation while keeping the original connector exterior/lumen and all native/gastric endpoint coordinates literal. All seven corrected connector cells pass closed genus-one topology, strict negative continuous local prism checks and zero detected self-crossings. Independent field/face/source replay passes. The welded complete native oesophagus–connector–accepted gastric macro wall is one closed genus-one wall with zero detected crossings and no internal junction rim. Its two physiological mouths remain lumen openings, with tissue-only terminal closure. The normalized allocation reuses that expensive macro check only after exact comparison of every outer/inner/index/edge/origin array.

Actual complete exterior and triangle-plane-clipped gastric-entry views were rendered in the already-open Blender, with prior scenes/files retained. The rendered preview arrays match the checked complete macro bank exactly. The existing connector still visibly has an abrupt upper waist/shoulder: projected principal widths at source station 8 are 7.783 × 7.470 mm versus 16.593 × 16.594 mm at the native station 9 interface. These source-plane diagnostics are not clinical lumen dimensions. The coarse gastric-entry contour also remains visible. The layer and mathematical join checks therefore do not establish anatomical acceptance. Repair the complete connector course/profile and body packing before ordinary-atlas replacement; do not disguise the waist by merely changing material or calling it a constriction.

Source and result evidence is in `reports/ta98-male/oesophageal-layered-connector{,-normalized,-normalized-replay,-render-review}.json`. The original failed allocation stays separate from the corrected bank. Both ordinary male/female atlas catalogue hashes remain unchanged. Rectal sigmoid/anal interfaces, pelvic packing and named rectal parts remain unfinished. Current female pharynx source bounds also stop well above the upper oesophageal source; cutting at the cricoid lower border alone cannot establish a complete upper alimentary junction.

## Upper connector waist repair plan

Use the retained pre-transfer native oesophageal sheet geometry from `male-gastric-transfer-install-candidate/male-before-install.json` as the actual profile reference. Separate its existing outer/luminal sheets by connectivity after excluding constant-height annular faces; record source triangle/edge interpolation for actual transverse profiles. Keep the accepted gastric mouth and current native upper mouth literal. Replace the connector's upper narrowed profiles with the native paired profiles at their actual world heights. The lower continuation blends the gastric paired profile into the native profile through a common proper plane rotation and a smooth center course with endpoint tangents. Every depth boundary follows the same source-constrained profile field; no new opening solid, duplicated oesophagus, or stomach modification. Check continuous ordered prisms first, then actual joined wall topology/self-crossings, rendered sections, and all eligible active native neighbors including inner vessel/nerve contacts. A wider tube that crosses adjacent tissue is a rejected proposal, even if the mathematical wall is valid. Preserve the original and the normalized-layer banks separately.

## Current local installation and female rectal donor fitting — 9 October 2026

The user authorized installation of both current oesophageal courses as-is, with known neighboring-organ defects retained. Both ordinary local TA98 catalogues now contain those courses; public delivery is unchanged.

For rectum, use the corrected native female shared carrier as a possible male donor, retaining its source-face ownership and wall boundary fields. Compare it with each sex's fixed sacrum/coccyx, bladder/prostate or cervix/vagina, puborectalis and external sphincters. The closed male FJ2571 source is an alternative to replace, never an additional rectum. The female HRA rectum is a separately registered alternative, never an adjoining native tissue.

The proximal interface is the existing 31-edge caecocolic rim. Reconcile its seven source edge subdivisions against the native female 24-edge inlet, using one coordinate per common seam. The downstream eleven-edge rectal interface remains shared with the native anal source; the source's distal perianal margin must not automatically acquire digestive wall terminology. Preserve the fixed skeletal and muscle guides. All wall cells derive from one exterior carrier; continuous longitudinal muscle replaces colonic taeniae below the rectosigmoid junction. Use finite tissue-only annular caps with the same diagonal on both sides, conditional on adjoining tissues, preserving the lumen. Differentiate the coats procedurally without copied textures. Covering distribution, calibrated thickness, microscopic glands and transverse fold anatomy require distinct evidence; do not turn a generic outer shell into an exact representation of them.

Before promotion: replay source coordinates and ancestry; test every wall prism and surface pairing, topology/components/tunnels, adjoining seams, pelvic packing and vessels against the innermost boundary; inspect actual exterior and cuts. Keep failed trials local and record their cause. Rebase the structural forest against the newly installed oesophageal catalogues while preserving historical evidence.

### Installed local wall result

Both ordinary TA98 atlases now replace their prior rectum with eight finite coats in one shared caecocolic/rectal coordinate and normal bank. The female donor is resized and fitted for the male, preserving the larger male transverse and anterior–posterior extent while matching the fixed pelvic endpoints. The source rectum and female HRA alternative are suppressed through explicit replacement scopes; no second rectum is added. Existing part order and unrelated meshes, visibility defaults, skeleton and muscle guides are preserved.

The installed Float32 fields pass all eight continuous interval tests and 45 rectal sheet-pair checks per sex. Reciprocal 31-edge sigmoid caps use the identical common diagonal and hide at assembled joins; no disk fills the lumen. Main-atlas selection exposes eight peel stages and the TA98 major muscular/mucosal subtrees. Actual Blender exterior and sagittal cuts, plus the male main-atlas cut and peel, were inspected. Both section and peel changes preserve the exact observed camera. Tests include ordinary legacy URL prefix restoration, reciprocal interfaces, peeling, camera intent and explosion plane/viewport containment.

This is an installed improvement for local inspection, **not a perfected anatomical acceptance**. Full-body screens still record vessels crossing the innermost surface and other bowel/soft-tissue conflicts. Named rectal folds, glands, attachment structures, regional peritoneal allocation and calibrated tissue thickness remain unfinished. The native perianal margin is not reclassified as digestive coats. Evidence: `reports/installed-oesophageal-models-as-is.json`, `reports/installed-rectal-donor-models.json`, `reports/installed-rectal-fields-validation.json`, `reports/rectal-donor-assembly-female-validation.json` and `reports/rectal-donor-assembly-male-validation.json`. The complete current dependency forests were refreshed to include the shared rectal assembly and its retired alternatives.

## Shared proximal rectal course continuation — 9 October 2026

The current forest contains the complete installed caecocolic/rectal bank and retained source alternatives. Before further movement, distinguish the default-hidden BodyParts3D distal-ileum pieces and old whole Skin from the ordinary displayed anatomy; they remain separate references and are not deleted to improve a collision count. The male live lumen still crosses the installed jejunal wall, left iliopsoas fascia and puboanalis source. Female native vessel contacts remain separate wall-entry/innermost-boundary obligations.

Test a bounded posterior displacement at the male sigmoid–rectum interface, propagated along the sigmoid and proximal rectum by source-carrier geodesic distance. Move every stored tissue boundary through the same field and retain the literal common seam. Fix all controls outside the bounded support, the distal rectal/anal mouth, skeleton, pelvic-floor muscles and vascular sources. Reject inverted wall prisms, new carrier/self crossings, failed tissue/section topology or newly intruding skeletal/urinary surfaces. Preserve source IDs, face ownership, source alternatives, existing shared cap diagonals and runtime ordering. Repeat packing against actual ordinary tissues and classify hidden alternatives separately. This continuation grants no anatomical acceptance by a reduced count alone.

### Paracolic source coordination inventory

The source name combines regional terminology; both applicable mesenteric node subtrees are inventoried before fitting. No exact new TA98 mapping or physical attachment is inferred from that name. The complete source moves through the common spatial field; its tissue geometry, source classification and surrounding registration remain under review.

| TA98 | Parent | Term | Current explicit source |
|---|---|---|---|
| A13.3.03.033 | A13.3.03.014 | superior mesenteric nodes | Missing direct representation |
| A13.3.03.034 | A13.3.03.033 | juxta-intestinal mesenteric nodes | Missing direct representation |
| A13.3.03.035 | A13.3.03.033 | central superior mesenteric nodes | Missing direct representation |
| A13.3.03.036 | A13.3.03.033 | ileocolic nodes | Missing direct representation |
| A13.3.03.037 | A13.3.03.033 | precaecal nodes | Missing direct representation |
| A13.3.03.038 | A13.3.03.033 | retrocaecal nodes | Missing direct representation |
| A13.3.03.039 | A13.3.03.033 | appendicular nodes | Missing direct representation |
| A13.3.03.040 | A13.3.03.033 | mesocolic nodes | Missing direct representation |
| A13.3.03.041 | A13.3.03.040 | paracolic nodes | Missing direct representation |
| A13.3.03.042 | A13.3.03.040 | colic nodes | Missing direct representation |

### Superior rectal arterial source course — before fitting

The source is a single open tubular carrier, with two twelve-edge mouths. It has no modeled vascular inner surface or inferred additional branch; fitting changes its course only. Preserve the original source mouth coordinates and a proximal fixed collar. Nearest parent-mesh vertex distances range from 0.12 to 0.49 mm; this does not establish a conforming tissue join. All applicable inferior-mesenteric authority descendants are inventoried below. No additional exact term is claimed from this fit.

| TA98 | Parent | Term | Existing explicit mapping |
|---|---|---|---|
| A12.2.12.069 | A12.2.12.001 | inferior mesenteric artery | Native `ZA:Inferior mesenteric artery`; FMA14750 is a reference set, not the native whole trunk |
| A12.2.12.070 | A12.2.12.069 | ascending artery | Candidate `ZA:Ascending branch of left colic artery`; distinct authority correspondence remains unverified |
| A12.2.12.071 | A12.2.12.069 | left colic artery | Native `ZA:Left colic artery`, separate ascending/descending branch carriers; full branch interfaces unverified |
| A12.2.12.072 | A12.2.12.069 | sigmoid arteries | Native `ZA:Sigmoid arteries`; complete branch tissue and lumen walls unrepresented |
| A12.2.12.073 | A12.2.12.069 | superior rectal artery | Active native `ZA:Superior anorectal artery`; FMA14832 currently references hidden alternative `FJ3656` |

The fitted native superior anorectal arterial carrier will replace the hidden BodyParts3D alternative as the selected source of FMA14832 / A12.2.12.073 in the local male catalogue. Keep FJ3656's original packed attributes recoverable and redirect its selections through an explicit replacement scope; do not render two registered alternatives together. This is a source correspondence and course repair, not a new complete vascular wall or anatomical acceptance. Other arterial reference sets and unrelated source meshes remain unchanged. Both source end profiles and four complete sections at each end remain fixed. The fit must clear the entire epithelial coat as well as the lumen, and introduce no new nonrectal structural contacts. Rebuild the current movement forest and digestive Visit bindings after installation.

The course is now installed in the male master and ordinary local wall catalogue. Source/epithelial surface crossings decreased from 52 to zero; vertex, face-centroid and edge-midpoint samples find no lumen intrusion. Both mouths and their four-section collars are exact, with all 2,304 original faces and source UVs retained. Whole-source self-crossings remain zero. No new nonrectal structural contacts were detected; inherited contacts with vessels and pubococcygeus remain. Independent packed replay checks both installations and proves unrelated parts/materials and rectal banks unchanged. The suppressed alternative's default-hide entry is removed, preventing its replacement scope from hiding the active native artery.

Actual source/fitted Blender renders and the installed atlas surface/sagittal section were inspected. Section and Visit-panel changes preserved the exact camera. The complete digestive itinerary remains available in both sexes after route regeneration. Current male ownership forest was regenerated from the installed catalogue. Nineteen arterial/rectal delivery and Visit tests, TypeScript, camera-intent and explosion-plane/viewport tests pass. This increment does not complete rectal folds, peritoneal allocation, vascular tissue junctions or the remaining female vascular/lumen conflicts. See the local `installed-male-rectal-arterial-course-validation.json` and installed atlas image.

### Female superior rectal artery — course plan before repair

The same five inferior-mesenteric authority terms listed above apply to the female arterial tree. The native sources are `LOCAL:female:inferior_mesenteric_artery:1685`, `LOCAL:female:ascending_branch_left_colic_artery:1671`, `LOCAL:female:descending_branch_left_colic_artery:1672`, `LOCAL:female:sigmoid_arteries:1669`, and `LOCAL:female:superior_rectal_artery:1687`. The two left-colic branch carriers do not establish a separately modeled complete left-colic trunk or a verified distinct ascending-artery authority correspondence. The superior rectal source has 1,434 faces, one exact-coordinate connected component and one six-edge proximal opening. Two exactly zero-area source faces cause the two nonmanifold edges and three invalid vertex links; removing only those empty faces leaves 1,432 real faces, one component, one proximal opening, and no nonmanifold, degenerate or invalid-link defects. Whole-source self-crossings are zero. Distal source terminations remain closed and are not interpreted as verified physiological vascular junctions.

The whole native source remains the A12.2.12.073 representation, not an invented set of new arterial branches or wall layers. Preserve all nonzero source faces, source UVs and materials. The proximal ring has six exact shared coordinates with the unchanged inferior mesenteric artery; keep that ring and its nearby collar fixed. Fit only the interfering distal courses using a smooth displacement field shared by coincident source coordinates. Exclude the entire rectal epithelial coat and lumen, retain section width as far as source sampling establishes it, and reject new self-crossings, collapsed faces or worsened neighboring-structure contacts. Inspect source/fitted geometry and actual rectal sections before installation. Existing uterine and vaginal plexus defects remain separate repair obligations; do not move a whole plexus merely to conceal an intrusion.

The first smooth trial cleared the rectal lumen but introduced iliac-vein and vaginal contacts; it is retained as rejected evidence. The revised compact field constrains the actual neighboring tissue sheets and keeps their original contact pairs unchanged. Bowel-wall contacts are classified according to the user's vascular requirement: sigmoid and rectal coats may host the existing arterial carrier, while neither epithelial coat nor the joined lumen may be penetrated. This does not establish a conforming vascular tissue junction. Test the combined sigmoid–rectal epithelial surface with their literal shared mouth coordinates, avoiding a measurement cap at that internal interface. Record all outer-coat contacts rather than treating them as anatomical acceptance; reject any new non-bowel tissue or vascular contact. All six native parent join coordinates, source collar and untouched organs remain fixed.

At the user's direction to finish in the best present form, fitting stops here. Both local TA98 models retain their installed eight-coat rectal assemblies and corrected superior rectal arterial courses. The final female fit moves at most 1.58 mm, preserves its six native parent coordinates/collar, and clears both epithelial boundaries and the combined sigmoid–rectal lumen (63 inner-surface crossings to zero; 244 combined-lumen vertex/face/edge samples to zero). No new non-bowel contacts or source self-crossings were detected. All three female catalogues deliver this exact candidate with original UVs/materials and unchanged surrounding geometry. The retired HRA rectum's inherited hide entry is removed so it cannot hide the replacement wall. Digestive Visit bindings and the female ownership forest are refreshed. Actual Blender source/fitted/context and atlas views were checked. Remaining folds, microscopic features, vascular junctions and other plexus defects are recorded as unaccepted; no further rectal fitting is authorized by this close-out.
