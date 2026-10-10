# Cricoid cartilage: TA98 inventory and construction boundaries

Inventory before geometry work, 10 October 2026. Authority: bundled TA98
metadata/parent index, originally derived from the TA98 terminology database;
all ten descendants of A06.2.03.001 are retained below. Anatomical reference:
[NCBI laryngeal anatomy](https://www.ncbi.nlm.nih.gov/sites/books/NBK538202/).

| TA98 | Term | Parent | Representation | Current source / missing representation |
|---|---|---|---|---|
| A06.2.03.001 | Cricoid cartilage | A06.2.02.001 | Whole solid cartilage ring around the airway | Male ZA:Cricoid cartilage; female LOCAL:female:cricoid_cartilage:1072 |
| A06.2.03.002 | Arch of cricoid cartilage | A06.2.03.001 | Anterior cartilage tissue region | Not separately partitioned in either source |
| A06.2.03.003 | Lamina of cricoid cartilage | A06.2.03.001 | Broad posterior cartilage tissue region | Not separately partitioned in either source |
| A06.2.03.004 | Arytenoid articular surface | A06.2.03.001 | Paired surface regions, not separate solids | Facets unpartitioned; named arytenoids exist in both body-specific frames |
| A06.2.03.005 | Thyroid articular surface | A06.2.03.001 | Paired lateral surface regions, not separate solids | Facets unpartitioned; named thyroid cartilage exists in both frames |
| A06.2.03.006 | Cricothyroid joint | A06.2.03.001 | Bilateral adjoining articulation assembly | Named cartilage relationship exists; joint tissues/mechanics unrepresented |
| A06.2.03.007 | Capsule of cricothyroid joint | A06.2.03.006 | Finite connective tissue covering | Missing in both models |
| A06.2.03.008 | Ceratocricoid ligament | A06.2.03.006 | Ligament tissue | Missing in both models |
| A06.2.03.009 | Median cricothyroid ligament | A06.2.03.001 | Separate anterior ligament tissue | Male ZA:Median cricothyroid ligament; female LOCAL:female:cricothyroid_ligament:1078; tissue extent unverified |
| A06.2.03.010 | Cricotracheal ligament | A06.2.03.001 | Separate adjoining ligament tissue | No verified distinct source in this construction |

Neither source establishes perichondrium, microscopic cartilage organization,
joint capsules, entheses or an airway mucosal wall. The airway opening is empty
space, not a cartilage plug. Never cap a regional boundary across that opening.

The female whole surface has one closed genus-one component, 362 triangles and
zero detected self-crossings. The male whole source has 1,882 triangles, four
nonmanifold edges, three opposite-winding duplicate face pairs and seven
detected crossing pairs. The defects are localized at an anterior midline
support. Repair the complete male carrier before partitioning either model.

First bounded repair: remove both sides of each exact opposite-winding
duplicate support (source faces 846/1727, 920/1859, 937/1875), only if complete
topology confirms they are unsupported two-sided fins rather than real tissue
walls. All remaining source faces, positions, normals, UVs and materials must
stay literal. Recheck one connected ring, intended airway tunnel, closed
surface, zero crossings, source-normal aliases, actual whole render and cuts.
Compare neighboring anatomy: unchanged retained surfaces must not acquire
new contacts. Do not move thyroid/arytenoid cartilages, ligaments, trachea,
oesophagus, skeleton or rectum for this repair.

Later arch/lamina work must first locate the anatomical lateral transitions
against the posterior plate and the adjoining articular landmarks. A blind
bounding-box split is insufficient. Partition one repaired carrier without
duplicating tissue faces; retain exact source coordinates/normals at shared
boundaries and use tissue-only conditional rims for disconnected presentation.
Whole selection retains the original identity and camera. Each regional part
gets its own hierarchy visibility/selection controls. Unverified boundaries
remain provisional until actual views and sections are reviewed. No new
images or duplicated full-organ buffers are planned.
