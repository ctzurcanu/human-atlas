# Male hepatic fitting dependencies

Source-bound movement-review inventory; no anatomy moved or accepted.
Full hepatic authority descendant inventory: [TA98-HEPATIC-PATHS.md](TA98-HEPATIC-PATHS.md).

| Source | TA98 | Role | Suppressed |
|---|---|---|---|
| ZA:Anterior lateral segment of liver (VI) | A05.8.01.050 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Anterior medial segment of liver (V) | A05.8.01.047 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Bile duct | A05.8.02.013 | vascular or ductal attachment; course/endpoints need review | True |
| ZA:Common hepatic artery | A12.2.12.015 | vascular or ductal attachment; course/endpoints need review | False |
| ZA:Cystic node | A13.3.03.031 | lymphatic tissue; hepatic/biliary adjacency requires review | False |
| ZA:Gallbladder | A05.8.02.001 | adjoining hollow organ; hepatic fossa and cystic duct interfaces require review | True |
| ZA:Hepatic portal vein | A12.3.12.001 | vascular or ductal attachment; course/endpoints need review | False |
| ZA:Hepatic veins | A12.3.09.005 | vascular or ductal attachment; course/endpoints need review | False |
| ZA:Inferior vena cava (abdominal part) | unmapped | vascular or ductal attachment; course/endpoints need review | False |
| ZA:Inferior vena cava (thoracic part) | unmapped | vascular or ductal attachment; course/endpoints need review | False |
| ZA:Left anterior lateral segment of liver (III) | A05.8.01.040 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Left medial segment of liver (IV) | A05.8.01.042 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Left posterior lateral segment of liver (II) | A05.8.01.039 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Lesser omentum | A10.1.02.101 | cross-organ attachment; gastric and hepatic endpoints require review | False |
| ZA:Liver:Ligament | A05.8.01.001 | original ligament material region; individual ligament identities unverified | False |
| ZA:Liver:Organ | A05.8.01.001 | whole-organ source exterior; open boundaries remain unverified | False |
| ZA:Posterior lateral segment of liver (VII) | A05.8.01.051 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Posterior medial segment of liver (VIII) | A05.8.01.048 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Posterior segment of liver (I) | A05.8.01.044 | overlapping source territory; preserve alignment, not an exclusive wall partition | False |
| ZA:Proper hepatic artery | A12.2.12.029 | vascular or ductal attachment; course/endpoints need review | False |

Exact shared source-coordinate constraints:

| Sources | Shared coordinates |
|---|---|
| ZA:Anterior lateral segment of liver (VI) ↔ ZA:Anterior medial segment of liver (V) | 31 |
| ZA:Anterior medial segment of liver (V) ↔ ZA:Posterior medial segment of liver (VIII) | 9 |
| ZA:Anterior medial segment of liver (V) ↔ ZA:Posterior segment of liver (I) | 16 |
| ZA:Liver:Ligament ↔ ZA:Liver:Organ | 100 |
| ZA:Posterior lateral segment of liver (VII) ↔ ZA:Posterior medial segment of liver (VIII) | 1 |

Coordinate coincidence does not prove a joint or establish a wall, lumen or exclusive territory. Suppressed donors remain suppressed. Cross-organ attachments are review constraints, not automatic rigid movement edges.

Before fitting, resolve the support region at the gastric impression and preserve hepatic inflow/outflow, vena-caval, ligament and omental interfaces. Missing ligament identities, branch ostia and the open liver source remain explicit limitations.

Current catalogue SHA256: 9b5d7adb28b2f6fb327eeee38d414aff71944aa0bc67eeb662d4c186f101759e. Local detailed evidence: reports/ta98-male/hepatic-fit-dependencies.json.
