# Humerus source interface repair

Inventory includes every authoritative descendant before geometry modification. The bone/articular-cartilage boundary is an adjoining tissue interface, not an inferred wall layer.

| TA98 | Parent | Class | Male mapped parts | Female mapped parts |
| --- | --- | --- | --- | --- |
| A02.4.04.001 humerus | A02.4.03.001 | whole bone | ZA:Humerus.l:Bone, ZA:Humerus.l:Cartilage, ZA:Humerus.r:Bone, ZA:Humerus.r:Cartilage | LOCAL:female:l_humerus:957, LOCAL:female:r_humerus:952 |
| A02.4.04.002 caput humeri | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.003 collum anatomicum | A02.4.04.001 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.004 collum chirurgicum | A02.4.04.001 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.005 tuberculum majus | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.006 tuberculum minus | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.007 sulcus intertubercularis | A02.4.04.001 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.008 crista tuberculi majoris | A02.4.04.005 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.009 crista tuberculi minoris | A02.4.04.006 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.010 corpus humeri | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.011 facies anteromedialis | A02.4.04.010 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.012 facies anterolateralis | A02.4.04.010 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.013 facies posterior | A02.4.04.010 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.014 sulcus nervi radialis | A02.4.04.013 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.015 margo medialis | A02.4.04.010 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.016 crista supraepicondylaris medialis | A02.4.04.015 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.017 processus supracondylaris | A02.4.04.015 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.018 margo lateralis | A02.4.04.010 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.019 crista supraepicondylaris lateralis | A02.4.04.018 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.020 tuberositas deltoidea | A02.4.04.010 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.021 condylus humeri | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.022 capitulum humeri | A02.4.04.021 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.023 trochlea humeri | A02.4.04.021 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.024 fossa olecrani | A02.4.04.021 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.025 fossa coronoidea | A02.4.04.021 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.026 fossa radialis | A02.4.04.021 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.027 epicondylus medialis | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |
| A02.4.04.028 sulcus nervi ulnaris | A02.4.04.027 | boundary/landmark | Missing mapping | Missing mapping |
| A02.4.04.029 epicondylus lateralis | A02.4.04.001 | bone region or surface feature | Missing mapping | Missing mapping |

A small proximal source interface has four open edges and two nonmanifold edges per side, with invalid vertex links. All original positions, material ownership and shoulder/elbow/attachment landmarks must constrain its repair. No bone pose, source shape or model has been changed.

Evidence: `reports/ta98-male/humerus-interface-repair-inventory.json`.

## Source defect and candidate repair

The original compound boundary is not usable for inside tests: four open edges, two nonmanifold edges and three invalid vertex links per side. Edge ownership shows the defect is entirely within the Bone material, not the cartilage junction as initially suspected. Two triangles (compound source faces 1985 and 1987) each have one vertex incident on only that face and lie on an edge already shared by two proper bone faces. They are dangling surface flaps, not closed tissue or a nutrient opening.

Removing only these two triangles yields one closed genus-zero humeral exterior per side, with no self-crossings, holes, nonmanifold edges, winding faults or invalid vertex links. Every original coordinate and retained face is exact; cartilage is unchanged. Removed surface area is approximately 0.792 mm² per side. Independent source/edge/topology replay passes. Nothing has been installed; attribute packing and actual rendering remain pending before a local model update.

This usable boundary localizes the long-head motor offshoot: 106 sampled controls/centroids/midpoints lie inside the bone, with maximum sampled depth about 0.516 mm. Distal muscle-target markers are outside bone. That supports preserving the target and junction while correcting the local transit, rather than moving the skeleton. Evidence: `humerus-dangling-face-repair-{plan,validation}.json` and `musculocutaneous-long-head-repaired-bone-transit.json`.
