# Whole male model audit

Generated from the current local catalogue and every geometry chunk. No geometry changed.

Term correspondence comes from the existing four-source TA98 inventory (`reports/ta98-male/audit.json`); this audit joins it to newly measured geometry. It does not independently certify those anatomical matches.
Terminology reference: [FIPAT TA98](https://ifaa.unifr.ch/Public/EntryPage/ShowTA98.html). Catalogue, inventory and chunk SHA-256 hashes are saved for reproducibility.

**7,474 TA98 terms; 5,655 active meshes; 11,269,791 active triangles.**
Uncompressed geometry chunks: 368.6 MiB. 444 meshes hidden by default.

## TA98 correspondence

| Recorded status | Terms |
| --- | ---: |
| component-subset | 34 |
| non-surface-entity | 1,407 |
| containing-geometry-review-needed | 4,314 |
| existing-direct-correspondence | 1,324 |
| source-surface-review-pending | 227 |
| component-group | 36 |
| source-group-unresolved | 19 |
| no-corresponding-geometry-found | 8 |
| reconstructed-region-boundary-provisional | 1 |
| source-entity-type-review-needed | 3 |
| female-anatomy-outside-male-model | 101 |

## Coverage by TA98 chapter

| Chapter | Terms | Direct correspondence | Candidate surfaces | Containing geometry only | Other statuses |
| --- | ---: | ---: | ---: | ---: | ---: |
| General anatomy | 271 | 6 | 1 | 39 | 225 |
| Skeletal system | 981 | 65 | 0 | 388 | 528 |
| Articular system | 383 | 160 | 8 | 193 | 22 |
| Muscular system | 654 | 325 | 27 | 276 | 26 |
| Alimentary system | 490 | 56 | 18 | 288 | 128 |
| Respiratory system | 248 | 56 | 13 | 108 | 71 |
| Thoracic cavity | 28 | 1 | 0 | 13 | 14 |
| Urinary system | 102 | 6 | 0 | 79 | 17 |
| Genital systems | 300 | 9 | 2 | 126 | 163 |
| Abdominopelvic cavity | 93 | 4 | 1 | 49 | 39 |
| Endocrine glands | 34 | 8 | 2 | 18 | 6 |
| Cardiovascular system | 1300 | 314 | 118 | 816 | 52 |
| Lymphoid system | 202 | 68 | 1 | 119 | 14 |
| Nervous system | 1819 | 194 | 29 | 1442 | 154 |
| Sense organs | 506 | 48 | 5 | 310 | 143 |
| Integument | 63 | 4 | 2 | 50 | 7 |

## Terms without a source or containing-geometry match

This means the current concordance found no match. It does not prove the structure is absent from all source surfaces.

| TA98 | Term | Parent |
| --- | --- | --- |
| A02.4.10.006 | base of phalanx | None |
| A02.4.10.007 | shaft of phalanx | None |
| A02.4.10.008 | head of phalanx | None |
| A02.4.10.009 | trochlea of phalanx | A02.4.10.008 |
| A02.5.02.001 | pelvis | None |
| A02.5.02.003 | pubic arch | A02.5.02.001 |
| A02.5.02.005 | greater pelvis | A02.5.02.001 |
| A02.5.02.006 | lesser pelvis | A02.5.02.001 |

Every term, parent, selected mesh, source candidate and containing-term lead is recorded in `whole-model-audit.json`; the term CSV is suitable for manual review.

## Geometry observations

| Observation | Active meshes |
| --- | ---: |
| boundaryEdges | 4,083 |
| nonmanifoldEdges | 205 |
| inconsistentWindingEdges | 54 |
| connectedComponents | 1,410 |
| zeroAreaTriangles | 73 |
| collapsedTrianglesAt1Micron | 78 |
| duplicateTrianglesAt1Micron | 109 |
| Exact coincident geometry groups | 3 |

These counts identify review work. They do not imply that every open named surface should be closed.

## Coincident geometry to inspect

| Pieces | Sources | Hidden by default |
| --- | --- | --- |
| Nasal bone (left) · Suture-4 / Nasal bone (right) · Suture-4 | male-detail / male-detail | False / False |
| Zonular fibres (left) / Zonular fibres (right) | male-detail / male-detail | False / False |
| Caudate lobe branch of right hepatic artery / Caudate lobe branch of right hepatic artery | male / male | True / True |

## Review order

1. Resolve degenerate faces, repeated faces, nonmanifold edges and winding conflicts on whole solid objects after reviewing their anatomical identity.
2. Review hidden recovered structures in the male coordinate frame, with their adjacent organs and source provenance visible.
3. Inventory all descendants and adjoining interfaces for each organ before partitioning it. Represent openings as openings and keep unsupported layers or microstructures missing.
4. Inspect rendered surfaces and cross-sections, then validate self-intersections, tunnels and adjoining interfaces before accepting each organ.

## Meshes requiring topology review

| Mesh | Nonmanifold edges | Zero-area faces | Collapsed faces at 1 µm | Duplicate faces | Components |
| --- | ---: | ---: | ---: | ---: | ---: |
| Spinal dura | 5053 | 10100 | 10100 | 0 | 36 |
| Longus colli muscle (right) | 3391 | 0 | 0 | 2252 | 3 |
| Longus colli muscle (left) | 3386 | 0 | 0 | 2247 | 3 |
| Ascending pharyngeal artery (left) | 2580 | 0 | 0 | 1728 | 2 |
| Ascending pharyngeal artery (right) | 2580 | 0 | 0 | 1728 | 2 |
| Anterior cruciate ligament (left) | 1182 | 0 | 0 | 788 | 1 |
| Anterior cruciate ligament (right) | 1182 | 0 | 0 | 788 | 1 |
| Common palmar digital branches of ulnar nerve (left) | 852 | 0 | 0 | 576 | 2 |
| Common palmar digital branches of ulnar nerve (right) | 852 | 0 | 0 | 576 | 2 |
| Parieto-occipital sulcus (left) · Interlobar sulci | 313 | 0 | 0 | 225 | 2 |
| Parieto-occipital sulcus (right) · Interlobar sulci | 313 | 0 | 0 | 225 | 2 |
| Cavernous sinus (left) | 246 | 0 | 0 | 0 | 2 |
| Cavernous sinus (right) | 246 | 0 | 0 | 0 | 2 |
| Parieto-occipital sulcus (left) · Brain-Inner | 241 | 0 | 0 | 177 | 2 |
| Parieto-occipital sulcus (right) · Brain-Inner | 241 | 0 | 0 | 177 | 2 |
| Pancreatic duct | 172 | 0 | 168 | 0 | 35 |
| Falx cerebri | 136 | 0 | 0 | 1 | 1 |
| Middle cerebral artery (left) | 133 | 165 | 162 | 59 | 1 |
| Middle cerebral artery (right) | 120 | 151 | 148 | 49 | 1 |
| Vertebra T8 · Bone-5 | 114 | 0 | 0 | 82 | 1 |
| Palmal digital vei (right) | 84 | 90 | 90 | 3 | 25 |
| Palmal digital vei (left) | 84 | 90 | 90 | 3 | 25 |
| Left corona ciliaris | 64 | 7 | 50 | 26 | 44 |
| Cricopharyngeal ligament (left) | 64 | 0 | 0 | 0 | 1 |
| Cricopharyngeal ligament (right) | 64 | 0 | 0 | 0 | 1 |
| Deltoid fascia (left) | 63 | 0 | 0 | 0 | 1 |
| Deltoid fascia (right) | 63 | 0 | 0 | 0 | 1 |
| Radial nerve (inferior lateral brachial cutaneous n) (right) | 60 | 66 | 66 | 3 | 5 |
| Radial nerve (inferior lateral brachial cutaneous n) (left) | 60 | 66 | 66 | 3 | 5 |
| Occipital bone · Bone-6 | 59 | 0 | 38 | 4 | 3 |

## Limitations

- One-micron coordinate quantization joins texture/material seam vertices; very small features can collapse at this tolerance.
- Boundary edges and disconnected components are observations, not proof of anatomical defects. Named surfaces and shared tissue pieces can be intentionally open.
- Exact geometry hashes ignore vertex order, face order and winding at one-micron resolution; coincident alternatives require review before removal.
- A terminology match or ancestor mesh does not establish the shape, wall layers, microscopic structures or complete descendant anatomy.
- Whole-model self-intersections, inter-object penetration, tunnels, adjoining interfaces and visual cross-sections remain unverified; the separately reviewed gastric tract has its own validation report.
- No geometry was cut, moved or reconstructed. Local sources remain local; redistribution rights have not been established.

## Reproduce

Run `scripts/audit-ta98-male-model.py` with Python, NumPy and SciPy installed. In the current workspace:

```sh
/Users/christiantzurcanu/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 scripts/audit-ta98-male-model.py
```

`--reuse-geometry` refreshes the reports only after verifying the saved catalogue and every geometry-chunk hash.
