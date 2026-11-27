# TA98 identities for the registered collection fascia, fat and connective parts

Opt-in review route: `?model=local-ta98&fascia=registered-collection-review` loads `fascia-fit-v13` (geometry of `fascia-fit-v12`; identities only). Not installed as default and not anatomically accepted.

Build (refuses to overwrite an existing stage; delete `fascia-fit-v13` first): `npm run build:registered-collection-fascia-identities`
Test: `npm run test:registered-collection-fascia` (skips when the local stages are absent).
Source: `scripts/bind-registered-collection-fascia-identities.py`.

## Status (v13)
All 1,840 fascia/fat/connective parts are bound to 365 distinct TA98 terms (exact where TA98 names the same structure, otherwise `parent`: the nearest broader TA98 term), plus 81 muscle, tendon and bursa concepts of the base catalogue. Nothing remains unbound. Geometry, chunks and bounds are asserted unchanged against v12 (script and test).

## Tiers, in order
1. `exact-normalized-name`: source name equals a TA98 English name or synonym after removing side and punctuation (kind `exact`).
2. `exact-name-ignoring-hyphenation`: same, ignoring hyphens and spacing (`sacroiliac` / `sacro-iliac`) (kind `exact`).
3. `reviewed-class-rule`: explicit regex to a TA98 code that was read from the authority. `exact` only where TA98 names the same structure; `parent` is a verified broader placement (an instance of a class, or an articular cartilage placed under the joint it lines).
4. `generic-class-ignoring-side-ordinal-digit-number`: ordinals, fingers/toes and vertebral levels removed, then a unique class match (kind `parent`).
Muscle concepts use `MUSCLE_SPEC` and unique authority FMA cross-references with an identical name.

Fuzzy matching is not used: it labelled a knee capsule with the kidney's fibrous capsule. Ambiguous or absent names are listed under `deferredParts` in `binding-ledger.json`.

## Deliberately unbound
- Articular cartilage of carpals, tarsals, metacarpals, metatarsals, phalanges, radius, ulna, femur shaft and ilium: TA98 has no cartilage term and the joint is not determined by the name alone.
- Facial fat compartments (nasolabial, jowl, temporal, cheek, orbital, subplatysmal): no TA98 terms. The system view files them under Integumentary › Subcutaneous tissue.
- Ligaments with no TA98 counterpart: cutaneous retaining ligaments of the face, tear trough, orbicularis retaining, ligamentous adhesions, brow/orbital thickenings.
- Ambiguous TA98 names (two codes, e.g. buccopharyngeal fascia, iliac fascia) and structures whose registered collection identity is unclear (trochlea, tarsal plates, superior transverse ligament).

## System view
`app/anatomy-hierarchy.ts` follows a bound TA98 chapter before its name rules. Fat pads therefore follow their organ (infrapatellar under the knee joint, buccal under the cheek); epicardial/paracardial fat go to the heart; unbound fat stays under Integumentary › Subcutaneous tissue.


## Update 2026-11-16: third reviewed rule set
+81 parts bound (1,398 parts, 309 distinct terms, 442 deferred): intermetacarpal and intermetatarsal dorsal/palmar/plantar/interosseous ligament classes (A03.5.11.402-404, A03.6.10.702-704), navicular-cuneiform and intercuneiform ligaments (A03.6.10.508/514/519), interosseous talocalcaneal (exact), cricothyroid ligament, lumbosacral, sacrococcygeal, xiphisternal and interchondral placements under their joints, pretracheal, prevertebral and superficial investing cervical fascia, temporal fascia, round ligament of the liver. Class placements are kind=parent. Geometry unchanged.

## Update 2026-11-17: nearest broader term for every remaining part
Policy chosen by the project owner: parts TA98 does not name are classified under the nearest broader TA98 term (kind=parent), never given the broader term's identity. Fifth reviewed rule set (56 codes, every one read from the authority): remaining joint cartilages under their joint group (wrist/hand joints, foot joints, elbow, distal radio-ulnar, ankle, hip/free-limb synovial joints), thumb and metacarpophalangeal ligaments, carpal/ulnar collateral and radio-ulnar ligaments, knee collaterals and anterolateral ligament, bursae, crural/gluteal/iliac/palmar/plantar/muscle fasciae, orbital and facial fasciae and septa under fascia of head, cutaneous (skin) ligaments, superficial fat pads under subcutaneous tissue (A16.0.03.001), deep fat pads under fascia of head or cervical fascia, epicardial/paracardial fat under serous pericardium/pericardium. Exact additions: superior/inferior tarsal plate (tarsus), cartilage of auditory tube, orbital trochlea, spinoglenoid ligament (= inferior transverse scapular ligament), plantar calcaneonavicular ligament. Ambiguous placements chosen: iliac fascia (iliopsoas fascia), buccopharyngeal fascia (fascia of head), lateral/medial collateral ligament (knee joint), femur/ilium cartilage (synovial joints of the free lower limb). The earlier "Deliberately unbound" list above is superseded.
