# Drugs hierarchy

Drugs appears before Physical Exercise and Chakras in the Additional hierarchies menu. It loads only when selected.

## Imported snapshot

Source: ChEMBL_37, released 2026-05-01. All ATC entries supplied by that ChEMBL release are retained; this is a pinned classification snapshot, not a live WHO index.

| Data | Count |
| --- | ---: |
| ATC level 1 | 14 |
| ATC level 2 | 90 |
| ATC level 3 | 248 |
| ATC level 4 | 841 |
| ATC level 5 | 5,579 |
| Drug records | 5,467 |
| Drug records without ATC codes | 1,816 |
| Curated mechanism records | 4,106 |
| Human target genes | 1,039 |
| Target genes with anatomy mappings | 205 |

## Categories

| ATC group | Open in viewer |
| --- | --- |
| A | [Alimentary tract and metabolism](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AA&context=0.18) |
| B | [Blood and blood forming organs](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AB&context=0.18) |
| C | [Cardiovascular system](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AC&context=0.18) |
| D | [Dermatologicals](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AD&context=0.18) |
| G | [Genito urinary system and sex hormones](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AG&context=0.18) |
| H | [Systemic hormonal preparations, excl. sex hormones and insulins](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AH&context=0.18) |
| J | [Antiinfectives for systemic use](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AJ&context=0.18) |
| L | [Antineoplastic and immunomodulating agents](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AL&context=0.18) |
| M | [Musculo-skeletal system](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AM&context=0.18) |
| N | [Nervous system](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AN&context=0.18) |
| P | [Antiparasitic products, insecticides and repellents](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AP&context=0.18) |
| R | [Respiratory system](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AR&context=0.18) |
| S | [Sensory organs](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AS&context=0.18) |
| V | [Various](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=ATC%3AV&context=0.18) |
| UNCLASSIFIED | [Unclassified](http://localhost:3016/?model=male-detail&tree=guest%3Adrugs&bio=UNCLASSIFIED&context=0.18) |

## Mapping rules

- The five ATC levels lead to the associated ChEMBL drug records, their curated mechanisms, human target genes, and available anatomical items or sets. Multiple ATC assignments and shared targets retain all their parent branches.
- Drug records include substances with ATC assignments and records classified as phase 4 in ChEMBL. Historical or withdrawn records remain present; this is not a list of currently authorized products. Salts and combination records retain their source identities.
- A gene is linked only for a molecular mechanism with an explicitly human target. Components use source HGNC cross references, or an unambiguous UniProt-to-HGNC match from the pinned HGNC snapshot. Nonhuman pathogen targets are not assigned to human homologues.
- Anatomy comes from the existing gene-location evidence. These are mapped target-gene locations; they do not establish drug action, measured target expression, or clinical effects throughout an entire organ mesh. Cell Types links refer to existing marker associations.
- Unmapped categories, drugs and mechanisms remain searchable without anatomy checkboxes. Missing mappings are not filled with whole systems or guessed organs.
- Share URLs, slides, connection snapshots and embeds use the same built-in hierarchy ID, `guest:drugs`. Search accepts drug names and synonyms, ChEMBL IDs, ATC codes, target names and HGNC IDs.

## Sources

All source credits are in [ATTRIBUTION.md](ATTRIBUTION.md#drugs). Exact request URLs and source hashes are in [scripts/drug-sources.json](scripts/drug-sources.json). Mechanism evidence, original references, mapping gaps and input hashes are in [public/assets/drug-provenance.json](public/assets/drug-provenance.json).

## Rebuild and validate

```sh
npm run build:drugs -- --download  # fetch missing pinned source snapshots
npm run build:drugs               # rebuild from the hash-verified cache
npm run build:drugs -- --check
npm run test:drugs
npm run check
npm run build
```

The ignored `work/drugs` directory stores source snapshots. UniProt matching also uses the pinned HGNC source in `work/biology`, downloaded by `npm run build:biology -- --download`. Deployable data is in `public/assets/drugs.json`.
