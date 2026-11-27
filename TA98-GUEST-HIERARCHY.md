# TA98 guest hierarchy

Choose **Layers → Additional hierarchies → TA98**. The guest file is
`public/assets/ta98.json`, in the existing `human-atlas-hierarchy/v2` format.
Preview: <http://localhost:3016/?model=male-detail&tree=guest:ta98>.
Female preview: <http://localhost:3016/?model=female&tree=guest:ta98>.

All 7,474 entries from the TA98 SQLite extract are retained under 16 chapter
navigation containers. Names, Latin synonyms, entity keys, and published TA98
codes are searchable. Original partonomic parents are preserved, including
cross-chapter relationships. Chapter entry points refer to the same canonical
nodes rather than duplicating anatomical identities. The three hand-phalanx
parent gaps use the same published TA98 correction as the existing concordance.

Expand a term to see its code, Latin name, entity type, and model coverage.
Click a mapped term to select its available parts and mapped descendants.
Its checkbox hides or restores that set. Unmodeled terms remain searchable
and expandable to their coverage notes, with no anatomy checkbox.

Mappings include public detailed/complete male, standard male, female, local
male/female and local reference catalogues. The complete male view shares the
detailed male catalogue. Local IDs resolve only when that local model is loaded.
Exact source name/synonym matches supplement the existing reviewed concordance,
particularly for aggregate organ concepts composed of several meshes.
Ambiguous names and fuzzy similarity are not used as evidence.

The [coverage report](reports/ta98/README.md) contains every term with linked
atlas IDs. [coverage.json](reports/ta98/coverage.json) includes exact matches,
parent placements, match methods, complete descendant selection concept sets,
source entity types and unmatched atlas concepts.
[coverage.csv](reports/ta98/coverage.csv) is the spreadsheet version.

Parent placement means a modeled subdivision belongs under a term; it does not
give that subdivision the parent's anatomical identity. A selection made from
descendants is available coverage, not proof that the complete structure exists.
Annotation/attachment meshes and suppressed parts are excluded.

Missing mappings require source and geometry review before Blender work.
TA98 includes directional terms, regions, cavities, openings and other
immaterial/non-physical entries that cannot all be represented as standalone
solid meshes. This guest hierarchy does not fabricate any geometry.

## Rebuild and validate

Run `npm run build:ta98` followed by `npm run test:ta98`. The builder uses
`/tmp/human-atlas-ta98.sqlite`, downloading the existing project's source extract
when absent. To use another verified copy, run
`python3 scripts/build-ta98-hierarchy.py /path/to/ta98.sqlite`.
Regenerate after changes to atlas catalogues or the reviewed concordance.

The full file exceeds the viewer's 500 KB arbitrary remote guest limit, so it is
registered as a trusted built-in guest with the existing 20 MB limit. The same
origin `/assets/ta98.json` URL can also be used in Load from URL. Deployment of
these changes is needed before the GitHub Pages viewer has this hierarchy.

Sources: [IFAA TA98 navigation](https://ifaa.unifr.ch/Public/EntryPage/ShowTA98.html),
[TA98 SQLite extract](https://github.com/mhalle/ta98-sqlite),
and the existing [mapping audit](TA98-HIERARCHY-AUDIT.md).


## Standard browser catalogue identity correction (2026-10-03)

The Systems, Regions and Depth views now read explicit model `ta98Term`/`ta98Kind` assignments for modeled leaves and reconstruct their lineage from the existing reviewed direct-parent index. Previously these views relied on shared source-name metadata alone; explicit stomach regions and many constructed regional placements could lose ancestry or appear in the wrong system. Stomach cardia, fundus, body and pyloric part now retain stomach ancestry in both main models. Parent placements include their broader containment term but have no false exact TA98/THA/Latin badge.

`audit-ta98-runtime-hierarchy-identities.mjs` independently checks every available piece, every explicit term against the authority code index and every lineage against direct parent edges. It binds source/catalogue/metadata hashes in `reports/ta98-completeness/runtime-hierarchy-identities.json`. Compared with the prior shared-metadata fallback, 258 male/95 female leaves have changed semantic identity or lineage; 21 male leaves change system placement. Side-specific source FMA and source UBERON identifiers are preserved when their exact source mapping agrees with the explicit identity. Unknown/conflicting explicit terms cannot resurrect stale source identifiers.

The taxonomy validator now runs full Systems/Regions/Depth one-parent and selection-membership checks on both runtime catalogues, alongside contradictory-identity, unknown-term, parent-placement and stomach tests. Related modeled-parent resolution continues to use exact catalogue identities. No meshes, textures, contacts or accepted anatomical representations are added by this browser correction. Local donor registration and full anatomical acceptance remain separate obligations.
